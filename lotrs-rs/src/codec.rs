//! Canonical serialization — port of `lotrs-py/codec.py`.
//!
//! Fixed-width packing for uniform / bit-dropped components;
//! Golomb–Rice coding for Gaussian-distributed components.  All decoders
//! enforce canonical form: out-of-range coefficients, nonzero padding
//! bits, and trailing bytes produce a [`CodecError`].
//!
//! The byte layout is identical to the Python reference, so anything
//! encoded by `lotrs-py/codec.py` round-trips bit-for-bit through Rust.
//!
//! **All public decoders are total: they return `Result<_, CodecError>`
//! and never panic on malformed input.** That property is crucial for
//! the scheme's `verify()` entry point, which surfaces *any* decoder
//! error as a plain `false` (see `lotrs::verify` once wired up).

use crate::params::LoTRSParams;
use crate::ring::Ring;

/// Codec-level error.  Carries no details about the underlying
/// deserialization failure — the caller is expected to surface a
/// binary accept / reject signal only.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CodecError {
    LengthMismatch,
    OutOfRange,
    NonZeroPadding,
    TrailingBytes,
    Truncated,
    UnaryOverflow,
    NonCanonical,
}

impl core::fmt::Display for CodecError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            CodecError::LengthMismatch => "length mismatch",
            CodecError::OutOfRange => "value out of range",
            CodecError::NonZeroPadding => "non-zero padding bits",
            CodecError::TrailingBytes => "trailing bytes",
            CodecError::Truncated => "truncated data",
            CodecError::UnaryOverflow => "unary code overflow",
            CodecError::NonCanonical => "non-canonical encoding",
        })
    }
}

pub type Result<T> = core::result::Result<T, CodecError>;

// =========================================================================
//  Signature structure
// =========================================================================

/// Non-interactive proof component of a LoTRS signature.
#[derive(Clone, Debug)]
pub struct Proof {
    /// `B_bin^(1)` — high-bit part of the binary-proof commitment (R_qhat).
    pub b_bin_hi: Vec<Vec<u64>>,
    /// `w̃^(1)` — high-bit slices of the DualMS commitment, indices `j >= 1`.
    /// Length is `kappa - 1`; each entry is a `Vec<k>` of R_q polynomials.
    /// For the only supported case `kappa == 1` this is empty.
    pub w_tilde_hi: Vec<Vec<Vec<u64>>>,
    /// Challenge polynomial `x ∈ C`.
    pub x: Vec<u64>,
    /// `f_1` — flattened `kappa * (beta - 1)` polynomials in R_q.
    pub f1: Vec<Vec<u64>>,
    /// `z_b` — `n_hat + k_hat` polynomials in R_qhat.
    pub z_b: Vec<Vec<u64>>,
}

/// Complete aggregated LoTRS signature.
#[derive(Clone, Debug)]
pub struct Signature {
    pub pi: Proof,
    /// `z̃` — `l` polynomials in R_q.
    pub z_tilde: Vec<Vec<u64>>,
    /// `r̃` — `l'` polynomials in R_q.
    pub r_tilde: Vec<Vec<u64>>,
    /// `ẽ` — `k` polynomials in R_q.
    pub e_tilde: Vec<Vec<u64>>,
}

// =========================================================================
//  Bit stream helpers
// =========================================================================

pub struct BitWriter {
    /// Packed little-endian within each byte; LSB first.
    bytes: Vec<u8>,
    bit_len: usize,
}

impl BitWriter {
    pub fn new() -> Self {
        Self {
            bytes: Vec::new(),
            bit_len: 0,
        }
    }

    pub fn write_bits(&mut self, value: u64, n: u32) {
        for i in 0..n {
            let bit = ((value >> i) & 1) as u8;
            let byte_idx = self.bit_len >> 3;
            let bit_idx = (self.bit_len & 7) as u8;
            if byte_idx >= self.bytes.len() {
                self.bytes.push(0);
            }
            self.bytes[byte_idx] |= bit << bit_idx;
            self.bit_len += 1;
        }
    }

    pub fn write_unary(&mut self, value: u64) {
        for _ in 0..value {
            self.write_bits(1, 1);
        }
        self.write_bits(0, 1);
    }

    pub fn pad_to_byte(&mut self) {
        while self.bit_len & 7 != 0 {
            self.write_bits(0, 1);
        }
    }

    pub fn into_bytes(mut self) -> Vec<u8> {
        self.pad_to_byte();
        self.bytes
    }

    pub fn bit_len(&self) -> usize {
        self.bit_len
    }
}

pub struct BitReader<'a> {
    data: &'a [u8],
    pos: usize, // bit position
}

impl<'a> BitReader<'a> {
    pub fn new(data: &'a [u8]) -> Self {
        Self { data, pos: 0 }
    }

    pub fn read_bits(&mut self, n: u32) -> Result<u64> {
        let mut v = 0u64;
        for i in 0..n {
            let byte_idx = self.pos >> 3;
            let bit_idx = (self.pos & 7) as u8;
            if byte_idx >= self.data.len() {
                return Err(CodecError::Truncated);
            }
            let bit = ((self.data[byte_idx] >> bit_idx) & 1) as u64;
            v |= bit << i;
            self.pos += 1;
        }
        Ok(v)
    }

    /// Read a unary-coded value; reject if it exceeds `max_val` (DoS guard).
    pub fn read_unary(&mut self, max_val: u64) -> Result<u64> {
        let mut count = 0u64;
        loop {
            let b = self.read_bits(1)?;
            if b == 0 {
                return Ok(count);
            }
            count += 1;
            if count > max_val {
                return Err(CodecError::UnaryOverflow);
            }
        }
    }

    /// Advance to the next byte boundary; reject nonzero padding.
    pub fn check_padding(&mut self) -> Result<()> {
        while self.pos & 7 != 0 {
            if self.read_bits(1)? != 0 {
                return Err(CodecError::NonZeroPadding);
            }
        }
        Ok(())
    }

    pub fn check_exhausted(&mut self) -> Result<()> {
        self.check_padding()?;
        if (self.pos >> 3) != self.data.len() {
            return Err(CodecError::TrailingBytes);
        }
        Ok(())
    }

    pub fn consumed_bytes(&self) -> usize {
        (self.pos + 7) >> 3
    }
}

// =========================================================================
//  Fixed-width polynomial packing
// =========================================================================

/// Pack `d` integer coefficients at `dx` bits each, shifting by `+offset`
/// so negative values become non-negative.  Byte-aligned.
pub fn pack_fixed(coeffs: &[i64], dx: u32, offset: i64) -> Result<Vec<u8>> {
    let mut w = BitWriter::new();
    let limit = 1i64 << dx;
    for &c in coeffs {
        let v = c + offset;
        if v < 0 || v >= limit {
            return Err(CodecError::OutOfRange);
        }
        w.write_bits(v as u64, dx);
    }
    w.pad_to_byte();
    Ok(w.into_bytes())
}

/// Unpack `d` coefficients at `dx` bits each.  Rejects nonzero padding.
/// Returned values are signed (offset removed).
pub fn unpack_fixed(data: &[u8], d: usize, dx: u32, offset: i64) -> Result<(Vec<i64>, usize)> {
    let mut r = BitReader::new(data);
    let mut coeffs = Vec::with_capacity(d);
    for _ in 0..d {
        let v = r.read_bits(dx)? as i64;
        coeffs.push(v - offset);
    }
    r.check_padding()?;
    Ok((coeffs, r.consumed_bytes()))
}

// =========================================================================
//  Golomb-Rice polynomial packing
// =========================================================================

/// Optimal Rice parameter for a discrete Gaussian with standard deviation
/// `sigma`.  Matches Python's `optimal_rice_k`.
pub fn optimal_rice_k(sigma: f64) -> u32 {
    if sigma < 1.0 {
        return 0;
    }
    let v = (1.1774_f64 * sigma).log2().floor() as i64;
    v.max(0) as u32
}

/// Golomb–Rice encode signed coefficients.  Byte-aligned after all `d`.
pub fn pack_rice(coeffs: &[i64], rice_k: u32, bound: i64) -> Result<Vec<u8>> {
    let mut w = BitWriter::new();
    let low_mask: u64 = if rice_k == 0 { 0 } else { (1u64 << rice_k) - 1 };
    for &c in coeffs {
        let abs_c = c.unsigned_abs();
        if abs_c as i64 > bound {
            return Err(CodecError::OutOfRange);
        }
        let low = abs_c as u64 & low_mask;
        let high = abs_c as u64 >> rice_k;
        w.write_bits(low, rice_k);
        w.write_unary(high);
        if abs_c != 0 {
            w.write_bits(if c < 0 { 1 } else { 0 }, 1);
        }
    }
    w.pad_to_byte();
    Ok(w.into_bytes())
}

/// Golomb–Rice decode `d` signed coefficients.  Rejects out-of-range,
/// nonzero padding, and pathological unary runs.
pub fn unpack_rice(data: &[u8], d: usize, rice_k: u32, bound: i64) -> Result<(Vec<i64>, usize)> {
    let max_high = (bound >> rice_k) as u64 + 1;
    let mut r = BitReader::new(data);
    let mut coeffs = Vec::with_capacity(d);
    for _ in 0..d {
        let low = r.read_bits(rice_k)?;
        let high = r.read_unary(max_high)?;
        let abs_c = (high << rice_k) | low;
        if abs_c as i64 > bound {
            return Err(CodecError::OutOfRange);
        }
        if abs_c == 0 {
            coeffs.push(0);
        } else {
            let sign = r.read_bits(1)?;
            let v = abs_c as i64;
            coeffs.push(if sign == 1 { -v } else { v });
        }
    }
    r.check_padding()?;
    Ok((coeffs, r.consumed_bytes()))
}

// =========================================================================
//  Challenge encoding
// =========================================================================

fn bits_for_positions(d: usize) -> u32 {
    core::cmp::max(1, 64 - (d as u64 - 1).leading_zeros())
}

/// Encode a challenge polynomial (w nonzero ±1 coefficients).
/// Layout: `w` positions (sorted ascending, `ceil(log2 d)` bits each)
///   + `w` sign bits (0 = +1, 1 = −1), byte-aligned.
///
/// Returns [`CodecError::OutOfRange`] if any coefficient is not in
/// `{-1, 0, 1}` or if the number of non-zero coefficients is not `w`.
/// Without this check the encoder would silently project every
/// non-zero integer to `±1`, so `unpack_challenge(pack_challenge(p))`
/// would disagree with `p` for malformed input.
pub fn pack_challenge(poly: &[i64], w: usize, d: usize) -> Result<Vec<u8>> {
    if poly.len() != d {
        return Err(CodecError::LengthMismatch);
    }
    let pos_bits = bits_for_positions(d);
    let mut positions: Vec<usize> = Vec::with_capacity(w);
    for i in 0..d {
        match poly[i] {
            0 => {}
            1 | -1 => positions.push(i),
            _ => return Err(CodecError::OutOfRange),
        }
    }
    if positions.len() != w {
        return Err(CodecError::OutOfRange);
    }
    let mut wr = BitWriter::new();
    for &p in &positions {
        wr.write_bits(p as u64, pos_bits);
    }
    for &p in &positions {
        wr.write_bits(if poly[p] < 0 { 1 } else { 0 }, 1);
    }
    wr.pad_to_byte();
    Ok(wr.into_bytes())
}

/// Decode a challenge polynomial.  Rejects duplicate / out-of-order
/// positions and non-binary sign bits.
pub fn unpack_challenge(data: &[u8], w: usize, d: usize) -> Result<(Vec<i64>, usize)> {
    let pos_bits = bits_for_positions(d);
    let mut r = BitReader::new(data);

    let mut positions = Vec::with_capacity(w);
    for _ in 0..w {
        let p = r.read_bits(pos_bits)? as usize;
        if p >= d {
            return Err(CodecError::OutOfRange);
        }
        if let Some(&last) = positions.last() {
            if p <= last {
                return Err(CodecError::NonCanonical);
            }
        }
        positions.push(p);
    }
    let mut signs = Vec::with_capacity(w);
    for _ in 0..w {
        signs.push(r.read_bits(1)? as u8);
    }
    r.check_padding()?;

    let mut poly = vec![0i64; d];
    for (&p, &s) in positions.iter().zip(signs.iter()) {
        poly[p] = if s == 1 { -1 } else { 1 };
    }
    Ok((poly, r.consumed_bytes()))
}

// =========================================================================
//  High-level codec wrapping a parameter set
// =========================================================================

/// Canonical encoder / decoder for LoTRS objects.  All widths and Rice
/// parameters are derived from the attached [`LoTRSParams`].
#[derive(Clone)]
pub struct LoTRSCodec {
    pub par: LoTRSParams,
    pub dx_pk: u32,
    pub dx_bbin: u32,
    pub dx_whi: u32,
    pub rice_f1: u32,
    pub bound_f1: i64,
    pub rice_zb: u32,
    pub bound_zb: i64,
    pub rice_zt: u32,
    pub bound_zt: i64,
    pub rice_rt: u32,
    pub bound_rt: i64,
    pub rice_et: u32,
    pub bound_et: i64,
}

impl LoTRSCodec {
    pub fn new(par: LoTRSParams) -> Self {
        let dx_pk = 64 - (par.q - 1).leading_zeros();
        let dx_bbin = (64 - (par.q_hat - 1).leading_zeros()) - par.K_B;
        let dx_whi = dx_pk - par.K_w;

        let rice_f1 = optimal_rice_k(par.sigma_a());
        let bound_f1 = (6.0 * par.sigma_a()).ceil() as i64;

        let rice_zb = optimal_rice_k(par.sigma_b());
        let bound_zb = (6.0 * par.sigma_b()).ceil() as i64;

        // e_tilde has effective width sigma_0 + sigma_0_prime
        // (triangle inequality); see the matching note in lotrs.rs verify.
        let sigma_z = par.sigma_0() * (par.T as f64).sqrt();
        let sigma_r = par.sigma_0_prime() * (par.T as f64).sqrt();
        let sigma_e = (par.sigma_0() + par.sigma_0_prime()) * (par.T as f64).sqrt();

        let t = par.tail_t;
        let rice_zt = optimal_rice_k(sigma_z);
        let bound_zt = (t * par.sigma_0() * ((par.T * par.d * par.l) as f64).sqrt()).ceil() as i64;

        let rice_rt = optimal_rice_k(sigma_r);
        let bound_rt =
            (t * par.sigma_0_prime() * ((par.T * par.d * par.l_prime) as f64).sqrt()).ceil() as i64;

        let rice_et = optimal_rice_k(sigma_e);
        let bound_et =
            (t * (par.sigma_0() + par.sigma_0_prime()) * ((par.T * par.d * par.k) as f64).sqrt())
                .ceil() as i64;

        Self {
            par,
            dx_pk,
            dx_bbin,
            dx_whi,
            rice_f1,
            bound_f1,
            rice_zb,
            bound_zb,
            rice_zt,
            bound_zt,
            rice_rt,
            bound_rt,
            rice_et,
            bound_et,
        }
    }

    // ---- public parameters ---------------------------------------------

    pub fn pp_encode(&self, pp: &[u8]) -> Result<Vec<u8>> {
        if pp.len() != 32 {
            return Err(CodecError::LengthMismatch);
        }
        Ok(pp.to_vec())
    }
    pub fn pp_decode(&self, data: &[u8]) -> Result<Vec<u8>> {
        if data.len() != 32 {
            return Err(CodecError::LengthMismatch);
        }
        Ok(data.to_vec())
    }

    // ---- secret key (seed form) ---------------------------------------

    pub fn sk_encode(&self, seed: &[u8]) -> Result<Vec<u8>> {
        if seed.len() != 32 {
            return Err(CodecError::LengthMismatch);
        }
        Ok(seed.to_vec())
    }
    pub fn sk_decode(&self, data: &[u8]) -> Result<Vec<u8>> {
        if data.len() != 32 {
            return Err(CodecError::LengthMismatch);
        }
        Ok(data.to_vec())
    }

    // ---- public key ----------------------------------------------------

    /// Encode a public key — `k` polynomials in `R_q`, each `d` coeffs
    /// packed at `dx_pk` bits without offset (unsigned canonical form).
    pub fn pk_encode(&self, pk: &[Vec<u64>]) -> Result<Vec<u8>> {
        let par = &self.par;
        if pk.len() != par.k {
            return Err(CodecError::LengthMismatch);
        }
        let mut out = Vec::new();
        for poly in pk {
            if poly.len() != par.d {
                return Err(CodecError::LengthMismatch);
            }
            // reinterpret u64 → i64 with offset = 0
            let signed: Vec<i64> = poly
                .iter()
                .map(|&c| {
                    debug_assert!(c < par.q);
                    c as i64
                })
                .collect();
            let chunk = pack_fixed(&signed, self.dx_pk, 0)?;
            out.extend_from_slice(&chunk);
        }
        Ok(out)
    }

    // ---- signature -----------------------------------------------------

    /// Encode a full aggregated LoTRS signature into a canonical byte stream.
    ///
    /// Matches `lotrs-py/codec.py::LoTRSCodec.sig_encode` byte-for-byte.
    /// Returns [`CodecError::LengthMismatch`] if any component vector or
    /// any individual polynomial has the wrong length for the bound
    /// parameter set — without this check the encoder would produce
    /// non-canonical bytes for malformed inputs.
    pub fn sig_encode(&self, sig: &Signature) -> Result<Vec<u8>> {
        let par = &self.par;
        let r_q = Ring::new(par.q, par.d);
        let r_qhat = Ring::new(par.q_hat, par.d);

        // helper — rejects any polynomial whose coefficient length isn't d
        let check_len = |p: &Vec<u64>| -> Result<()> {
            if p.len() != par.d {
                Err(CodecError::LengthMismatch)
            } else {
                Ok(())
            }
        };

        let mut out = Vec::new();

        // B_bin^(1) -- fixed width over R_qhat, centred offset
        if sig.pi.b_bin_hi.len() != par.n_hat {
            return Err(CodecError::LengthMismatch);
        }
        let offset_bbin = 1i64 << (self.dx_bbin - 1);
        for poly in &sig.pi.b_bin_hi {
            check_len(poly)?;
            let centred = r_qhat.centered(poly);
            out.extend_from_slice(&pack_fixed(&centred, self.dx_bbin, offset_bbin)?);
        }

        // w_tilde^(1) -- (kappa - 1) vectors of k ring elements over R_q
        if sig.pi.w_tilde_hi.len() != par.kappa - 1 {
            return Err(CodecError::LengthMismatch);
        }
        let offset_whi = 1i64 << (self.dx_whi - 1);
        for j_vec in &sig.pi.w_tilde_hi {
            if j_vec.len() != par.k {
                return Err(CodecError::LengthMismatch);
            }
            for poly in j_vec {
                check_len(poly)?;
                let centred = r_q.centered(poly);
                out.extend_from_slice(&pack_fixed(&centred, self.dx_whi, offset_whi)?);
            }
        }

        // x — challenge
        check_len(&sig.pi.x)?;
        let x_centred = r_q.centered(&sig.pi.x);
        out.extend_from_slice(&pack_challenge(&x_centred, par.w, par.d)?);

        // f1 — Rice, R_q
        if sig.pi.f1.len() != par.kappa * (par.beta - 1) {
            return Err(CodecError::LengthMismatch);
        }
        for poly in &sig.pi.f1 {
            check_len(poly)?;
            out.extend_from_slice(&pack_rice(
                &r_q.centered(poly),
                self.rice_f1,
                self.bound_f1,
            )?);
        }

        // z_b — Rice, R_qhat
        if sig.pi.z_b.len() != par.n_hat + par.k_hat {
            return Err(CodecError::LengthMismatch);
        }
        for poly in &sig.pi.z_b {
            check_len(poly)?;
            out.extend_from_slice(&pack_rice(
                &r_qhat.centered(poly),
                self.rice_zb,
                self.bound_zb,
            )?);
        }

        // z_tilde — Rice, R_q
        if sig.z_tilde.len() != par.l {
            return Err(CodecError::LengthMismatch);
        }
        for poly in &sig.z_tilde {
            check_len(poly)?;
            out.extend_from_slice(&pack_rice(
                &r_q.centered(poly),
                self.rice_zt,
                self.bound_zt,
            )?);
        }

        // r_tilde — Rice, R_q
        if sig.r_tilde.len() != par.l_prime {
            return Err(CodecError::LengthMismatch);
        }
        for poly in &sig.r_tilde {
            check_len(poly)?;
            out.extend_from_slice(&pack_rice(
                &r_q.centered(poly),
                self.rice_rt,
                self.bound_rt,
            )?);
        }

        // e_tilde — Rice, R_q
        if sig.e_tilde.len() != par.k {
            return Err(CodecError::LengthMismatch);
        }
        for poly in &sig.e_tilde {
            check_len(poly)?;
            out.extend_from_slice(&pack_rice(
                &r_q.centered(poly),
                self.rice_et,
                self.bound_et,
            )?);
        }

        Ok(out)
    }

    /// Decode and validate a signature.  **Total function** — returns
    /// [`CodecError`] for every malformed-input class (length mismatch,
    /// out-of-range coefficient, nonzero padding bit, trailing bytes,
    /// oversized Rice unary run, …) without panicking.
    pub fn sig_decode(&self, data: &[u8]) -> Result<Signature> {
        let par = &self.par;
        let mut pos = 0usize;

        // helper — consume one fixed-width polynomial with the centred offset
        let fixed_poly_bytes = |dx: u32| ((dx as usize) * par.d + 7) / 8;

        // helper — consume a fixed chunk and advance pos
        let consume_fixed =
            |data: &[u8], pos: &mut usize, n: usize, dx: u32, q: u64| -> Result<Vec<Vec<u64>>> {
                let poly_bytes = fixed_poly_bytes(dx);
                let offset = 1i64 << (dx - 1);
                let mut out = Vec::with_capacity(n);
                for _ in 0..n {
                    if *pos + poly_bytes > data.len() {
                        return Err(CodecError::Truncated);
                    }
                    let chunk = &data[*pos..*pos + poly_bytes];
                    let (coeffs, _) = unpack_fixed(chunk, par.d, dx, offset)?;
                    // reduce into [0, q) — centred-form decoding may yield
                    // values outside [0, q) only if dx is wider than needed;
                    // we rem_euclid them into canonical form here.
                    let mut poly = Vec::with_capacity(par.d);
                    for c in coeffs {
                        let v = c.rem_euclid(q as i64) as u64;
                        poly.push(v);
                    }
                    out.push(poly);
                    *pos += poly_bytes;
                }
                Ok(out)
            };

        // B_bin^(1)
        let b_bin_hi = consume_fixed(data, &mut pos, par.n_hat, self.dx_bbin, par.q_hat)?;

        // w_tilde^(1) — (kappa - 1) vectors of k polys
        let n_whi = (par.kappa - 1) * par.k;
        let w_hi_flat = consume_fixed(data, &mut pos, n_whi, self.dx_whi, par.q)?;
        let mut w_tilde_hi: Vec<Vec<Vec<u64>>> = Vec::with_capacity(par.kappa - 1);
        let mut idx = 0usize;
        for _ in 0..(par.kappa - 1) {
            w_tilde_hi.push(w_hi_flat[idx..idx + par.k].to_vec());
            idx += par.k;
        }

        // x — challenge
        if pos >= data.len() {
            return Err(CodecError::Truncated);
        }
        let (x_centred, consumed) = unpack_challenge(&data[pos..], par.w, par.d)?;
        let r_q = Ring::new(par.q, par.d);
        let x = r_q.from_centered(&x_centred);
        pos += consumed;

        // helper — consume Rice-coded polys
        let consume_rice = |data: &[u8],
                            pos: &mut usize,
                            n: usize,
                            rice_k: u32,
                            bound: i64,
                            ring: &Ring|
         -> Result<Vec<Vec<u64>>> {
            let mut out = Vec::with_capacity(n);
            for _ in 0..n {
                if *pos > data.len() {
                    return Err(CodecError::Truncated);
                }
                let (coeffs, consumed) = unpack_rice(&data[*pos..], par.d, rice_k, bound)?;
                out.push(ring.from_centered(&coeffs));
                *pos += consumed;
            }
            Ok(out)
        };

        let r_qhat = Ring::new(par.q_hat, par.d);

        let f1 = consume_rice(
            data,
            &mut pos,
            par.kappa * (par.beta - 1),
            self.rice_f1,
            self.bound_f1,
            &r_q,
        )?;
        let z_b = consume_rice(
            data,
            &mut pos,
            par.n_hat + par.k_hat,
            self.rice_zb,
            self.bound_zb,
            &r_qhat,
        )?;
        let z_tilde = consume_rice(data, &mut pos, par.l, self.rice_zt, self.bound_zt, &r_q)?;
        let r_tilde = consume_rice(
            data,
            &mut pos,
            par.l_prime,
            self.rice_rt,
            self.bound_rt,
            &r_q,
        )?;
        let e_tilde = consume_rice(data, &mut pos, par.k, self.rice_et, self.bound_et, &r_q)?;

        if pos != data.len() {
            return Err(CodecError::TrailingBytes);
        }

        Ok(Signature {
            pi: Proof {
                b_bin_hi,
                w_tilde_hi,
                x,
                f1,
                z_b,
            },
            z_tilde,
            r_tilde,
            e_tilde,
        })
    }

    // ---- public key ----------------------------------------------------

    /// Decode and validate a public key.  Returns `CodecError` for every
    /// malformed encoding class.
    pub fn pk_decode(&self, data: &[u8]) -> Result<Vec<Vec<u64>>> {
        let par = &self.par;
        let poly_bytes = ((self.dx_pk as usize * par.d) + 7) / 8;
        let expected = par.k * poly_bytes;
        if data.len() != expected {
            return Err(CodecError::LengthMismatch);
        }
        let mut pk = Vec::with_capacity(par.k);
        for i in 0..par.k {
            let chunk = &data[i * poly_bytes..(i + 1) * poly_bytes];
            let (coeffs, _) = unpack_fixed(chunk, par.d, self.dx_pk, 0)?;
            let mut poly = Vec::with_capacity(par.d);
            for c in coeffs {
                if c < 0 || (c as u64) >= par.q {
                    return Err(CodecError::OutOfRange);
                }
                poly.push(c as u64);
            }
            pk.push(poly);
        }
        Ok(pk)
    }
}

// -------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::TEST_PARAMS;

    #[test]
    fn fixed_roundtrip() {
        let xs = vec![-5i64, 3, 0, 7, -1, 2];
        let enc = pack_fixed(&xs, 4, 8).unwrap();
        let (dec, _) = unpack_fixed(&enc, xs.len(), 4, 8).unwrap();
        assert_eq!(dec, xs);
    }

    #[test]
    fn rice_roundtrip_zeros_and_small() {
        let xs = vec![0i64, 1, -1, 2, -3, 5, 0];
        let bound = 10;
        let enc = pack_rice(&xs, 2, bound).unwrap();
        let (dec, _) = unpack_rice(&enc, xs.len(), 2, bound).unwrap();
        assert_eq!(dec, xs);
    }

    #[test]
    fn rice_rejects_out_of_range() {
        let xs = vec![50i64];
        let res = pack_rice(&xs, 2, 10);
        assert!(matches!(res, Err(CodecError::OutOfRange)));
    }

    #[test]
    fn challenge_roundtrip() {
        let d = 32;
        let mut poly = vec![0i64; d];
        poly[0] = 1;
        poly[3] = -1;
        poly[7] = 1;
        poly[31] = -1;
        let enc = pack_challenge(&poly, 4, d).unwrap();
        let (dec, _) = unpack_challenge(&enc, 4, d).unwrap();
        assert_eq!(dec, poly);
    }

    #[test]
    fn challenge_rejects_non_pm_one_coefficients() {
        // poly[3] = 2 must be rejected, not silently normalised to +1.
        let d = 32;
        let mut poly = vec![0i64; d];
        poly[0] = 1;
        poly[3] = 2;
        poly[7] = -1;
        poly[31] = 1;
        assert!(matches!(
            pack_challenge(&poly, 4, d),
            Err(CodecError::OutOfRange)
        ));
    }

    #[test]
    fn challenge_rejects_wrong_length() {
        // poly length != d
        let short = vec![1i64, 0, -1];
        assert!(matches!(
            pack_challenge(&short, 2, 32),
            Err(CodecError::LengthMismatch)
        ));
    }

    #[test]
    fn challenge_rejects_unsorted_positions() {
        // hand-build a buffer with positions [5, 3, …] (not ascending)
        let d = 32;
        let pos_bits = bits_for_positions(d);
        let mut w = BitWriter::new();
        w.write_bits(5, pos_bits);
        w.write_bits(3, pos_bits);
        w.write_bits(0, 1);
        w.write_bits(0, 1); // sign bits
        w.pad_to_byte();
        let buf = w.into_bytes();
        let res = unpack_challenge(&buf, 2, d);
        assert!(matches!(res, Err(CodecError::NonCanonical)));
    }

    #[test]
    fn pk_roundtrip_synthetic() {
        let par = TEST_PARAMS;
        let codec = LoTRSCodec::new(par);
        // synthesise a pk with coefficients near [0, q)
        let pk: Vec<Vec<u64>> = (0..par.k)
            .map(|i| {
                (0..par.d)
                    .map(|j| ((i as u64 * 31 + j as u64 * 997) % par.q))
                    .collect()
            })
            .collect();
        let enc = codec.pk_encode(&pk).unwrap();
        let dec = codec.pk_decode(&enc).unwrap();
        assert_eq!(dec, pk);
    }

    #[test]
    fn sig_encode_rejects_wrong_poly_length() {
        use crate::params::TEST_PARAMS;
        let par = TEST_PARAMS;
        let codec = LoTRSCodec::new(par);
        // Build a structurally valid but internally malformed signature:
        // b_bin_hi has the right number of entries, but one polynomial
        // is shorter than d — encoder must refuse it.
        let short = vec![0u64; par.d - 1]; // too short
        let full = vec![0u64; par.d];
        let mut b_bin_hi = vec![full.clone(); par.n_hat];
        b_bin_hi[0] = short;
        let sig = Signature {
            pi: Proof {
                b_bin_hi,
                w_tilde_hi: vec![], // kappa == 1
                x: vec![0u64; par.d],
                f1: vec![vec![0u64; par.d]; par.kappa * (par.beta - 1)],
                z_b: vec![vec![0u64; par.d]; par.n_hat + par.k_hat],
            },
            z_tilde: vec![vec![0u64; par.d]; par.l],
            r_tilde: vec![vec![0u64; par.d]; par.l_prime],
            e_tilde: vec![vec![0u64; par.d]; par.k],
        };
        assert!(matches!(
            codec.sig_encode(&sig),
            Err(CodecError::LengthMismatch)
        ));
    }

    #[test]
    fn pk_rejects_length_mismatch() {
        let codec = LoTRSCodec::new(TEST_PARAMS);
        let res = codec.pk_decode(&[0u8; 5]);
        assert!(matches!(res, Err(CodecError::LengthMismatch)));
    }
}
