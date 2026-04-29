//! XOF-based samplers — port of `lotrs-py/sample.py`.
//!
//! Every sampler here reads bytes from a caller-supplied SHAKE256 XOF.
//! Output is fully deterministic from the XOF seed.  This matches the
//! Python reference byte-for-byte, so the same `(seed, tag…)` pair
//! yields the same polynomial in Rust and Python.
//!
//! **Gaussian sampling — two backends.**  CDT tables are built from
//! mpmath high-precision CDFs in the Python reference.  Rust does not
//! ship a high-precision library by default, so we consume CDTs as
//! precomputed `&[u128]` slices provided by the caller (see
//! [`crate::cdt`]).  The CDT sampling algorithm (binary search +
//! separate sign byte) in [`xof_sample_gaussian`] is identical to
//! `xof_sample_gaussian` in Python.
//!
//! For the two mask widths `sigma_0` / `sigma_0_prime` at
//! `BENCH_4OF32` / `BENCH_PARAMS` / `PRODUCTION_PARAMS` the CDT would
//! need `10^7`–`10^8` entries, so we use the FACCT-style integer
//! pipeline — [`prepare_facct`] builds the per-sigma constants once,
//! [`xof_sample_gaussian_facct`] is the hot path.  Spec at
//! `../lotrs-facct-sampler.md`; agreement with the Python reference is
//! pinned by the cross-language KAT at `tests/sampler_kat.rs`.
//!
//! **Backend selection is explicit per parameter set** —
//! `LoTRSParams.mask_sampler` (`MaskSamplerKind::Cdt | Facct`) is read
//! once in `LoTRS::try_new`; no runtime `sigma`-threshold dispatch.

use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

const SHAKE256_RATE: usize = 136;
const XOF_BATCH_BLOCKS: usize = 4;
const XOF_BATCH_BYTES: usize = SHAKE256_RATE * XOF_BATCH_BLOCKS;

/// SHAKE256 XOF with an indefinite output stream.  Wraps `sha3::Shake256`
/// together with a lazily-initialized reader so callers can interleave
/// `read(n)` calls after initial absorbing.
pub struct Xof {
    reader: <Shake256 as ExtendableOutput>::Reader,
    buf: [u8; XOF_BATCH_BYTES],
    pos: usize,
    len: usize,
}

/// A tag for XOF domain separation — mirrors the Python `*tags` variadic.
/// Use [`Xof::new`] with any slice of [`Tag`] values.
#[derive(Clone, Copy, Debug)]
pub enum Tag<'a> {
    Bytes(&'a [u8]),
    Str(&'a str),
    /// 4-byte little-endian encoding.
    Int(u32),
}

impl Xof {
    /// Create a SHAKE256 XOF by absorbing `seed` followed by each `tag`,
    /// matching `lotrs-py/sample.py::make_xof(seed, *tags)`.
    pub fn new(seed: &[u8], tags: &[Tag<'_>]) -> Self {
        let mut h = Shake256::default();
        h.update(seed);
        for tag in tags {
            match *tag {
                Tag::Bytes(b) => h.update(b),
                Tag::Str(s) => h.update(s.as_bytes()),
                Tag::Int(i) => h.update(&i.to_le_bytes()),
            }
        }
        Self {
            reader: h.finalize_xof(),
            buf: [0u8; XOF_BATCH_BYTES],
            pos: 0,
            len: 0,
        }
    }

    #[inline]
    fn refill(&mut self) {
        self.reader.read(&mut self.buf);
        self.pos = 0;
        self.len = self.buf.len();
    }

    /// Read the next `n` bytes of XOF output into `out`.
    pub fn read_into(&mut self, out: &mut [u8]) {
        let mut written = 0usize;
        while written < out.len() {
            if self.pos == self.len {
                self.refill();
            }
            let avail = self.len - self.pos;
            let want = out.len() - written;
            let take = core::cmp::min(avail, want);
            out[written..written + take].copy_from_slice(&self.buf[self.pos..self.pos + take]);
            self.pos += take;
            written += take;
        }
    }

    /// Read the next `n` bytes of XOF output into a fresh `Vec`.
    pub fn read(&mut self, n: usize) -> Vec<u8> {
        let mut out = vec![0u8; n];
        self.read_into(&mut out);
        out
    }

    /// Read the next byte.  Convenience wrapper.
    pub fn read_byte(&mut self) -> u8 {
        if self.pos == self.len {
            self.refill();
        }
        let b = self.buf[self.pos];
        self.pos += 1;
        b
    }

    #[inline]
    pub fn read_u64_le(&mut self, n_bytes: usize) -> u64 {
        debug_assert!(n_bytes <= 8);
        let mut buf = [0u8; 8];
        self.read_into(&mut buf[..n_bytes]);
        u64::from_le_bytes(buf)
    }

    #[inline]
    pub fn read_u128_le(&mut self, n_bytes: usize) -> u128 {
        debug_assert!(n_bytes <= 16);
        let mut buf = [0u8; 16];
        self.read_into(&mut buf[..n_bytes]);
        u128::from_le_bytes(buf)
    }

    #[inline]
    pub fn read_array<const N: usize>(&mut self) -> [u8; N] {
        let mut out = [0u8; N];
        self.read_into(&mut out);
        out
    }
}

// ---- little helpers ------------------------------------------------------

#[inline]
fn bits_to_hold(n: u64) -> u32 {
    if n == 0 {
        0
    } else {
        64 - n.leading_zeros()
    }
}

#[inline]
fn read_le_u64(xof: &mut Xof, n_bytes: usize) -> u64 {
    xof.read_u64_le(n_bytes)
}

// ---- polynomial samplers ------------------------------------------------

/// Uniform polynomial in `[0, q)^d` via rejection sampling.  Mirrors
/// `xof_sample_uniform` in Python.
pub fn xof_sample_uniform(xof: &mut Xof, q: u64, d: usize) -> Vec<u64> {
    let bits_q = bits_to_hold(q - 1);
    let byte_q = ((bits_q + 7) / 8) as usize;
    let mask_q: u64 = if bits_q == 64 {
        u64::MAX
    } else {
        (1u64 << bits_q) - 1
    };

    let mut coeffs = Vec::with_capacity(d);
    while coeffs.len() < d {
        let x = read_le_u64(xof, byte_q) & mask_q;
        if x < q {
            coeffs.push(x);
        }
    }
    coeffs
}

/// Uniform polynomial with coefficients in `[-eta, eta]`.  Internally
/// returns unsigned representatives in `[0, q)` is not done here —
/// callers typically want the centred form directly, so we return
/// `Vec<i64>`.  Mirrors `xof_sample_short`.
pub fn xof_sample_short(xof: &mut Xof, eta: i64, d: usize) -> Vec<i64> {
    debug_assert!(eta >= 0);
    let width = (2 * eta + 1) as u64;
    let bits = bits_to_hold(width);
    let nbytes = ((bits + 7) / 8) as usize;
    let mask = if bits == 64 {
        u64::MAX
    } else {
        (1u64 << bits) - 1
    };

    let mut coeffs = Vec::with_capacity(d);
    while coeffs.len() < d {
        let x = read_le_u64(xof, nbytes) & mask;
        if x < width {
            coeffs.push(x as i64 - eta);
        }
    }
    coeffs
}

/// Uniform ternary polynomial (coefficients in `{-1, 0, 1}`).  Reads 2
/// bits per coefficient with rejection of the bit-pair `11`.  Mirrors
/// `xof_sample_ternary`.
pub fn xof_sample_ternary(xof: &mut Xof, d: usize) -> Vec<i64> {
    let mut coeffs = Vec::with_capacity(d);
    while coeffs.len() < d {
        let b = xof.read_byte();
        for shift in (0..8).step_by(2) {
            if coeffs.len() >= d {
                break;
            }
            let val = (b >> shift) & 0b11;
            if val < 3 {
                coeffs.push(val as i64 - 1);
            }
        }
    }
    coeffs
}

/// Sample the challenge polynomial `x ∈ C = { x : ‖x‖_∞ = 1, ‖x‖_1 = w }`.
/// Returns signed coefficients in `{-1, 0, 1}` of length `d`, with exactly
/// `w` non-zero entries.  Matches `xof_sample_challenge` byte-for-byte.
pub fn xof_sample_challenge(xof: &mut Xof, w: usize, d: usize) -> Vec<i64> {
    let pos_bits = core::cmp::max(1, bits_to_hold((d - 1) as u64));
    let pos_bytes = ((pos_bits + 7) / 8) as usize;
    let pos_mask: u64 = if pos_bits == 64 {
        u64::MAX
    } else {
        (1u64 << pos_bits) - 1
    };

    let mut positions: Vec<usize> = Vec::with_capacity(w);
    let mut seen = vec![false; d];
    while positions.len() < w {
        let raw = read_le_u64(xof, pos_bytes) & pos_mask;
        if raw >= d as u64 {
            continue;
        } // reject to avoid bias
        let idx = raw as usize;
        if !seen[idx] {
            seen[idx] = true;
            positions.push(idx);
        }
    }

    let sign_nbytes = (w + 7) / 8;
    let mut sign_bytes = vec![0u8; sign_nbytes];
    xof.read_into(&mut sign_bytes);
    // Python: int.from_bytes(sign_bytes, "little") so bit i is
    //    (sign_bytes[i/8] >> (i%8)) & 1
    let mut coeffs = vec![0i64; d];
    for (i, &pos) in positions.iter().enumerate() {
        let bit = (sign_bytes[i / 8] >> (i % 8)) & 1;
        coeffs[pos] = if bit == 1 { -1 } else { 1 };
    }
    coeffs
}

/// Discrete Gaussian `D_sigma` via CDT look-up.  The CDT must be
/// precomputed (scaled to `2^lam`, last entry clamped to `2^lam` so the
/// binary search always terminates).  Matches `xof_sample_gaussian`.
pub fn xof_sample_gaussian(xof: &mut Xof, cdt: &[u128], lam: u32, d: usize) -> Vec<i64> {
    assert!(lam % 8 == 0, "lam must be a whole number of bytes");
    let lam_bytes = (lam / 8) as usize;
    debug_assert!(lam_bytes <= 16, "lam <= 128 bits supported");

    let mut coeffs = Vec::with_capacity(d);
    for _ in 0..d {
        let u = xof.read_u128_le(lam_bytes);

        // separate sign byte: only LSB is used
        let sign = xof.read_byte() & 1;

        // binary search: smallest k with cdt[k] > u
        let mut lo = 0usize;
        let mut hi = cdt.len() - 1;
        while lo < hi {
            let mid = (lo + hi) >> 1;
            if cdt[mid] <= u {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        let k = lo as i64;
        coeffs.push(if k > 0 && sign == 1 { -k } else { k });
    }
    coeffs
}

// =========================================================================
//  FACCT-style large-sigma Gaussian sampler
// =========================================================================
//
//  Port of the FACCT-style path in `lotrs-py/sample.py`, driven by the
//  spec in `../lotrs-facct-sampler.md`.  The target law is the truncated
//  discrete Gaussian on `[-ceil(14*sigma), +ceil(14*sigma)]`; runtime is
//  integer-only.  Acceptance uses the decomposition
//
//      u = x^2 / (2 sigma^2)  =  k * ln 2 + r          (0 <= r < ln 2)
//      exp(-u)  =  2^-k * exp(-r)
//
//  and the two independent Bernoulli trials
//
//      Ber(2^-k)   :  k fresh XOF bits must all be 0
//      Ber(exp(-r)) :  degree-20 Q(64) Horner polynomial vs. a 64-bit coin
//
//  Critical portability points:
//
//  * The Horner accumulator must be a signed 128-bit integer.  Several
//    polynomial coefficients are negative and intermediate sums routinely
//    cross zero; using `u128` would logical-shift them and corrupt the
//    threshold (see the cross-language KAT for a regression test).
//  * The 21 polynomial coefficients are shipped verbatim (`FACCT_EXP_POLY_Q64`).
//    The Python side also freezes these and refuses mismatched prep args,
//    so both ends agree integer-for-integer independent of any float path.

/// Q(u_bits) precision of `u_q` / `ln2_q` / `r_q`.
pub const FACCT_U_BITS: u32 = 48;
/// Q(prob_bits) precision of Bernoulli thresholds.
pub const FACCT_PROB_BITS: u32 = 64;
/// Taylor degree for `exp(-r)` on `[0, ln 2)`.
pub const FACCT_EXP_POLY_DEGREE: usize = 20;

/// Frozen Q(64) coefficients of `(-1)^i / i!` for `i = 0..=20`.  These
/// match `lotrs-py/sample.py::FACCT_EXP_POLY_Q64` exactly, so both
/// implementations compute the same threshold from the same `r_q`.
///
/// `coeff[0] = 2^64` and `coeff[1] = -2^64` fit in `i128` but not in
/// any 64-bit type.
pub const FACCT_EXP_POLY_Q64: &[i128] = &[
    18_446_744_073_709_551_616_i128,
    -18_446_744_073_709_551_616_i128,
    9_223_372_036_854_775_808_i128,
    -3_074_457_345_618_258_432_i128,
    768_614_336_404_564_608_i128,
    -153_722_867_280_912_928_i128,
    25_620_477_880_152_156_i128,
    -3_660_068_268_593_165_i128,
    457_508_533_574_146_i128,
    -50_834_281_508_238_i128,
    5_083_428_150_824_i128,
    -462_129_831_893_i128,
    38_510_819_324_i128,
    -2_962_370_717_i128,
    211_597_908_i128,
    -14_106_527_i128,
    881_658_i128,
    -51_862_i128,
    2_881_i128,
    -152_i128,
    8_i128,
];

/// Prepared fixed-point state for the FACCT-style sampler at a given
/// sigma.  One instance per sigma is enough; callers typically cache
/// these alongside the parameter set.
#[derive(Clone, Debug)]
pub struct FacctParams {
    /// `ceil(14 * sigma)` — truncation bound.  Proposals are uniform on
    /// `[-tail, +tail]`.
    pub tail: u64,
    /// `round(sigma * sigma)` — integer part of sigma².
    pub sigma_sq_int: u128,
    /// `round(ln(2) * 2^U_BITS)` — Q(48) representation of `ln 2`.
    pub ln2_q: u128,
}

/// Prepare FACCT fixed-point constants for a given `sigma`.
///
/// Returns `None` if `sigma` is not a finite positive `f64`, matching
/// Python's `_validate_sigma`: NaN, ±infinity, zero, and negative
/// values are all cleanly rejected rather than panicking or silently
/// saturating.
pub fn prepare_facct(sigma: f64) -> Option<FacctParams> {
    if !(sigma.is_finite() && sigma > 0.0) {
        return None;
    }
    let tail = (14.0 * sigma).ceil() as u64;
    let sigma_sq_int = (sigma * sigma).round() as u128;
    // `2^FACCT_U_BITS = 2^48` is exact in f64; ln(2) is a constant; so
    // `ln(2) * 2^48` is evaluated in f64 with full precision and then
    // rounded to an integer.  Matches Python's
    //   int(round(math.log(2.0) * (1 << u_bits))).
    let ln2_q = ((2.0_f64.ln()) * ((1u128 << FACCT_U_BITS) as f64)).round() as u128;
    Some(FacctParams {
        tail,
        sigma_sq_int,
        ln2_q,
    })
}

/// `Ber(2^-k)` — return `true` with probability `2^-k` by requiring `k`
/// XOF bits to be zero.  Matches the trimmed-byte partial-chunk layout
/// in `_xof_bernoulli_power_of_two` in the Python reference.
fn xof_bernoulli_power_of_two(xof: &mut Xof, k: usize) -> bool {
    let mut k = k;
    while k >= 64 {
        if xof.read_u64_le(8) != 0 {
            return false;
        }
        k -= 64;
    }
    if k == 0 {
        return true;
    }
    let nbytes = (k + 7) / 8;
    let v = xof.read_u64_le(nbytes);
    let mask = (1u64 << k) - 1;
    (v & mask) == 0
}

/// Horner-evaluate the Q(64) Taylor polynomial for `exp(-r)` at
/// `r_q / 2^U_BITS` and return a `u64` Bernoulli threshold.  The value
/// is clamped into `[0, 2^64)`.
///
/// **`acc` must be `i128`**, not `u128` — several coefficients are
/// negative and intermediate sums can cross zero.  An unsigned
/// accumulator would logical-shift those values and corrupt the result.
fn facct_exp_neg_threshold(r_q: u128) -> u64 {
    let u_bits = FACCT_U_BITS;
    let poly = FACCT_EXP_POLY_Q64;
    let mut acc: i128 = poly[poly.len() - 1];
    // r_q < ln2_q < 2^48, so `acc * (r_q as i128)` is bounded by
    // |acc| * 2^48.  |acc| stays <= 2^64 during evaluation (bound
    // established empirically and checked in tests), so the product
    // is bounded by 2^112 and fits comfortably in i128.
    let r_q_i = r_q as i128;
    for i in (0..poly.len() - 1).rev() {
        acc = poly[i] + ((acc * r_q_i) >> u_bits);
    }
    if acc < 0 {
        0
    } else if acc >= (1i128 << FACCT_PROB_BITS) {
        u64::MAX
    } else {
        acc as u64
    }
}

/// Draw `d` independent samples from the truncated discrete Gaussian at
/// `sigma`, using the FACCT-style integer pipeline.  Matches
/// `xof_sample_gaussian_facct` in `lotrs-py/sample.py` byte-for-byte
/// given identical XOF input.
pub fn xof_sample_gaussian_facct(xof: &mut Xof, prepared: &FacctParams, d: usize) -> Vec<i64> {
    let tail = prepared.tail;
    let tail_i = tail as i64;
    let sigma_sq_int = prepared.sigma_sq_int;
    let ln2_q = prepared.ln2_q;
    let u_bits = FACCT_U_BITS;

    let mut out = Vec::with_capacity(d);
    while out.len() < d {
        // Signed-uniform proposal on  [-tail, +tail]  (2*tail + 1 values).
        let unsigned = xof_randbelow(xof, 2 * tail + 1);
        let x = (unsigned as i64) - tail_i;

        // u_q  =  floor(x^2 * 2^u_bits / (2 * sigma^2_int))
        // x^2 fits in u128 (tail up to ~10^9, tail^2 ~ 10^18 < 2^64).
        let x_sq = (x as i128 * x as i128) as u128;
        let u_q = (x_sq << u_bits) / (2u128 * sigma_sq_int);

        // k, r_q from  u_q = k * ln2_q + r_q.
        let k = (u_q / ln2_q) as usize;
        let r_q = u_q - (k as u128) * ln2_q;

        if !xof_bernoulli_power_of_two(xof, k) {
            continue;
        }
        let threshold = facct_exp_neg_threshold(r_q);
        let coin = xof.read_u64_le(8);
        if coin >= threshold {
            continue;
        }
        out.push(x);
    }
    out
}

/// Uniform integer in `[0, bound)` via rejection on a power-of-2 mask.
/// Matches `_xof_randbelow` in the Python reference.
fn xof_randbelow(xof: &mut Xof, bound: u64) -> u64 {
    assert!(bound > 0);
    let bits = core::cmp::max(1, 64 - (bound - 1).leading_zeros());
    let nbytes = ((bits + 7) / 8) as usize;
    let mask: u64 = if bits == 64 {
        u64::MAX
    } else {
        (1u64 << bits) - 1
    };
    loop {
        let x = xof.read_u64_le(nbytes) & mask;
        if x < bound {
            return x;
        }
    }
}

// =========================================================================
//  Unified Gaussian sampler dispatch
// =========================================================================

/// Backend-agnostic prepared Gaussian sampler — either a compile-time
/// CDT (small sigma) or runtime-built FACCT parameters (large sigma).
pub enum GaussianSampler {
    /// Exact discrete-Gaussian CDT scaled to `lam` bits.
    Cdt { cdt: &'static [u128], lam: u32 },
    /// FACCT-style integer-only pipeline.
    Facct(FacctParams),
}

impl GaussianSampler {
    /// Draw `d` samples from the prepared Gaussian.  Transparent over
    /// the CDT / FACCT variants.
    pub fn sample(&self, xof: &mut Xof, d: usize) -> Vec<i64> {
        match self {
            Self::Cdt { cdt, lam } => xof_sample_gaussian(xof, cdt, *lam, d),
            Self::Facct(p) => xof_sample_gaussian_facct(xof, p, d),
        }
    }
}

// =========================================================================
//  Rejection sampling  (Fig. 1)
// =========================================================================

/// `Rej(z, v, phi, K)` — baseline rejection sampling.  Returns `true`
/// on accept (corresponds to output 0 in the paper, meaning *continue*)
/// and `false` on reject (output 1, meaning *abort / restart*).
///
/// The f64 computations match Python's `math.log` / `math.exp` bit-for-
/// bit on any IEEE 754 platform, which is what gives deterministic
/// accept / reject decisions across the Rust and Python reference.
pub fn rej(xof: &mut Xof, z_flat: &[i64], v_flat: &[i64], phi: f64, k_bound: f64) -> bool {
    debug_assert_eq!(z_flat.len(), v_flat.len());
    let sigma = phi * k_bound;
    let sigma_sq = sigma * sigma;
    let mu_phi = (12.0 / phi + 1.0 / (2.0 * phi * phi)).exp();

    let inner_zv = inner_product_i128(z_flat, v_flat);
    let norm_v_sq = norm_sq_i128(v_flat);

    let log_p = (-2.0 * inner_zv as f64 + norm_v_sq as f64) / (2.0 * sigma_sq);
    let log_threshold = log_p - mu_phi.ln();

    let u = sample_u01(xof);
    if u.ln() >= log_threshold {
        return false;
    }

    let max_abs = z_flat.iter().map(|&c| c.unsigned_abs()).max().unwrap_or(0);
    if (max_abs as f64) > 6.0 * sigma {
        return false;
    }
    true
}

/// `RejOp(z, v, phi, K)` — optimised rejection sampling that exploits
/// `<z, v> ≥ 0` by refusing `<z, v> < 0` outright.  Same return
/// convention as [`rej`].
pub fn rej_op(xof: &mut Xof, z_flat: &[i64], v_flat: &[i64], phi: f64, k_bound: f64) -> bool {
    debug_assert_eq!(z_flat.len(), v_flat.len());
    let inner_zv = inner_product_i128(z_flat, v_flat);
    if inner_zv < 0 {
        return false;
    }

    let sigma = phi * k_bound;
    let sigma_sq = sigma * sigma;
    let mu_phi = (1.0 / (2.0 * phi * phi)).exp();

    let norm_v_sq = norm_sq_i128(v_flat);
    let log_p = (-2.0 * inner_zv as f64 + norm_v_sq as f64) / (2.0 * sigma_sq);
    let log_threshold = log_p - mu_phi.ln();

    let u = sample_u01(xof);
    if u.ln() >= log_threshold {
        return false;
    }

    let max_abs = z_flat.iter().map(|&c| c.unsigned_abs()).max().unwrap_or(0);
    if (max_abs as f64) > 6.0 * sigma {
        return false;
    }
    true
}

#[inline]
fn inner_product_i128(a: &[i64], b: &[i64]) -> i128 {
    a.iter()
        .zip(b.iter())
        .map(|(&x, &y)| x as i128 * y as i128)
        .sum()
}

#[inline]
fn norm_sq_i128(v: &[i64]) -> i128 {
    v.iter()
        .map(|&x| {
            let xi = x as i128;
            xi * xi
        })
        .sum()
}

/// Sample `u ∈ (0, 1)` as a 53-bit mantissa in `[0, 1)`, with a hard
/// floor of `2^-53` that matches Python's treatment of the `u == 0.0`
/// case.  Reads 8 XOF bytes.
#[inline]
fn sample_u01(xof: &mut Xof) -> f64 {
    let u_raw = xof.read_u64_le(8);
    let mantissa = u_raw & ((1u64 << 53) - 1);
    let u = mantissa as f64 / ((1u64 << 53) as f64);
    if u == 0.0 {
        (-53f64).exp2()
    } else {
        u
    }
}

// -------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn determinism() {
        let seed = [0u8; 32];
        let mut x1 = Xof::new(&seed, &[Tag::Bytes(b"test")]);
        let mut x2 = Xof::new(&seed, &[Tag::Bytes(b"test")]);
        assert_eq!(
            xof_sample_uniform(&mut x1, 7681, 8),
            xof_sample_uniform(&mut x2, 7681, 8)
        );
    }

    #[test]
    fn challenge_weight_and_linf() {
        let seed = [1u8; 32];
        let mut x = Xof::new(&seed, &[Tag::Bytes(b"chal")]);
        let ch = xof_sample_challenge(&mut x, 10, 64);
        let l1: i64 = ch.iter().map(|c| c.abs()).sum();
        let linf: i64 = ch.iter().map(|c| c.abs()).max().unwrap_or(0);
        assert_eq!(l1, 10);
        assert_eq!(linf, 1);
    }

    #[test]
    fn short_bounds() {
        let seed = [3u8; 32];
        let mut x = Xof::new(&seed, &[Tag::Bytes(b"short")]);
        let s = xof_sample_short(&mut x, 3, 64);
        for c in s {
            assert!(c.abs() <= 3);
        }
    }

    #[test]
    fn ternary_range() {
        let seed = [4u8; 32];
        let mut x = Xof::new(&seed, &[Tag::Bytes(b"ter")]);
        let t = xof_sample_ternary(&mut x, 64);
        for c in t {
            assert!(c.abs() <= 1);
        }
    }

    #[test]
    fn gaussian_uses_bundled_cdt() {
        // Smoke test: the Gaussian sampler accepts the CDT from crate::cdt
        // and produces reasonable-looking samples.  Full byte-for-byte
        // interop against Python is covered in tests/interop.rs.
        use crate::cdt::CDT_SIGMA_A_TEST;
        let mut x = Xof::new(&[5u8; 32], &[Tag::Bytes(b"gauss")]);
        let g = xof_sample_gaussian(&mut x, CDT_SIGMA_A_TEST, 128, 32);
        // 14 * sigma_a is the nominal tail — absolutely bounded by the
        // CDT size, so generated magnitudes cannot exceed it.
        let bound = CDT_SIGMA_A_TEST.len() as i64;
        for c in g {
            assert!(c.abs() < bound);
        }
    }

    #[test]
    fn facct_frozen_coefficients_match_spec() {
        // Sanity: there are exactly degree+1 coefficients.
        assert_eq!(FACCT_EXP_POLY_Q64.len(), FACCT_EXP_POLY_DEGREE + 1);
        // Head / tail anchor values.
        assert_eq!(FACCT_EXP_POLY_Q64[0], 1_i128 << 64);
        assert_eq!(FACCT_EXP_POLY_Q64[1], -(1_i128 << 64));
        assert_eq!(FACCT_EXP_POLY_Q64[FACCT_EXP_POLY_DEGREE], 8);
    }

    #[test]
    fn facct_poly_matches_exp_neg_r_to_2_to_the_minus_48() {
        // max error on [0, ln 2) should be < 3e-15 — matches the
        // Python test threshold.  Sweep 257 equally-spaced points.
        let scale = 1u128 << FACCT_U_BITS;
        let ln2 = 2.0_f64.ln();
        for i in 0..=256 {
            let r = ln2 * i as f64 / 256.0;
            let r_q = (r * scale as f64).round() as u128;
            let thr = facct_exp_neg_threshold(r_q);
            let got = thr as f64 / (1u128 << FACCT_PROB_BITS) as f64;
            let want = (-r).exp();
            let err = (got - want).abs();
            assert!(err < 3e-15, "at r={r}: err {err} (got {got}, want {want})");
        }
    }

    #[test]
    fn facct_bernoulli_pow_of_two_zero_always_accepts() {
        // k=0 must accept without consuming XOF bytes.
        let mut x = Xof::new(&[0u8; 32], &[]);
        let before = x.read(0); // empty read, no-op sanity
        assert!(xof_bernoulli_power_of_two(&mut x, 0));
        drop(before);
    }

    #[test]
    fn facct_determinism_at_moderate_sigma() {
        // Same seed + same prepared params => identical samples.
        let params = prepare_facct(100.0).expect("finite positive sigma");
        let seed = [0x37u8; 32];
        let mut a = Xof::new(&seed, &[Tag::Bytes(b"kat")]);
        let mut b = Xof::new(&seed, &[Tag::Bytes(b"kat")]);
        let s1 = xof_sample_gaussian_facct(&mut a, &params, 256);
        let s2 = xof_sample_gaussian_facct(&mut b, &params, 256);
        assert_eq!(s1, s2);
    }

    #[test]
    fn facct_samples_sit_within_truncation_bound() {
        let sigma = 100.0_f64;
        let params = prepare_facct(sigma).unwrap();
        let mut x = Xof::new(&[0x4eu8; 32], &[Tag::Bytes(b"bound")]);
        let s = xof_sample_gaussian_facct(&mut x, &params, 2048);
        let tail = params.tail as i64;
        for c in &s {
            assert!(c.abs() <= tail, "sample {c} exceeds tail {tail}");
        }
    }

    #[test]
    fn facct_samples_are_roughly_centred() {
        let sigma = 100.0_f64;
        let params = prepare_facct(sigma).unwrap();
        let mut x = Xof::new(&[0x5eu8; 32], &[Tag::Bytes(b"ctr")]);
        let s = xof_sample_gaussian_facct(&mut x, &params, 4096);
        let mean: f64 = s.iter().map(|&c| c as f64).sum::<f64>() / s.len() as f64;
        // Within 3σ/√N of zero.
        let bound = 3.0 * sigma / (s.len() as f64).sqrt();
        assert!(mean.abs() < bound, "mean {mean} not centred (±{bound})");
    }

    #[test]
    fn prepare_facct_rejects_non_finite_sigma() {
        // Matches Python's `_validate_sigma`: NaN, ±inf, zero, negative
        // all produce a clean `None` rather than a panic or saturated
        // output.
        assert!(prepare_facct(f64::NAN).is_none());
        assert!(prepare_facct(f64::INFINITY).is_none());
        assert!(prepare_facct(f64::NEG_INFINITY).is_none());
        assert!(prepare_facct(0.0).is_none());
        assert!(prepare_facct(-1.0).is_none());
        assert!(prepare_facct(100.0).is_some());
    }

    #[test]
    fn tag_int_is_4_byte_le() {
        // Python: h.update(tag.to_bytes(4, "little"))
        let mut xa = Xof::new(&[0u8; 32], &[Tag::Int(0x0403_0201)]);
        let mut xb = Xof::new(&[0u8; 32], &[Tag::Bytes(&[0x01, 0x02, 0x03, 0x04])]);
        assert_eq!(xa.read(16), xb.read(16));
    }
}
