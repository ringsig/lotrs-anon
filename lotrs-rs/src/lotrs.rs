//! LoTRS scheme — port of `lotrs-py/lotrs.py`.
//!
//! * [`LoTRS::setup`]  — public parameters (`pp` = 32-byte seed).
//! * [`LoTRS::keygen`] — per-signer key generation (byte-exact vs Python).
//! * [`LoTRS::kagg`]   — column-aggregated ring public keys.
//! * [`LoTRS::sign`]   — full two-round rejection-sampling ceremony,
//!   byte-exact with the Python reference under identical seeds.
//! * [`LoTRS::verify`] — full verifier (Fig. 6).  Returns `false` on
//!   every malformed-input class without surfacing error detail.
//!
//! The Gaussian mask backend (`MaskSampler::Cdt | Facct`) is resolved
//! once in [`LoTRS::try_new`] from the explicit `mask_sampler` field on
//! [`LoTRSParams`] — there is no runtime threshold dispatch.

use std::time::{Duration, Instant};

use crate::cdt;
use crate::codec::{LoTRSCodec, Proof as _Proof, Signature};
use crate::params::{LoTRSParams, MaskSamplerKind};
use crate::ring::{NttPolyMat, Poly, PolyMat, PolyVec, Ring};
use crate::sample::{
    prepare_facct, rej, rej_op, xof_sample_challenge, xof_sample_gaussian,
    xof_sample_gaussian_facct, xof_sample_short, xof_sample_uniform, FacctParams, Tag, Xof,
};

/// Wall-clock breakdown for a successful [`LoTRS::sign_with_timings`] call.
///
/// All four duration fields sum over every rejection-sampling attempt
/// the call made (including restarted attempts), so
/// `sign1 + sign2_rest + sign_bin + sagg` is approximately the total
/// `sign` wall-clock minus the one-shot context expansion.
///
/// The split is the obvious one: `sign_bin` is the binary ring proof
/// (Fig. 5 right column); everything else is the DualMS-style multi-
/// signature (sign1 commitments, sign2 response + rejection, sagg).
#[derive(Default, Clone, Copy, Debug)]
pub struct SignTimings {
    /// Per-signer round-1 commitments (sign1) summed over T signers
    /// and over every attempt.
    pub sign1: Duration,
    /// Round-2 work *excluding* the binary proof — w̃ aggregation,
    /// w̃₀ stability check, response z_u, BG rejection, r_u — summed
    /// over T signers and over every attempt.
    pub sign2_rest: Duration,
    /// The binary ring proof (`sign_bin`) — summed over T signers and
    /// over every attempt.  Each of the T signers runs `sign_bin`
    /// independently in this reference implementation; `sagg` then
    /// asserts they all produced the same `pi`.  An optimised multi-
    /// signer protocol would compute `sign_bin` once and broadcast,
    /// in which case the standalone RS cost is this value divided by T.
    pub sign_bin: Duration,
    /// Threshold aggregation across the T transcripts.
    pub sagg: Duration,
    /// Number of rejection-sampling attempts the call made before
    /// landing the accepted signature (≥1).
    pub attempts: u32,
}

/// Wall-clock breakdown for a successful
/// [`LoTRS::verify_decoded_with_timings`] call.
///
/// `verify_bin` covers the binary-proof checks (f1/f0/g0/g1 norm
/// bounds, A_hat_bin reconstruction & low-bit check); `verify_dualms`
/// covers everything else — the L2 bounds on z̃/r̃/ẽ, the
/// `A·z̃ + B·r̃ + ẽ` reconstruction, KAgg, the ring-keyed `pk_sum`
/// loop over N, the w̃₀ recovery, and the closing FS hash.  We attribute
/// `pk_sum` and the FS hash to DualMS because they're part of the
/// LHS=RHS check that closes the multi-signature.
#[derive(Default, Clone, Copy, Debug)]
pub struct VerifyTimings {
    pub verify_bin: Duration,
    pub verify_dualms: Duration,
}

/// Backend chosen for the two mask-sampling widths (`sigma_0`,
/// `sigma_0_prime`).  At `TEST_PARAMS` both sigmas are small and we use
/// the exact CDT.  At `BENCH_4OF32` / `BENCH_PARAMS` / `PRODUCTION_PARAMS`
/// the sigmas sit in the `10^6`–`10^7` range and a CDT would need
/// multi-GB tables, so we use the FACCT-style integer pipeline.
#[derive(Clone)]
enum MaskSampler {
    Cdt(&'static [u128]),
    Facct(FacctParams),
}

/// A LoTRS scheme context.  Holds the ring arithmetic backends, the
/// codec, and the pre-resolved Gaussian samplers for the given
/// parameter set.
pub struct LoTRS {
    pub par: LoTRSParams,
    pub r_q: Ring,
    pub r_qhat: Ring,
    pub codec: LoTRSCodec,

    // Gaussian samplers resolved at construction time.  sigma_a and
    // sigma_b are always CDT — they're small for every parameter set
    // we support.  sigma_0 and sigma_0_prime dispatch per-parameter-set
    // via `MaskSampler` because the CDT is infeasible for the large-
    // sigma mask widths at BENCH / PRODUCTION.
    cdt_sigma_a: &'static [u128],
    cdt_sigma_b: &'static [u128],
    sampler_0: MaskSampler,
    sampler_0_prime: MaskSampler,
}

impl LoTRS {
    /// Construct a scheme context for `par`.  Panics on an unsupported
    /// parameter set — if you want a total function, use
    /// [`LoTRS::try_new`] instead.
    pub fn new(par: LoTRSParams) -> Self {
        Self::try_new(par).expect("unsupported LoTRSParams")
    }

    /// Panic-free constructor.  Returns `None` for any parameter set
    /// that this implementation does not support:
    ///
    /// * `par.kappa != 1`
    /// * no shipped CDT table matches `par.sigma_a` or `par.sigma_b`
    ///
    /// Callers building a panic-free verifier must use this rather
    /// than [`LoTRS::new`].
    pub fn try_new(par: LoTRSParams) -> Option<Self> {
        if par.kappa != 1 {
            return None;
        }
        // sigma_a / sigma_b are always CDT.  CDTs are resolved by
        // (parameter-set name, width) — no length-only lookup, which
        // could pick a table built for a different distribution that
        // happened to share the same truncation length.
        let cdt_sigma_a = resolve_cdt(&par, CdtWidth::SigmaA)?;
        let cdt_sigma_b = resolve_cdt(&par, CdtWidth::SigmaB)?;

        // The mask backend is declared on the parameter set — no
        // threshold dispatch here.
        let sampler_0 = build_mask_sampler(&par, CdtWidth::Sigma0)?;
        let sampler_0_prime = build_mask_sampler(&par, CdtWidth::Sigma0Prime)?;

        Some(Self {
            r_q: Ring::new(par.q, par.d),
            r_qhat: Ring::new(par.q_hat, par.d),
            codec: LoTRSCodec::new(par),
            cdt_sigma_a,
            cdt_sigma_b,
            sampler_0,
            sampler_0_prime,
            par,
        })
    }

    // =================================================================
    //  Setup / KeyGen  (Fig. 4)
    // =================================================================

    /// `Setup(seed) -> pp`.  Public parameters are the 32-byte seed.
    pub fn setup(&self, seed: &[u8]) -> Vec<u8> {
        assert_eq!(seed.len(), 32);
        seed.to_vec()
    }

    /// Expand public matrix `A ∈ R_q^{k × l}` from the pp seed.
    #[allow(non_snake_case)]
    pub fn expand_A(&self, pp: &[u8]) -> PolyMat {
        let par = &self.par;
        let mut a = Vec::with_capacity(par.k);
        for i in 0..par.k {
            let mut row = Vec::with_capacity(par.l);
            for j in 0..par.l {
                let mut x = Xof::new(
                    pp,
                    &[Tag::Bytes(b"A"), Tag::Int(i as u32), Tag::Int(j as u32)],
                );
                row.push(xof_sample_uniform(&mut x, par.q, par.d));
            }
            a.push(row);
        }
        a
    }

    /// Expand binary-proof matrix `G ∈ R_qhat^{n_hat × G_cols}` from pp.
    #[allow(non_snake_case)]
    pub fn expand_G(&self, pp: &[u8]) -> PolyMat {
        let par = &self.par;
        let mut g = Vec::with_capacity(par.n_hat);
        for i in 0..par.n_hat {
            let mut row = Vec::with_capacity(par.G_cols());
            for j in 0..par.G_cols() {
                let mut x = Xof::new(
                    pp,
                    &[Tag::Bytes(b"G"), Tag::Int(i as u32), Tag::Int(j as u32)],
                );
                row.push(xof_sample_uniform(&mut x, par.q_hat, par.d));
            }
            g.push(row);
        }
        g
    }

    /// `KGen(pp, seed) -> (sk, pk)`.
    pub fn keygen(&self, pp: &[u8], seed: &[u8]) -> (Vec<Poly>, Vec<Poly>) {
        assert_eq!(seed.len(), 32);
        let par = &self.par;
        let rq = &self.r_q;

        let a = self.expand_A(pp);

        // all (l + k) short polys from one XOF tagged "keygen"
        let mut xof = Xof::new(seed, &[Tag::Bytes(b"keygen")]);
        let mut sk: Vec<Poly> = Vec::with_capacity(par.l + par.k);
        for _ in 0..(par.l + par.k) {
            let signed = xof_sample_short(&mut xof, par.eta as i64, par.d);
            sk.push(rq.from_centered(&signed));
        }

        // t = A · s[..l] + s[l..l+k]
        let s_left = &sk[..par.l];
        let s_right = &sk[par.l..];
        let pk_left: PolyVec = rq.mat_vec(&a, s_left);
        let pk: Vec<Poly> = pk_left
            .iter()
            .zip(s_right.iter())
            .map(|(l_row, r_row)| rq.add(l_row, r_row))
            .collect();

        (sk, pk)
    }

    // =================================================================
    //  Hash helpers
    // =================================================================

    /// Canonical 8-byte unsigned-LE serialization of the PK table
    /// (column-major, row-within-column, poly-within-row, coefficient-
    /// within-poly).  Matches `lotrs-py/lotrs.py::_pk_table_bytes`.
    pub fn pk_table_bytes(&self, pk_table: &[Vec<Vec<Poly>>]) -> Vec<u8> {
        let mut out = Vec::new();
        for col in pk_table {
            for pk in col {
                for poly in pk {
                    for &c in poly {
                        out.extend_from_slice(&c.to_le_bytes());
                    }
                }
            }
        }
        out
    }

    /// 256-bit digest of the canonical PK-table serialization.
    pub fn pk_table_hash(&self, pk_table: &[Vec<Vec<Poly>>]) -> [u8; 32] {
        let pk_bytes = self.pk_table_bytes(pk_table);
        let mut x = Xof::new(&pk_bytes, &[Tag::Bytes(b"pk")]);
        x.read_array::<32>()
    }

    /// `H_agg(H(PK), u) -> alpha_u ∈ R_q`  (challenge-set element).
    pub fn hash_agg(&self, pk_hash: &[u8; 32], u: u32) -> Poly {
        let par = &self.par;
        let mut x = Xof::new(pk_hash, &[Tag::Bytes(b"agg"), Tag::Int(u)]);
        let raw = xof_sample_challenge(&mut x, par.w, par.d);
        self.r_q.from_centered(&raw)
    }

    /// `H_com(H(PK), mu) -> B ∈ R_q^{k × l'}`.
    pub fn hash_com(&self, pk_hash: &[u8; 32], mu: &[u8]) -> PolyMat {
        let par = &self.par;
        // First XOF returns a 32-byte seed, from which the individual
        // matrix columns are expanded.
        let seed = {
            let mut x = Xof::new(pk_hash, &[Tag::Bytes(b"com"), Tag::Bytes(mu)]);
            x.read_array::<32>()
        };
        let mut b = Vec::with_capacity(par.k);
        for i in 0..par.k {
            let mut row = Vec::with_capacity(par.l_prime);
            for j in 0..par.l_prime {
                let mut x = Xof::new(
                    &seed,
                    &[Tag::Bytes(b"B"), Tag::Int(i as u32), Tag::Int(j as u32)],
                );
                row.push(xof_sample_uniform(&mut x, par.q, par.d));
            }
            b.push(row);
        }
        b
    }

    /// Fiat-Shamir hash.  Takes the message, the `A_hat^{(1)}` component,
    /// the `B_bin^{(1)}` component, the `w^{(1)}` vectors, and the
    /// 256-bit PK-table digest, and returns the challenge `x ∈ C ⊂ R_q`.
    ///
    /// Coefficients are hashed as 8-byte *signed* LE (matching Python's
    /// `int.to_bytes(8, "little", signed=True)`).  The scheme's moduli
    /// `q, q_hat < 2^40 < 2^63`, so the cast is lossless.
    pub fn hash_fs(
        &self,
        mu: &[u8],
        a_hi: &[Poly],
        b_hi: &[Poly],
        w_hi: &[Vec<Poly>],
        pk_hash: &[u8; 32],
    ) -> Poly {
        use sha3::{
            digest::{ExtendableOutput, Update, XofReader},
            Shake256,
        };

        let mut h = Shake256::default();
        // make_xof(b"", b"FS")
        h.update(b"FS");
        h.update(mu);
        for comp in a_hi {
            for &c in comp {
                h.update(&(c as i64).to_le_bytes());
            }
        }
        for comp in b_hi {
            for &c in comp {
                h.update(&(c as i64).to_le_bytes());
            }
        }
        for j_vec in w_hi {
            for comp in j_vec {
                for &c in comp {
                    h.update(&(c as i64).to_le_bytes());
                }
            }
        }
        h.update(pk_hash);

        // The resulting XOF is then used to sample the challenge.
        // We can't stream this directly through `Xof::new` because the
        // XOF state is already populated, so use the challenge sampler
        // over the reader directly.
        let mut reader = h.finalize_xof();
        let par = &self.par;
        let pos_bits = core::cmp::max(1, 64 - ((par.d as u64 - 1).leading_zeros()));
        let pos_bytes = ((pos_bits + 7) / 8) as usize;
        let pos_mask: u64 = if pos_bits == 64 {
            u64::MAX
        } else {
            (1u64 << pos_bits) - 1
        };

        let mut positions: Vec<usize> = Vec::with_capacity(par.w);
        let mut seen = vec![false; par.d];
        while positions.len() < par.w {
            let mut buf = [0u8; 8];
            reader.read(&mut buf[..pos_bytes]);
            let raw = u64::from_le_bytes(buf) & pos_mask;
            if raw >= par.d as u64 {
                continue;
            }
            let idx = raw as usize;
            if !seen[idx] {
                seen[idx] = true;
                positions.push(idx);
            }
        }
        let sign_nbytes = (par.w + 7) / 8;
        let mut sign_bytes = vec![0u8; sign_nbytes];
        reader.read(&mut sign_bytes);

        let mut coeffs = vec![0i64; par.d];
        for (i, &pos) in positions.iter().enumerate() {
            let bit = (sign_bytes[i / 8] >> (i % 8)) & 1;
            coeffs[pos] = if bit == 1 { -1 } else { 1 };
        }
        self.r_q.from_centered(&coeffs)
    }

    // =================================================================
    //  kagg  (Fig. 4 — KAgg)
    // =================================================================

    /// `KAgg(PK) -> [t̃_0 .. t̃_{N-1}]`.
    pub fn kagg(&self, pk_table: &[Vec<Vec<Poly>>]) -> Vec<PolyVec> {
        let pk_hash = self.pk_table_hash(pk_table);
        self.kagg_with_hash(pk_table, &pk_hash)
    }

    fn kagg_with_hash(&self, pk_table: &[Vec<Vec<Poly>>], pk_hash: &[u8; 32]) -> Vec<PolyVec> {
        let par = &self.par;
        // α_u depends only on (H(PK), u), so precompute once per u
        // rather than N·T times.
        let alphas: Vec<Poly> = (0..par.T)
            .map(|u| self.hash_agg(pk_hash, u as u32))
            .collect();
        let mut out = Vec::with_capacity(par.N());
        for i in 0..par.N() {
            let mut acc: PolyVec = self.r_q.vec_zero(par.k);
            for u in 0..par.T {
                // acc += α_u · pk_table[i][u]  (fused, no intermediate Vec)
                self.r_q
                    .vec_add_scaled(&mut acc, &alphas[u], &pk_table[i][u]);
            }
            out.push(acc);
        }
        out
    }

    // =================================================================
    //  Misc polynomial helpers used by verify
    // =================================================================

    /// Decompose `n` into `num_digits` base-`base` digits, least-
    /// significant first.  Matches `_base_digits` in Python.
    pub fn base_digits(n: usize, base: usize, num_digits: usize) -> Vec<usize> {
        let mut digits = Vec::with_capacity(num_digits);
        let mut nn = n;
        for _ in 0..num_digits {
            digits.push(nn % base);
            nn /= base;
        }
        digits
    }

    /// Reconstruct full `f_{j, i}` from `f_1` (the `i >= 1` entries)
    /// and the challenge `x`:  `f_{j, 0} = x − sum_{i >= 1} f_{j, i}`.
    pub fn reconstruct_f(&self, x: &[u64], f1_flat: &[Poly]) -> Vec<Vec<Poly>> {
        let par = &self.par;
        let mut idx = 0usize;
        let mut f_full: Vec<Vec<Poly>> = Vec::with_capacity(par.kappa);
        for _ in 0..par.kappa {
            let mut row: Vec<Option<Poly>> = vec![None; par.beta];
            let mut acc = self.r_q.zero();
            for i in 1..par.beta {
                row[i] = Some(f1_flat[idx].clone());
                acc = self.r_q.add(&acc, &f1_flat[idx]);
                idx += 1;
            }
            row[0] = Some(self.r_q.sub(x, &acc));
            let row: Vec<Poly> = row.into_iter().map(|o| o.unwrap()).collect();
            f_full.push(row);
        }
        f_full
    }

    /// `g_{j,i} = f_{j,i} * (x − f_{j,i})` computed in signed integer
    /// arithmetic and reduced into `R_qhat`.  Returns `(g0, g1)`
    /// where `g0[j] = g_{j,0}` and `g1` contains the flattened
    /// `g_{j, i>=1}` entries in order.
    pub fn compute_g_qhat(&self, x: &[u64], f_full: &[Vec<Poly>]) -> (PolyVec, PolyVec) {
        let par = &self.par;
        let x_signed = self.r_q.centered(x);
        let mut g0 = Vec::with_capacity(par.kappa);
        let mut g1 = Vec::new();
        for row in f_full {
            for (i, f) in row.iter().enumerate() {
                let f_signed = self.r_q.centered(f);
                let diff: Vec<i64> = (0..par.d).map(|c| x_signed[c] - f_signed[c]).collect();
                let g_signed = poly_mul_signed(&f_signed, &diff, par.d);
                let g_qhat = self.r_qhat.from_centered(&g_signed);
                if i == 0 {
                    g0.push(g_qhat);
                } else {
                    g1.push(g_qhat);
                }
            }
        }
        (g0, g1)
    }

    /// Flatten `g0, g1` into the order expected by G multiplication.
    pub fn interleave_g(&self, g0: &[Poly], g1: &[Poly]) -> PolyVec {
        let par = &self.par;
        let mut out = Vec::new();
        let mut idx1 = 0usize;
        for j in 0..par.kappa {
            out.push(g0[j].clone());
            for _ in 1..par.beta {
                out.push(g1[idx1].clone());
                idx1 += 1;
            }
        }
        out
    }

    /// Flatten `f_full[j][i]` into a single `Vec<Poly>`.
    pub fn flatten_f(&self, f_full: &[Vec<Poly>]) -> PolyVec {
        let mut out = Vec::new();
        for row in f_full {
            for p in row {
                out.push(p.clone());
            }
        }
        out
    }

    // =================================================================
    //  Verify  (Fig. 6)
    // =================================================================

    /// Entry-point verifier.  **Returns `false` on every malformed-
    /// input class** (length mismatch, out-of-range coefficient, nonzero
    /// padding bit, trailing bytes, oversized Rice unary runs, norm
    /// violation, Fiat–Shamir mismatch, …) without panicking.
    ///
    /// `pk_table_bytes[col][row]` is the byte-encoded public key at that
    /// column / row of the ring.  Every entry is decoded and validated
    /// internally; any decode failure surfaces as a plain `false`.
    pub fn verify(
        &self,
        pp: &[u8],
        mu: &[u8],
        sig_bytes: &[u8],
        pk_table_bytes: &[Vec<Vec<u8>>],
    ) -> bool {
        self.verify_with_timings(pp, mu, sig_bytes, pk_table_bytes)
            .0
    }

    /// Same as [`LoTRS::verify`] but additionally returns a
    /// [`VerifyTimings`] breakdown.  The returned timings cover only
    /// the inner verifier; signature / pk decoding is not measured.
    /// On any decoding error the returned bool is `false` and the
    /// timings are zero.
    pub fn verify_with_timings(
        &self,
        pp: &[u8],
        mu: &[u8],
        sig_bytes: &[u8],
        pk_table_bytes: &[Vec<Vec<u8>>],
    ) -> (bool, VerifyTimings) {
        let zero = VerifyTimings::default();
        if pp.len() != 32 {
            return (false, zero);
        }

        let par = &self.par;
        if pk_table_bytes.len() != par.N() {
            return (false, zero);
        }

        let mut pk_table: Vec<Vec<Vec<Poly>>> = Vec::with_capacity(par.N());
        for col in pk_table_bytes {
            if col.len() != par.T {
                return (false, zero);
            }
            let mut col_pks = Vec::with_capacity(par.T);
            for pk_bytes in col {
                match self.codec.pk_decode(pk_bytes) {
                    Ok(pk) => col_pks.push(pk),
                    Err(_) => return (false, zero),
                }
            }
            pk_table.push(col_pks);
        }

        let sig = match self.codec.sig_decode(sig_bytes) {
            Ok(s) => s,
            Err(_) => return (false, zero),
        };

        self.verify_decoded_with_timings(pp, mu, &sig, &pk_table)
    }

    /// Internal verifier acting on decoded inputs.  Same `false` on
    /// every failure semantics as [`verify`].
    pub fn verify_decoded(
        &self,
        pp: &[u8],
        mu: &[u8],
        sig: &Signature,
        pk_table: &[Vec<Vec<Poly>>],
    ) -> bool {
        let mut t = VerifyTimings::default();
        self.verify_impl(pp, mu, sig, pk_table, &mut t).is_ok()
    }

    /// Same as [`LoTRS::verify_decoded`] but additionally returns a
    /// [`VerifyTimings`] breakdown.  The split is binary-proof checks
    /// (`verify_bin`) vs everything else (`verify_dualms`, including
    /// the closing FS hash).  Sub-times are populated even on a `false`
    /// return; on early failure, blocks past the failure point are 0.
    pub fn verify_decoded_with_timings(
        &self,
        pp: &[u8],
        mu: &[u8],
        sig: &Signature,
        pk_table: &[Vec<Vec<Poly>>],
    ) -> (bool, VerifyTimings) {
        let mut t = VerifyTimings::default();
        let ok = self.verify_impl(pp, mu, sig, pk_table, &mut t).is_ok();
        (ok, t)
    }

    fn verify_impl(
        &self,
        pp: &[u8],
        mu: &[u8],
        sig: &Signature,
        pk_table: &[Vec<Vec<Poly>>],
        timings: &mut VerifyTimings,
    ) -> core::result::Result<(), ()> {
        let par = &self.par;
        let rq = &self.r_q;
        let rqh = &self.r_qhat;

        // --- norm bounds on f1, z_b ------------------------------------
        // Binary-proof bucket.  Bound comparisons mirror Python's
        // `int > float` semantics so we stay on the float side.  f0
        // is bound-checked after reconstruction below.
        let bin_t0 = Instant::now();
        let f1_inf_bound = par.B_f1();
        for poly in &sig.pi.f1 {
            if (rq.inf_norm(poly) as f64) > f1_inf_bound {
                timings.verify_bin += bin_t0.elapsed();
                return Err(());
            }
        }
        let zb_inf_bound = 6.0 * par.phi_b * par.B_b();
        for poly in &sig.pi.z_b {
            if (rqh.inf_norm(poly) as f64) > zb_inf_bound {
                timings.verify_bin += bin_t0.elapsed();
                return Err(());
            }
        }
        timings.verify_bin += bin_t0.elapsed();

        // --- norm bounds on z_tilde, r_tilde, e_tilde (l2) -------------
        // DualMS bucket.  Use the triangle bound for e_u = z''_u + r''_u:
        // width sigma_0 + sigma_0'. This is intentionally looser than the
        // compact max-width expression and admits honest signatures with
        // comfortable slack.
        let dm_t0 = Instant::now();
        let t = par.tail_t;
        let b_z_sq = (par.sigma_0() * t).powi(2) * (par.T * par.d * par.l) as f64;
        let b_r_sq = (par.sigma_0_prime() * t).powi(2) * (par.T * par.d * par.l_prime) as f64;
        let b_e_sq =
            ((par.sigma_0() + par.sigma_0_prime()) * t).powi(2) * (par.T * par.d * par.k) as f64;
        if (rq.vec_l2_norm_sq(&sig.z_tilde) as f64) > b_z_sq {
            timings.verify_dualms += dm_t0.elapsed();
            return Err(());
        }
        if (rq.vec_l2_norm_sq(&sig.r_tilde) as f64) > b_r_sq {
            timings.verify_dualms += dm_t0.elapsed();
            return Err(());
        }
        if (rq.vec_l2_norm_sq(&sig.e_tilde) as f64) > b_e_sq {
            timings.verify_dualms += dm_t0.elapsed();
            return Err(());
        }
        timings.verify_dualms += dm_t0.elapsed();

        // --- reconstruct full f, f0 bound, g0/g1 bounds, A_hat_bin ----
        // All of this is the binary ring proof.
        let bin_t0 = Instant::now();
        let f_full = self.reconstruct_f(&sig.pi.x, &sig.pi.f1);

        // B_f0 = 1 + tau_{f_0} * sqrt(beta-1) * sigma_a.  (Fig. 6 line 4.)
        let f0_inf_bound = par.B_f0();
        for j in 0..par.kappa {
            if (rq.inf_norm(&f_full[j][0]) as f64) > f0_inf_bound {
                timings.verify_bin += bin_t0.elapsed();
                return Err(());
            }
        }

        let (g0, g1) = self.compute_g_qhat(&sig.pi.x, &f_full);
        if (rqh.vec_inf_norm(&g0) as f64) > par.B_g0() {
            timings.verify_bin += bin_t0.elapsed();
            return Err(());
        }
        if (rqh.vec_inf_norm(&g1) as f64) > par.B_g1() {
            timings.verify_bin += bin_t0.elapsed();
            return Err(());
        }

        // --- reconstruct A_hat_bin from binary proof -------------------
        let g_all = self.interleave_g(&g0, &g1);
        let f_all = self.flatten_f(&f_full);

        // f values are small — embed from signed into R_qhat
        let f_in_qhat: PolyVec = f_all
            .iter()
            .map(|c| rqh.from_centered(&rq.centered(c)))
            .collect();

        let g_mat = self.expand_G(pp);
        let g_mat_ntt = rqh.mat_to_ntt(&g_mat);
        // vec_for_G = concat(z_b, f_in_qhat, g_all)
        let mut vec_for_g: PolyVec = Vec::new();
        vec_for_g.extend_from_slice(&sig.pi.z_b);
        vec_for_g.extend_from_slice(&f_in_qhat);
        vec_for_g.extend_from_slice(&g_all);
        let gv = match g_mat_ntt
            .as_ref()
            .and_then(|m| rqh.mat_vec_ntt(m, &vec_for_g))
        {
            Some(v) => v,
            None => rqh.mat_vec(&g_mat, &vec_for_g),
        };

        // A_hat_bin[i] = Gv[i] - scale · (x_hat · B_bin_hi[i])
        let x_in_qhat = rqh.from_centered(&rq.centered(&sig.pi.x));
        let scale_val = 1i64 << par.K_B;
        let mut a_hat_bin: PolyVec = Vec::with_capacity(par.n_hat);
        for i in 0..par.n_hat {
            let prod = rqh.mul(&x_in_qhat, &sig.pi.b_bin_hi[i]);
            let scaled = rqh.scale(scale_val, &prod);
            a_hat_bin.push(rqh.sub(&gv[i], &scaled));
        }

        // --- decompose A_hat_bin; check |lo| bound ---------------------
        let mut a_hat_hi: PolyVec = Vec::with_capacity(par.n_hat);
        let lo_bound = (1i64 << (par.K_A - 1)) - (par.w as i64) * (1i64 << (par.K_B - 1));
        for comp in &a_hat_bin {
            let (hi, lo) = rqh.centered_decompose(comp, par.K_A);
            let lo_signed = rqh.centered(&lo);
            for c in &lo_signed {
                if c.unsigned_abs() as i64 > lo_bound {
                    timings.verify_bin += bin_t0.elapsed();
                    return Err(());
                }
            }
            a_hat_hi.push(hi);
        }
        timings.verify_bin += bin_t0.elapsed();

        // --- lhs = A_bar z̃ + B_bar r̃ + ẽ, pk_sum, w_hat_0, FS -------
        // DualMS bucket.  We attribute pk_sum (the `Σ_i (Π_j f) · t̃_i`
        // ring-keyed loop) to DualMS because it's part of the LHS=RHS
        // multi-signature closing check, not a separate ring-membership
        // step.  Same for the closing FS rehash.
        let dm_t0 = Instant::now();
        let pk_hash = self.pk_table_hash(pk_table);
        let b_mat = self.hash_com(&pk_hash, mu);
        let b_mat_ntt = rq.mat_to_ntt(&b_mat);
        let agg_keys = self.kagg_with_hash(pk_table, &pk_hash);
        let a_mat = self.expand_A(pp);
        let a_mat_ntt = rq.mat_to_ntt(&a_mat);

        // full_z = z_tilde || 0^k, full_r = r_tilde || 0^k
        let zero_k: PolyVec = rq.vec_zero(par.k);
        let mut full_z: PolyVec = sig.z_tilde.clone();
        full_z.extend(zero_k.iter().cloned());
        let mut full_r: PolyVec = sig.r_tilde.clone();
        full_r.extend(zero_k.iter().cloned());

        // A_bar · full_z = A·z_tilde + 0   (since right half of full_z is 0)
        // B_bar · full_r = B·r_tilde + 0
        let lhs_a = match a_mat_ntt
            .as_ref()
            .and_then(|m| rq.mat_vec_ntt(m, &sig.z_tilde))
        {
            Some(v) => v,
            None => rq.mat_vec(&a_mat, &sig.z_tilde),
        };
        let lhs_b = match b_mat_ntt
            .as_ref()
            .and_then(|m| rq.mat_vec_ntt(m, &sig.r_tilde))
        {
            Some(v) => v,
            None => rq.mat_vec(&b_mat, &sig.r_tilde),
        };
        let mut lhs = rq.vec_add(&lhs_a, &lhs_b);
        lhs = rq.vec_add(&lhs, &sig.e_tilde);

        // --- reconstruct w_hat_0 ---------------------------------------
        let mut x_pow = vec![rq.one()];
        for _ in 0..par.kappa {
            let last = x_pow.last().unwrap().clone();
            x_pow.push(rq.mul(&last, &sig.pi.x));
        }

        // pk_sum = sum_i (prod_j f_{j, i_j}) · t_tilde_i
        let mut pk_sum = rq.vec_zero(par.k);
        for i in 0..par.N() {
            let digits = Self::base_digits(i, par.beta, par.kappa);
            let mut selector = rq.one();
            for j in 0..par.kappa {
                selector = rq.mul(&selector, &f_full[j][digits[j]]);
            }
            let scaled = rq.vec_scale(&selector, &agg_keys[i]);
            pk_sum = rq.vec_add(&pk_sum, &scaled);
        }

        let mut w_hat_0 = rq.vec_sub(&pk_sum, &lhs);
        for j in 1..par.kappa {
            let scale_kw = 1i64 << par.K_w;
            let scaled_poly = rq.vec_scale(&x_pow[j], &sig.pi.w_tilde_hi[j - 1]);
            let term = rq.vec_scale_int(scale_kw, &scaled_poly);
            w_hat_0 = rq.vec_sub(&w_hat_0, &term);
        }

        // decompose w_hat_0 → w_hat_0^(1)
        let mut w_hat_0_hi = Vec::with_capacity(par.k);
        for comp in &w_hat_0 {
            let (hi, _) = rq.centered_decompose(comp, par.K_w);
            w_hat_0_hi.push(hi);
        }

        // --- Fiat-Shamir check -----------------------------------------
        let mut w_hi_for_hash: Vec<Vec<Poly>> = Vec::with_capacity(par.kappa);
        w_hi_for_hash.push(w_hat_0_hi);
        for j_vec in &sig.pi.w_tilde_hi {
            w_hi_for_hash.push(j_vec.clone());
        }

        let x_check = self.hash_fs(
            mu,
            &a_hat_hi,
            &sig.pi.b_bin_hi,
            &w_hi_for_hash,
            &pk_hash,
        );
        timings.verify_dualms += dm_t0.elapsed();

        if x_check == sig.pi.x {
            Ok(())
        } else {
            Err(())
        }
    }

    // =================================================================
    //  Signer support
    // =================================================================

    /// Expand the Gaussian mask coefficients `a_{j, i}` used by the
    /// binary proof.  For `i ∈ [1, beta)` they are freshly sampled;
    /// `a_{j, 0}` is set to `-sum_{i ≥ 1} a_{j, i}` so the full row
    /// sums to zero.  Returns signed coefficient lists (not reduced
    /// mod any prime) — matches `_expand_a_coeffs` in Python.
    pub fn expand_a_coeffs(&self, rho: &[u8], attempt: u32) -> Vec<Vec<Vec<i64>>> {
        let par = &self.par;
        let mut out: Vec<Vec<Vec<i64>>> = Vec::with_capacity(par.kappa);
        for j in 0..par.kappa {
            let mut row: Vec<Vec<i64>> = vec![Vec::new(); par.beta];
            // i >= 1: freshly sampled Gaussians.
            for i in 1..par.beta {
                let mut x = Xof::new(
                    rho,
                    &[
                        Tag::Bytes(b"a"),
                        Tag::Int(j as u32),
                        Tag::Int(i as u32),
                        Tag::Int(attempt),
                    ],
                );
                row[i] = xof_sample_gaussian(&mut x, self.cdt_sigma_a, par.lam, par.d);
            }
            // a_{j, 0} = -sum_{i >= 1} a_{j, i}
            let mut acc = vec![0i64; par.d];
            for i in 1..par.beta {
                for c in 0..par.d {
                    acc[c] += row[i][c];
                }
            }
            row[0] = acc.iter().map(|&v| -v).collect();
            out.push(row);
        }
        out
    }

    /// Selector polynomial coefficients `p_{i, j}` — at `kappa == 1`
    /// this degenerates to `p_{i, 0} = a_{0, i_0}` (see Remark 2
    /// in the paper).  Returns `p[i][0]` as `R_q` elements.
    pub fn compute_p_coeffs(&self, a_coeffs: &[Vec<Vec<i64>>], _ell: usize) -> Vec<Vec<Poly>> {
        let par = &self.par;
        assert_eq!(par.kappa, 1);
        let mut p: Vec<Vec<Poly>> = Vec::with_capacity(par.N());
        for i in 0..par.N() {
            // For kappa == 1, i decomposes as [i] in base-beta.  So
            //   p[i][0] = a_{0, i}.
            p.push(vec![self.r_q.from_centered(&a_coeffs[0][i])]);
        }
        p
    }

    /// `A_bar · y` where `A_bar = [A | I_k]` and `y` is length `l + k`.
    /// Equivalent to `A · y[..l] + y[l..l+k]` but kept as a dedicated
    /// helper so the sign1 flow reads straightforwardly.
    fn mat_vec_augmented(
        &self,
        ring: &Ring,
        a: &PolyMat,
        a_ntt: Option<&NttPolyMat>,
        y_full: &[Poly],
        left_cols: usize,
    ) -> PolyVec {
        let left = &y_full[..left_cols];
        let right = &y_full[left_cols..];
        let a_y = match a_ntt.and_then(|m| ring.mat_vec_ntt(m, left)) {
            Some(v) => v,
            None => ring.mat_vec(a, left),
        };
        a_y.iter()
            .zip(right.iter())
            .map(|(x, r)| ring.add(x, r))
            .collect()
    }

    // =================================================================
    //  Sign_1  (Fig. 4)
    // =================================================================

    /// Per-signer round 1.  Produces the commitments `w_{u, j}` and
    /// returns a state struct carrying everything `sign2` will need.
    ///
    /// `ctx` bundles every pp-/μ-/PK-dependent structure that is
    /// constant across all `T` signers and every restart attempt —
    /// `A`, `G`, `B`, the PK-table digest, and all `α_u`.  sign1
    /// no longer re-expands those on each call.
    fn sign1(
        &self,
        ctx: &SigningContext,
        sk_u: &[Poly],
        row_u: usize,
        ell: usize,
        pk_table: &[Vec<Vec<Poly>>],
        rho: &[u8],
        attempt: u32,
    ) -> (Sign1State, Vec<PolyVec>) {
        let par = &self.par;
        let rq = &self.r_q;
        let kappa = par.kappa;

        let alpha_u = &ctx.alphas[row_u];

        let a_coeffs = self.expand_a_coeffs(rho, attempt);
        let p_coeffs = self.compute_p_coeffs(&a_coeffs, ell);

        let a_mat = &ctx.a_mat;
        let b_mat = &ctx.b_mat;

        // Sample y, r per round
        let mut y_list: Vec<PolyVec> = Vec::with_capacity(kappa);
        let mut r_list: Vec<PolyVec> = Vec::with_capacity(kappa);
        for j in 0..kappa {
            let mut xof_y = Xof::new(
                rho,
                &[
                    Tag::Bytes(b"y"),
                    Tag::Int(row_u as u32),
                    Tag::Int(j as u32),
                    Tag::Int(attempt),
                ],
            );
            let mut y_j: PolyVec = Vec::with_capacity(par.l + par.k);
            for _ in 0..(par.l + par.k) {
                // sigma_0 mask width — CDT at TEST, FACCT at BENCH / PRODUCTION.
                let signed = match &self.sampler_0 {
                    MaskSampler::Cdt(cdt) => xof_sample_gaussian(&mut xof_y, cdt, par.lam, par.d),
                    MaskSampler::Facct(p) => xof_sample_gaussian_facct(&mut xof_y, p, par.d),
                };
                y_j.push(rq.from_centered(&signed));
            }
            y_list.push(y_j);

            let mut xof_r = Xof::new(
                rho,
                &[
                    Tag::Bytes(b"r"),
                    Tag::Int(row_u as u32),
                    Tag::Int(j as u32),
                    Tag::Int(attempt),
                ],
            );
            let mut r_j: PolyVec = Vec::with_capacity(par.l_prime + par.k);
            for _ in 0..(par.l_prime + par.k) {
                // sigma_0_prime mask width — same CDT-vs-FACCT split as sigma_0.
                let signed = match &self.sampler_0_prime {
                    MaskSampler::Cdt(cdt) => xof_sample_gaussian(&mut xof_r, cdt, par.lam, par.d),
                    MaskSampler::Facct(p) => xof_sample_gaussian_facct(&mut xof_r, p, par.d),
                };
                r_j.push(rq.from_centered(&signed));
            }
            r_list.push(r_j);
        }

        // Commitments w_{u, j}
        let mut coms: Vec<PolyVec> = Vec::with_capacity(kappa);
        for j in 0..kappa {
            let term_a =
                self.mat_vec_augmented(rq, a_mat, ctx.a_mat_ntt.as_ref(), &y_list[j], par.l);
            let term_b =
                self.mat_vec_augmented(rq, b_mat, ctx.b_mat_ntt.as_ref(), &r_list[j], par.l_prime);

            // pk_term = α_u · Σ_i (p_{i,j} · t_{u,i})
            // — accumulate in place rather than allocating a fresh PolyVec per i.
            let mut pk_term: PolyVec = rq.vec_zero(par.k);
            for i in 0..par.N() {
                rq.vec_add_scaled(&mut pk_term, &p_coeffs[i][j], &pk_table[i][row_u]);
            }
            // scale by α_u in place
            for poly in &mut pk_term {
                *poly = rq.mul(alpha_u, poly);
            }

            // w_uj = term_a + term_b + pk_term, accumulating in place
            let mut w_uj = term_a;
            rq.vec_add_assign(&mut w_uj, &term_b);
            rq.vec_add_assign(&mut w_uj, &pk_term);
            coms.push(w_uj);
        }

        let state = Sign1State {
            ell,
            sk_u: sk_u.to_vec(),
            y_list,
            r_list,
            rho: rho.to_vec(),
            attempt,
            row_u,
        };
        (state, coms)
    }

    // =================================================================
    //  Sign_2  (Fig. 5)
    // =================================================================

    /// Per-signer round 2.  Returns `None` if any rejection triggers a
    /// restart; otherwise a per-signer partial transcript.
    ///
    /// `sign_bin_elapsed` is incremented by the wall-clock time spent
    /// inside the [`LoTRS::sign_bin`] call so the caller can attribute
    /// it separately to the ring (binary-proof) bucket.
    fn sign2(
        &self,
        ctx: &SigningContext,
        state: &Sign1State,
        mu: &[u8],
        all_coms: &[Vec<PolyVec>],
        sign_bin_elapsed: &mut Duration,
    ) -> Option<Sign2Transcript> {
        let par = &self.par;
        let rq = &self.r_q;
        let kappa = par.kappa;

        // Aggregate commitments  w_tilde_j = sum_u w_{u', j}
        let mut w_tilde: Vec<PolyVec> = Vec::with_capacity(kappa);
        for j in 0..kappa {
            let mut acc = rq.vec_zero(par.k);
            for u_prime in 0..par.T {
                rq.vec_add_assign(&mut acc, &all_coms[u_prime][j]);
            }
            w_tilde.push(acc);
        }

        // Decompose w_tilde_j = 2^{K_w} w_tilde_j^(1) + w_tilde_j^(0)
        let mut w_tilde_hi: Vec<PolyVec> = Vec::with_capacity(kappa);
        let mut w_tilde_lo: Vec<PolyVec> = Vec::with_capacity(kappa);
        for j in 0..kappa {
            let (hi_j, lo_j): (PolyVec, PolyVec) = w_tilde[j]
                .iter()
                .map(|c| rq.centered_decompose(c, par.K_w))
                .unzip();
            w_tilde_hi.push(hi_j);
            w_tilde_lo.push(lo_j);
        }

        // Binary selection proof — timed separately so the caller can
        // bucket it under the ring (RS) cost rather than DualMS.
        let bin_t0 = Instant::now();
        let pi_opt = self.sign_bin(ctx, state.ell, mu, &w_tilde_hi, &state.rho, state.attempt);
        *sign_bin_elapsed += bin_t0.elapsed();
        let pi = pi_opt?;

        let x = pi.x.clone();

        // w_tilde_0 stability check
        //   M_w = sum_{j >= 1} ||x^j||_1 * 2^{K_{w,j}-1}
        //   (For kappa == 1 the sum is empty and M_w = 0.)
        let mut m_w: u128 = 0;
        if kappa > 1 {
            let mut x_pow = vec![rq.one()];
            for _ in 0..(kappa - 1) {
                let last = x_pow.last().unwrap().clone();
                x_pow.push(rq.mul(&last, &x));
            }
            for j in 1..kappa {
                m_w += rq.l1_norm(&x_pow[j]) * (1u128 << (par.K_w - 1));
            }
        }
        let threshold_w: u128 = (1u128 << (par.K_w - 1)).saturating_sub(m_w);
        if (rq.vec_inf_norm(&w_tilde_lo[0]) as u128) > threshold_w {
            return None;
        }

        // Main response z_u = x^kappa * alpha_u * s_u - sum x^j y_{u, j}
        let mut x_pow = vec![rq.one()];
        for _ in 0..kappa {
            let last = x_pow.last().unwrap().clone();
            x_pow.push(rq.mul(&last, &x));
        }
        let alpha_u = &ctx.alphas[state.row_u];
        let coeff = rq.mul(&x_pow[kappa], alpha_u);
        let mut z_u = rq.vec_scale(&coeff, &state.sk_u);
        for j in 0..kappa {
            let term = rq.vec_scale(&x_pow[j], &state.y_list[j]);
            z_u = rq.vec_sub(&z_u, &term);
        }

        // Shift v = x^kappa alpha_u s_u - sum_{j >= 1} x^j y_{u, j}
        // (Same as z_u but without the j == 0 subtraction.)
        let mut shift = rq.vec_scale(&coeff, &state.sk_u);
        for j in 1..kappa {
            let t = rq.vec_scale(&x_pow[j], &state.y_list[j]);
            shift = rq.vec_sub(&shift, &t);
        }

        // Rejection sampling on z_u
        let z_u_signed: Vec<i64> = flatten_centered(rq, &z_u);
        let shift_signed: Vec<i64> = flatten_centered(rq, &shift);
        let mut xof_rej = Xof::new(
            &state.rho,
            &[
                Tag::Bytes(b"rej_z"),
                Tag::Int(state.row_u as u32),
                Tag::Int(state.attempt),
            ],
        );
        if !rej(&mut xof_rej, &z_u_signed, &shift_signed, par.phi, par.B_0()) {
            return None;
        }

        // Auxiliary response r_u = - sum x^j r_{u, j}
        let mut r_u: PolyVec = rq.vec_zero(par.l_prime + par.k);
        for j in 0..kappa {
            let t = rq.vec_scale(&x_pow[j], &state.r_list[j]);
            r_u = rq.vec_sub(&r_u, &t);
        }

        Some(Sign2Transcript { pi, z_u, r_u })
    }

    // =================================================================
    //  Sign_bin  (Fig. 5, right column)
    // =================================================================

    fn sign_bin(
        &self,
        ctx: &SigningContext,
        ell: usize,
        mu: &[u8],
        w_tilde_hi: &[PolyVec],
        rho: &[u8],
        attempt: u32,
    ) -> Option<_Proof> {
        let par = &self.par;
        let rq = &self.r_q;
        let rqh = &self.r_qhat;
        let kappa = par.kappa;
        let beta = par.beta;
        let d = par.d;

        let g_mat = &ctx.g_mat;

        // one-hot b encoding of ell
        let ell_digits = Self::base_digits(ell, beta, kappa);
        let mut b_vecs: Vec<Vec<Poly>> = Vec::with_capacity(kappa);
        for j in 0..kappa {
            let mut row: Vec<Poly> = (0..beta).map(|_| rqh.zero()).collect();
            row[ell_digits[j]] = rqh.one();
            b_vecs.push(row);
        }

        let a_coeffs = self.expand_a_coeffs(rho, attempt);

        // a_flat, c_flat, d_flat over R_qhat
        let mut a_flat: PolyVec = Vec::with_capacity(kappa * beta);
        let mut c_flat: PolyVec = Vec::with_capacity(kappa * beta);
        let mut d_flat: PolyVec = Vec::with_capacity(kappa * beta);
        for j in 0..kappa {
            for i in 0..beta {
                let a_ji = rqh.from_centered(&a_coeffs[j][i]);
                let c_ji = if b_vecs[j][i] == rqh.one() {
                    rqh.neg(&a_ji)
                } else {
                    a_ji.clone()
                };
                let d_ji = rqh.neg(&rqh.mul(&a_ji, &a_ji));
                a_flat.push(a_ji);
                c_flat.push(c_ji);
                d_flat.push(d_ji);
            }
        }

        // r_b ∈ S_1^{n_hat + k_hat} (ternary)
        let mut xof_rb = Xof::new(rho, &[Tag::Bytes(b"rb"), Tag::Int(attempt)]);
        let mut r_b: PolyVec = Vec::with_capacity(par.n_hat + par.k_hat);
        for _ in 0..(par.n_hat + par.k_hat) {
            let signed = xof_sample_short(&mut xof_rb, 1, d);
            r_b.push(rqh.from_centered(&signed));
        }

        // r_a ∈ D_{sigma_b}^{n_hat + k_hat}
        let mut xof_ra = Xof::new(rho, &[Tag::Bytes(b"ra"), Tag::Int(attempt)]);
        let mut r_a: PolyVec = Vec::with_capacity(par.n_hat + par.k_hat);
        for _ in 0..(par.n_hat + par.k_hat) {
            let signed = xof_sample_gaussian(&mut xof_ra, self.cdt_sigma_b, par.lam, d);
            r_a.push(rqh.from_centered(&signed));
        }

        // b_flat for B_bin commitment
        let mut b_flat: PolyVec = Vec::with_capacity(kappa * beta);
        for j in 0..kappa {
            for i in 0..beta {
                b_flat.push(b_vecs[j][i].clone());
            }
        }

        // B_bin = G (r_b, b, c)^T
        let mut vec_commit: PolyVec = Vec::new();
        vec_commit.extend_from_slice(&r_b);
        vec_commit.extend_from_slice(&b_flat);
        vec_commit.extend_from_slice(&c_flat);
        let b_bin = match ctx
            .g_mat_ntt
            .as_ref()
            .and_then(|m| rqh.mat_vec_ntt(m, &vec_commit))
        {
            Some(v) => v,
            None => rqh.mat_vec(&g_mat, &vec_commit),
        };

        // Decompose B_bin
        let mut b_bin_hi: PolyVec = Vec::with_capacity(par.n_hat);
        let mut _b_bin_lo: PolyVec = Vec::with_capacity(par.n_hat);
        for comp in &b_bin {
            let (hi, lo) = rqh.centered_decompose(comp, par.K_B);
            b_bin_hi.push(hi);
            _b_bin_lo.push(lo);
        }

        // A_bin = G (r_a, a, d)^T
        let mut vec_a: PolyVec = Vec::new();
        vec_a.extend_from_slice(&r_a);
        vec_a.extend_from_slice(&a_flat);
        vec_a.extend_from_slice(&d_flat);
        let a_bin = match ctx
            .g_mat_ntt
            .as_ref()
            .and_then(|m| rqh.mat_vec_ntt(m, &vec_a))
        {
            Some(v) => v,
            None => rqh.mat_vec(&g_mat, &vec_a),
        };

        // Decompose A_bin, high part only
        let mut a_bin_hi: PolyVec = Vec::with_capacity(par.n_hat);
        for comp in &a_bin {
            let (hi, _lo) = rqh.centered_decompose(comp, par.K_A);
            a_bin_hi.push(hi);
        }

        // Fiat-Shamir challenge
        let x = self.hash_fs(mu, &a_bin_hi, &b_bin_hi, w_tilde_hi, &ctx.pk_hash);

        // z_b = r_a + x * r_b
        let x_hat = rqh.from_centered(&rq.centered(&x));
        let z_b = rqh.vec_add(&r_a, &rqh.vec_scale(&x_hat, &r_b));

        // RejOp on z_b
        let z_b_signed = flatten_centered(rqh, &z_b);
        let v_signed = flatten_centered(rqh, &rqh.vec_scale(&x_hat, &r_b));
        let mut xof_rej_b = Xof::new(rho, &[Tag::Bytes(b"rej_b"), Tag::Int(attempt)]);
        if !rej_op(&mut xof_rej_b, &z_b_signed, &v_signed, par.phi_b, par.B_b()) {
            return None;
        }

        // f_{j, i} = x * delta_{ell_j, i} + a_{j, i} (as signed lists)
        let x_signed = rq.centered(&x);
        let mut f_full: Vec<Vec<Vec<i64>>> = Vec::with_capacity(kappa);
        for j in 0..kappa {
            let mut row: Vec<Vec<i64>> = Vec::with_capacity(beta);
            for i in 0..beta {
                let mut f_ji = a_coeffs[j][i].clone();
                if i == ell_digits[j] {
                    for c in 0..d {
                        f_ji[c] += x_signed[c];
                    }
                }
                row.push(f_ji);
            }
            f_full.push(row);
        }

        // f1_signed for rejection sampling
        let mut f1_signed: Vec<i64> = Vec::new();
        for j in 0..kappa {
            for i in 1..beta {
                f1_signed.extend_from_slice(&f_full[j][i]);
            }
        }
        let mut b1_signed: Vec<i64> = Vec::new();
        for j in 0..kappa {
            for i in 1..beta {
                if i == ell_digits[j] {
                    b1_signed.extend_from_slice(&x_signed);
                } else {
                    b1_signed.extend(core::iter::repeat(0i64).take(d));
                }
            }
        }
        let mut xof_rej_a = Xof::new(rho, &[Tag::Bytes(b"rej_a"), Tag::Int(attempt)]);
        if !rej(&mut xof_rej_a, &f1_signed, &b1_signed, par.phi_a, par.B_a()) {
            return None;
        }

        // Infinity-norm restart checks on f_1 and f_0 (paper Fig. 5
        // Sign_bin line 12).  Violated only with probability eps_each =
        // eps_tot / 4 in an honest execution.
        let b_f1 = par.B_f1();
        if f1_signed.iter().any(|&c| (c.unsigned_abs() as f64) > b_f1) {
            return None;
        }
        let b_f0 = par.B_f0();
        for j in 0..kappa {
            if f_full[j][0]
                .iter()
                .any(|&c| (c.unsigned_abs() as f64) > b_f0)
            {
                return None;
            }
        }

        // g_{j, i} = f_{j, i} * (x - f_{j, i}) in R_qhat
        let mut g0_list: PolyVec = Vec::with_capacity(kappa);
        let mut g1_list: PolyVec = Vec::new();
        for j in 0..kappa {
            for i in 0..beta {
                let f_ji = &f_full[j][i];
                let diff: Vec<i64> = (0..d).map(|c| x_signed[c] - f_ji[c]).collect();
                let g_ji = poly_mul_signed(f_ji, &diff, d);
                let g_ji_reduced = rqh.from_centered(&g_ji);
                if i == 0 {
                    g0_list.push(g_ji_reduced);
                } else {
                    g1_list.push(g_ji_reduced);
                }
            }
        }
        if (rqh.vec_inf_norm(&g0_list) as f64) > par.B_g0() {
            return None;
        }
        if (rqh.vec_inf_norm(&g1_list) as f64) > par.B_g1() {
            return None;
        }

        // A_hat_bin low-bit stability check
        let mut f_all: PolyVec = Vec::with_capacity(kappa * beta);
        for j in 0..kappa {
            for i in 0..beta {
                f_all.push(rqh.from_centered(&f_full[j][i]));
            }
        }
        let mut g_all: PolyVec = Vec::with_capacity(kappa * beta);
        for j in 0..kappa {
            g_all.push(g0_list[j].clone());
            for i in 1..beta {
                let idx = j * (beta - 1) + (i - 1);
                g_all.push(g1_list[idx].clone());
            }
        }
        let mut vec_for_g: PolyVec = Vec::new();
        vec_for_g.extend_from_slice(&z_b);
        vec_for_g.extend_from_slice(&f_all);
        vec_for_g.extend_from_slice(&g_all);
        let gv = match ctx
            .g_mat_ntt
            .as_ref()
            .and_then(|m| rqh.mat_vec_ntt(m, &vec_for_g))
        {
            Some(v) => v,
            None => rqh.mat_vec(&g_mat, &vec_for_g),
        };

        let scale_kb = 1i64 << par.K_B;
        let mut a_hat_bin: PolyVec = Vec::with_capacity(par.n_hat);
        for i in 0..par.n_hat {
            let prod = rqh.mul(&x_hat, &b_bin_hi[i]);
            let scaled = rqh.scale(scale_kb, &prod);
            a_hat_bin.push(rqh.sub(&gv[i], &scaled));
        }
        let lo_bound = (1i64 << (par.K_A - 1)) - (par.w as i64) * (1i64 << (par.K_B - 1));
        for comp in &a_hat_bin {
            let (_hi, lo) = rqh.centered_decompose(comp, par.K_A);
            let lo_s = rqh.centered(&lo);
            for c in &lo_s {
                if c.unsigned_abs() as i64 > lo_bound {
                    return None;
                }
            }
        }

        // Final pi — f1 stored as R_q elements
        let mut f1_rq: PolyVec = Vec::with_capacity(kappa * (beta - 1));
        for j in 0..kappa {
            for i in 1..beta {
                f1_rq.push(rq.from_centered(&f_full[j][i]));
            }
        }

        // Paper (new): w_tilde_hi[0] is NOT included; verifier reconstructs.
        let w_tilde_hi_rest: Vec<PolyVec> = w_tilde_hi.iter().skip(1).cloned().collect();

        Some(_Proof {
            b_bin_hi,
            w_tilde_hi: w_tilde_hi_rest,
            x,
            f1: f1_rq,
            z_b,
        })
    }

    // =================================================================
    //  SAgg  (Fig. 4 — aggregation)
    // =================================================================

    fn sagg(&self, sigmas: &[Sign2Transcript]) -> Option<Signature> {
        let par = &self.par;
        let rq = &self.r_q;
        if sigmas.len() != par.T {
            return None;
        }

        // Consistency: all signers must agree on pi.
        let pi0 = &sigmas[0].pi;
        for s in &sigmas[1..] {
            if !proofs_equal(&s.pi, pi0) {
                return None;
            }
        }

        let mut z_tilde = rq.vec_zero(par.l);
        let mut r_tilde = rq.vec_zero(par.l_prime);
        let mut e_tilde = rq.vec_zero(par.k);

        for s in sigmas {
            let z_prime = &s.z_u[..par.l];
            let z_dprime = &s.z_u[par.l..];
            let r_prime = &s.r_u[..par.l_prime];
            let r_dprime = &s.r_u[par.l_prime..];
            z_tilde = rq.vec_add(&z_tilde, z_prime);
            r_tilde = rq.vec_add(&r_tilde, r_prime);
            e_tilde = rq.vec_add(&e_tilde, &rq.vec_add(z_dprime, r_dprime));
        }

        Some(Signature {
            pi: pi0.clone(),
            z_tilde,
            r_tilde,
            e_tilde,
        })
    }

    // =================================================================
    //  sign — convenience loop (Fig. 4 + Fig. 5 combined)
    // =================================================================

    /// Run the full interactive signing protocol to completion.
    /// Returns the encoded aggregated signature bytes on success, or an
    /// error if all `par.max_attempts` attempts rejected.
    ///
    /// The CDT / FACCT dispatch for each Gaussian width is resolved
    /// once in [`LoTRS::try_new`]; unsupported parameter sets are
    /// refused there rather than here.
    pub fn sign(
        &self,
        pp: &[u8],
        sks: &[Vec<Poly>],
        ell: usize,
        mu: &[u8],
        pk_table: &[Vec<Vec<Poly>>],
        signing_seed: &[u8],
    ) -> core::result::Result<Vec<u8>, &'static str> {
        self.sign_with_timings(pp, sks, ell, mu, pk_table, signing_seed)
            .map(|(bytes, _)| bytes)
    }

    /// Same as [`LoTRS::sign`] but additionally returns a [`SignTimings`]
    /// breakdown that splits sign1 / sign2-rest / sign_bin / sagg.
    /// Times are summed across every rejection-sampling attempt.
    pub fn sign_with_timings(
        &self,
        pp: &[u8],
        sks: &[Vec<Poly>],
        ell: usize,
        mu: &[u8],
        pk_table: &[Vec<Vec<Poly>>],
        signing_seed: &[u8],
    ) -> core::result::Result<(Vec<u8>, SignTimings), &'static str> {
        let par = &self.par;
        if sks.len() != par.T {
            return Err("wrong number of signing keys");
        }
        if signing_seed.len() != 32 {
            return Err("signing_seed must be 32 bytes");
        }
        if pk_table.len() != par.N() {
            return Err("pk_table has wrong N");
        }

        // Expand every pp/μ/PK-dependent structure once — reused across
        // every signer and every rejection-sampling attempt.
        let pk_hash = self.pk_table_hash(pk_table);
        let alphas: Vec<Poly> = (0..par.T)
            .map(|u| self.hash_agg(&pk_hash, u as u32))
            .collect();
        let a_mat = self.expand_A(pp);
        let g_mat = self.expand_G(pp);
        let b_mat = self.hash_com(&pk_hash, mu);
        let ctx = SigningContext {
            a_mat_ntt: self.r_q.mat_to_ntt(&a_mat),
            g_mat_ntt: self.r_qhat.mat_to_ntt(&g_mat),
            b_mat_ntt: self.r_q.mat_to_ntt(&b_mat),
            a_mat,
            g_mat,
            b_mat,
            pk_hash,
            alphas,
        };

        let mut timings = SignTimings::default();

        for attempt in 0..par.max_attempts {
            timings.attempts = timings.attempts.saturating_add(1);
            let mut rho_xof = Xof::new(signing_seed, &[Tag::Bytes(b"rho"), Tag::Int(attempt)]);
            let rho = rho_xof.read_array::<32>();

            // round 1
            let r1_t0 = Instant::now();
            let mut states: Vec<Sign1State> = Vec::with_capacity(par.T);
            let mut all_coms: Vec<Vec<PolyVec>> = Vec::with_capacity(par.T);
            for u in 0..par.T {
                let (st, com) = self.sign1(&ctx, &sks[u], u, ell, pk_table, &rho, attempt);
                states.push(st);
                all_coms.push(com);
            }
            timings.sign1 += r1_t0.elapsed();

            // round 2 — sign2 internally writes its sign_bin sub-time
            // into `bin_acc`; the rest of sign2's wall-clock lands in
            // `sign2_rest` after subtraction.
            let r2_t0 = Instant::now();
            let mut bin_acc = Duration::ZERO;
            let mut sigmas: Vec<Sign2Transcript> = Vec::with_capacity(par.T);
            let mut restart = false;
            for u in 0..par.T {
                match self.sign2(&ctx, &states[u], mu, &all_coms, &mut bin_acc) {
                    Some(s) => sigmas.push(s),
                    None => {
                        restart = true;
                        break;
                    }
                }
            }
            let r2_total = r2_t0.elapsed();
            timings.sign_bin += bin_acc;
            timings.sign2_rest += r2_total.saturating_sub(bin_acc);
            if restart {
                continue;
            }

            let agg_t0 = Instant::now();
            let sig = self.sagg(&sigmas).ok_or("sagg failed")?;
            timings.sagg += agg_t0.elapsed();

            let bytes = self
                .codec
                .sig_encode(&sig)
                .map_err(|_| "sig_encode failed")?;
            return Ok((bytes, timings));
        }

        Err("signing failed after max_attempts")
    }
}

// =====================================================================
//  Transcript types internal to the signer
// =====================================================================

/// Per-signer carry-state between sign1 and sign2.  All pp / μ / PK
/// derived values live on [`SigningContext`] and are passed
/// alongside this struct.
struct Sign1State {
    ell: usize,
    sk_u: Vec<Poly>,
    y_list: Vec<PolyVec>,
    r_list: Vec<PolyVec>,
    rho: Vec<u8>,
    attempt: u32,
    row_u: usize,
}

/// Pre-computed, pp / mu / pk_table dependent values used by every
/// signer and every restart attempt of a single [`LoTRS::sign`] call.
/// Expanding `A`, `G`, `B` and hashing the PK table is the single
/// largest fixed cost at BENCH / PRODUCTION — doing it once per sign
/// instead of T × attempts times is why a typical 8-attempt BENCH sign
/// drops from O(thousands) of SHAKE reads to O(hundreds) on this axis.
struct SigningContext {
    a_mat: PolyMat,
    g_mat: PolyMat,
    b_mat: PolyMat,
    a_mat_ntt: Option<NttPolyMat>,
    g_mat_ntt: Option<NttPolyMat>,
    b_mat_ntt: Option<NttPolyMat>,
    pk_hash: [u8; 32],
    /// `alphas[u] = H_agg(pk_hash, u)` for every row `u ∈ [0, T)`.
    alphas: Vec<Poly>,
}

struct Sign2Transcript {
    pi: _Proof,
    z_u: PolyVec,
    r_u: PolyVec,
}

// Polynomial equality including inner-vector equality.
fn proofs_equal(a: &_Proof, b: &_Proof) -> bool {
    a.x == b.x
        && a.b_bin_hi == b.b_bin_hi
        && a.w_tilde_hi == b.w_tilde_hi
        && a.f1 == b.f1
        && a.z_b == b.z_b
}

// Flatten a vector of unsigned polys into a single i64 list in centred form.
fn flatten_centered(ring: &Ring, v: &[Poly]) -> Vec<i64> {
    let mut out = Vec::with_capacity(v.len() * ring.d);
    for p in v {
        out.extend_from_slice(&ring.centered(p));
    }
    out
}

/// Which Gaussian width on [`LoTRSParams`] a CDT serves.  Used by the
/// CDT resolver so each shipped table is tied to a specific
/// (parameter-set, width) pair.
#[derive(Copy, Clone)]
enum CdtWidth {
    SigmaA,
    SigmaB,
    Sigma0,
    Sigma0Prime,
}

/// Resolve the shipped CDT that corresponds to `(par, width)`.  Returns
/// `None` for any parameter set we don't have a shipped table for —
/// including custom `LoTRSParams` whose sigmas might coincidentally
/// share a truncation length with a built-in table.  This replaces the
/// earlier length-only lookup, which could silently pick the wrong
/// distribution.
///
/// A defensive length sanity check runs on the chosen table so that
/// any future drift between the Python-generated tables and the
/// runtime-computed sigmas surfaces as `None` rather than as silently
/// wrong samples.
fn resolve_cdt(par: &LoTRSParams, width: CdtWidth) -> Option<&'static [u128]> {
    let chosen = match (par.name, width) {
        ("test-32", CdtWidth::SigmaA) => cdt::CDT_SIGMA_A_TEST,
        ("test-32", CdtWidth::SigmaB) => cdt::CDT_SIGMA_B_TEST,
        ("test-32", CdtWidth::Sigma0) => cdt::CDT_SIGMA_0_TEST,
        ("test-32", CdtWidth::Sigma0Prime) => cdt::CDT_SIGMA_0_PRIME_TEST,
        ("lotrs-bench-4of32", CdtWidth::SigmaA) => cdt::CDT_SIGMA_A_BENCH_4OF32,
        ("lotrs-bench-4of32", CdtWidth::SigmaB) => cdt::CDT_SIGMA_B_BENCH_PRODUCTION,
        ("lotrs-bench-16of32", CdtWidth::SigmaA) => cdt::CDT_SIGMA_A_BENCH,
        ("lotrs-bench-16of32", CdtWidth::SigmaB) => cdt::CDT_SIGMA_B_BENCH_PRODUCTION,
        ("lotrs-128", CdtWidth::SigmaA) => cdt::CDT_SIGMA_A_PRODUCTION,
        ("lotrs-128", CdtWidth::SigmaB) => cdt::CDT_SIGMA_B_BENCH_PRODUCTION,
        // mask widths at BENCH / PRODUCTION use FACCT, so no CDT here.
        _ => return None,
    };
    let sigma = match width {
        CdtWidth::SigmaA => par.sigma_a(),
        CdtWidth::SigmaB => par.sigma_b(),
        CdtWidth::Sigma0 => par.sigma_0(),
        CdtWidth::Sigma0Prime => par.sigma_0_prime(),
    };
    let expected_len = (14.0 * sigma).ceil() as usize + 1;
    if chosen.len() == expected_len {
        Some(chosen)
    } else {
        None
    }
}

/// Build the mask sampler declared on the parameter set.
fn build_mask_sampler(par: &LoTRSParams, width: CdtWidth) -> Option<MaskSampler> {
    let sigma = match width {
        CdtWidth::Sigma0 => par.sigma_0(),
        CdtWidth::Sigma0Prime => par.sigma_0_prime(),
        _ => unreachable!("mask sampler only applies to sigma_0 / sigma_0_prime"),
    };
    match par.mask_sampler {
        MaskSamplerKind::Cdt => resolve_cdt(par, width).map(MaskSampler::Cdt),
        MaskSamplerKind::Facct => prepare_facct(sigma).map(MaskSampler::Facct),
    }
}

// =====================================================================
//  Signed-integer negacyclic convolution  (Z, not Z_q).
//  Used for the quadratic proof terms  g = f · (x − f).
// =====================================================================

fn poly_mul_signed(a: &[i64], b: &[i64], d: usize) -> Vec<i64> {
    let mut c = vec![0i128; d];
    for i in 0..d {
        let ai = a[i] as i128;
        if ai == 0 {
            continue;
        }
        for j in 0..d {
            let prod = ai * b[j] as i128;
            let k = i + j;
            if k < d {
                c[k] += prod;
            } else {
                c[k - d] -= prod;
            }
        }
    }
    // For our parameter ranges this always fits in i64.
    c.into_iter().map(|v| v as i64).collect()
}

/// Panic-free top-level entry point.  Returns `false` on every invalid-
/// input class, *including* unsupported parameter sets (e.g. `kappa > 1`).
pub fn verify(
    par: LoTRSParams,
    pp: &[u8],
    mu: &[u8],
    sig_bytes: &[u8],
    pk_table_bytes: &[Vec<Vec<u8>>],
) -> bool {
    match LoTRS::try_new(par) {
        Some(scheme) => scheme.verify(pp, mu, sig_bytes, pk_table_bytes),
        None => false,
    }
}

// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::TEST_PARAMS;

    #[test]
    fn keygen_public_key_shape() {
        let scheme = LoTRS::new(TEST_PARAMS);
        let pp = scheme.setup(&[0u8; 32]);
        let (sk, pk) = scheme.keygen(&pp, &[1u8; 32]);
        assert_eq!(sk.len(), TEST_PARAMS.l + TEST_PARAMS.k);
        assert_eq!(pk.len(), TEST_PARAMS.k);
    }

    #[test]
    fn base_digits_are_little_endian() {
        let d = LoTRS::base_digits(5, 4, 2);
        assert_eq!(d, vec![1, 1]); // 5 in base 4, LE = [1, 1]
    }

    #[test]
    fn verify_rejects_wrong_pp_length() {
        let scheme = LoTRS::new(TEST_PARAMS);
        assert!(!scheme.verify(&[0u8; 5], b"hi", b"", &[]));
    }

    #[test]
    fn verify_rejects_empty_pk_table() {
        let scheme = LoTRS::new(TEST_PARAMS);
        assert!(!scheme.verify(&[0u8; 32], b"hi", b"", &[]));
    }

    #[test]
    fn top_level_verify_returns_false_for_unsupported_kappa() {
        // par.kappa != 1 must surface as `false`, not a panic.
        let mut bogus = TEST_PARAMS;
        bogus.kappa = 2;
        assert!(!super::verify(bogus, &[0u8; 32], b"x", b"", &[]));
        assert!(super::LoTRS::try_new(bogus).is_none());
    }

    #[test]
    fn sign_rejects_on_wrong_key_count() {
        let scheme = LoTRS::new(TEST_PARAMS);
        let r = scheme.sign(&[0u8; 32], &[], 0, b"hi", &[], &[0u8; 32]);
        assert!(r.is_err());
    }

    #[test]
    fn mask_sampler_dispatch_matches_sigma() {
        // TEST sigma_0 is small enough for a CDT; BENCH_4OF32 /
        // BENCH_PARAMS / PRODUCTION sigma_0 is multi-million and must
        // go to FACCT.
        use crate::params::{BENCH_4OF32, BENCH_PARAMS, PRODUCTION_PARAMS};
        let test_scheme = LoTRS::new(TEST_PARAMS);
        assert!(matches!(test_scheme.sampler_0, MaskSampler::Cdt(_)));
        assert!(matches!(test_scheme.sampler_0_prime, MaskSampler::Cdt(_)));
        for par in &[BENCH_4OF32, BENCH_PARAMS, PRODUCTION_PARAMS] {
            let sch = LoTRS::new(*par);
            assert!(matches!(sch.sampler_0, MaskSampler::Facct(_)));
            assert!(matches!(sch.sampler_0_prime, MaskSampler::Facct(_)));
        }
    }

    #[test]
    fn cdt_resolution_rejects_unknown_parameter_name() {
        // A custom parameter set copies the TEST sigmas but uses a
        // different name.  We have no shipped CDT for it, so the new
        // resolver must refuse rather than silently reuse TEST's
        // tables based on length matching.
        use crate::params::MaskSamplerKind;
        let mut custom = TEST_PARAMS;
        custom.name = "custom-does-not-exist";
        custom.mask_sampler = MaskSamplerKind::Cdt;
        assert!(LoTRS::try_new(custom).is_none());
    }

    #[test]
    fn cdt_resolution_rejects_tampered_sigma() {
        // If sigma is shifted out from under the table (e.g. via a
        // parameter tweak that forgot to regenerate the shipped CDT),
        // resolve_cdt's length sanity check catches the mismatch.
        let mut tampered = TEST_PARAMS;
        tampered.phi = TEST_PARAMS.phi * 2.0; // doubles sigma_0
        assert!(LoTRS::try_new(tampered).is_none());
    }
}

// Re-export the signature types for downstream users.
pub use crate::codec::Proof;
