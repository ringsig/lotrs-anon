//! Polynomial ring `R_q = Z_q[X]/(X^d + 1)` — port of `lotrs-py/ring.py`.
//!
//! Coefficients are unsigned `u64` in `[0, q)` (canonical form).  Signed
//! / centred form `[-q/2, q/2]` is used where noted (returned as `i64`).
//!
//! Two multiplication paths:
//!
//! * **Schoolbook** — O(d²) reference, always available.
//! * **CRT-NTT** — enabled when `d ∈ {128, 256}`, see [`crate::aux_ntt`].
//!
//! Polynomials are plain `Vec<u64>` (or `&[u64]` / `&mut [u64]` in
//! signatures).  Vectors of polynomials are `Vec<Vec<u64>>`; matrices
//! are row-major `Vec<Vec<Vec<u64>>>`.

use crate::aux_ntt::{CrtBackend, CrtNttMat};

pub type Poly = Vec<u64>;
pub type PolyVec = Vec<Poly>;
pub type PolyMat = Vec<PolyVec>;
pub type NttPolyMat = CrtNttMat;

/// Arithmetic context for `R_q`.  Cheap to construct when the NTT fast
/// path does not apply; O(d) memory + O(d) precompute otherwise.
#[derive(Clone)]
pub struct Ring {
    pub q: u64,
    pub d: usize,
    pub half_q: u64,
    backend: Option<CrtBackend>,
}

impl Ring {
    pub fn new(q: u64, d: usize) -> Self {
        Self {
            q,
            d,
            half_q: q / 2,
            backend: CrtBackend::new(q, d),
        }
    }

    // ---- element creation -----------------------------------------------

    #[inline]
    pub fn zero(&self) -> Poly {
        vec![0u64; self.d]
    }

    #[inline]
    pub fn one(&self) -> Poly {
        let mut r = vec![0u64; self.d];
        r[0] = 1;
        r
    }

    pub fn const_poly(&self, c: i64) -> Poly {
        let mut r = vec![0u64; self.d];
        r[0] = (c.rem_euclid(self.q as i64)) as u64;
        r
    }

    // ---- coefficient form conversions -----------------------------------

    pub fn reduce(&self, a: &[u64]) -> Poly {
        a.iter().map(|&c| c % self.q).collect()
    }

    /// Unsigned `[0, q)` → signed `[-q/2, q/2]`.
    pub fn centered(&self, a: &[u64]) -> Vec<i64> {
        let q = self.q as i64;
        let h = self.half_q;
        a.iter()
            .map(|&c| if c > h { c as i64 - q } else { c as i64 })
            .collect()
    }

    /// Signed → unsigned `[0, q)`.
    pub fn from_centered(&self, a: &[i64]) -> Poly {
        let q = self.q as i64;
        a.iter().map(|&c| c.rem_euclid(q) as u64).collect()
    }

    // ---- element-wise arithmetic ---------------------------------------

    #[inline]
    fn add_coeff(&self, a: u64, b: u64) -> u64 {
        let q = self.q;
        let s = a + b;
        let mask = 0u64.wrapping_sub((s >= q) as u64);
        s.wrapping_sub(q & mask)
    }

    #[inline]
    fn sub_coeff(&self, a: u64, b: u64) -> u64 {
        let q = self.q;
        let d = a.wrapping_sub(b);
        let mask = 0u64.wrapping_sub((a < b) as u64);
        d.wrapping_add(q & mask)
    }

    #[inline]
    fn neg_coeff(&self, c: u64) -> u64 {
        let q = self.q;
        let zmask = 0u64.wrapping_sub((c != 0) as u64);
        q.wrapping_sub(c) & zmask
    }

    pub fn add(&self, a: &[u64], b: &[u64]) -> Poly {
        (0..self.d).map(|i| self.add_coeff(a[i], b[i])).collect()
    }

    pub fn sub(&self, a: &[u64], b: &[u64]) -> Poly {
        (0..self.d).map(|i| self.sub_coeff(a[i], b[i])).collect()
    }

    pub fn neg(&self, a: &[u64]) -> Poly {
        a.iter().map(|&c| self.neg_coeff(c)).collect()
    }

    pub fn scale(&self, c: i64, a: &[u64]) -> Poly {
        let q = self.q;
        let c_mod = c.rem_euclid(q as i64) as u64;
        a.iter()
            .map(|&ai| ((c_mod as u128 * ai as u128) % q as u128) as u64)
            .collect()
    }

    pub fn mul(&self, a: &[u64], b: &[u64]) -> Poly {
        if let Some(ref bk) = self.backend {
            bk.mul(a, b)
        } else {
            self.mul_schoolbook(a, b)
        }
    }

    /// Schoolbook negacyclic convolution — always available, used by
    /// tests as the correctness oracle.
    pub fn mul_schoolbook(&self, a: &[u64], b: &[u64]) -> Poly {
        let d = self.d;
        let q = self.q as i128;
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
                    c[k] = (c[k] + prod) % q;
                } else {
                    c[k - d] = (c[k - d] - prod) % q;
                }
            }
        }
        c.into_iter().map(|v| v.rem_euclid(q) as u64).collect()
    }

    // ---- norms (on centred representation) ------------------------------

    pub fn inf_norm(&self, a: &[u64]) -> u64 {
        let h = self.half_q;
        let q = self.q;
        a.iter()
            .map(|&c| if c > h { q - c } else { c })
            .max()
            .unwrap_or(0)
    }

    pub fn l2_norm_sq(&self, a: &[u64]) -> u128 {
        let h = self.half_q;
        let q = self.q;
        a.iter()
            .map(|&c| {
                let cc = if c > h { (q - c) as u128 } else { c as u128 };
                cc * cc
            })
            .sum()
    }

    pub fn l1_norm(&self, a: &[u64]) -> u128 {
        let h = self.half_q;
        let q = self.q;
        a.iter()
            .map(|&c| if c > h { (q - c) as u128 } else { c as u128 })
            .sum()
    }

    // ---- vector operations ---------------------------------------------

    pub fn vec_zero(&self, n: usize) -> PolyVec {
        (0..n).map(|_| self.zero()).collect()
    }

    pub fn vec_add(&self, u: &[Poly], v: &[Poly]) -> PolyVec {
        u.iter()
            .zip(v.iter())
            .map(|(a, b)| self.add(a, b))
            .collect()
    }

    pub fn vec_sub(&self, u: &[Poly], v: &[Poly]) -> PolyVec {
        u.iter()
            .zip(v.iter())
            .map(|(a, b)| self.sub(a, b))
            .collect()
    }

    pub fn vec_neg(&self, v: &[Poly]) -> PolyVec {
        v.iter().map(|p| self.neg(p)).collect()
    }

    pub fn vec_scale(&self, c_poly: &[u64], v: &[Poly]) -> PolyVec {
        v.iter().map(|p| self.mul(c_poly, p)).collect()
    }

    pub fn vec_scale_int(&self, c: i64, v: &[Poly]) -> PolyVec {
        v.iter().map(|p| self.scale(c, p)).collect()
    }

    // ---- in-place / fused helpers --------------------------------------
    // These exist to cut allocation churn on the sign / verify hot paths.
    // The generic public vec_* functions stay (used widely), but the
    // BENCH-dominating loops in sign1 / kagg switch to in-place variants.

    /// `a += b` in place.
    #[inline]
    pub fn add_assign(&self, a: &mut [u64], b: &[u64]) {
        for i in 0..self.d {
            a[i] = self.add_coeff(a[i], b[i]);
        }
    }

    /// `a -= b` in place.
    #[inline]
    pub fn sub_assign(&self, a: &mut [u64], b: &[u64]) {
        for i in 0..self.d {
            a[i] = self.sub_coeff(a[i], b[i]);
        }
    }

    /// `u += v` element-wise in place.
    pub fn vec_add_assign(&self, u: &mut [Poly], v: &[Poly]) {
        for (a, b) in u.iter_mut().zip(v.iter()) {
            self.add_assign(a, b);
        }
    }

    /// `u -= v` element-wise in place.
    pub fn vec_sub_assign(&self, u: &mut [Poly], v: &[Poly]) {
        for (a, b) in u.iter_mut().zip(v.iter()) {
            self.sub_assign(a, b);
        }
    }

    /// `acc += c · v` element-wise (ring multiply then add).  Fused helper
    /// used by `kagg` and the `pk_term` loop inside `sign1`.
    pub fn vec_add_scaled(&self, acc: &mut [Poly], c: &[u64], v: &[Poly]) {
        for i in 0..v.len() {
            let prod = self.mul(c, &v[i]);
            self.add_assign(&mut acc[i], &prod);
        }
    }

    /// `acc -= c · v` element-wise.
    pub fn vec_sub_scaled(&self, acc: &mut [Poly], c: &[u64], v: &[Poly]) {
        for i in 0..v.len() {
            let prod = self.mul(c, &v[i]);
            self.sub_assign(&mut acc[i], &prod);
        }
    }

    pub fn inner(&self, u: &[Poly], v: &[Poly]) -> Poly {
        let mut acc = self.zero();
        for i in 0..u.len() {
            let prod = self.mul(&u[i], &v[i]);
            self.add_assign(&mut acc, &prod);
        }
        acc
    }

    pub fn vec_inf_norm(&self, v: &[Poly]) -> u64 {
        v.iter().map(|p| self.inf_norm(p)).max().unwrap_or(0)
    }

    pub fn vec_l2_norm_sq(&self, v: &[Poly]) -> u128 {
        v.iter().map(|p| self.l2_norm_sq(p)).sum()
    }

    pub fn mat_vec(&self, m: &[PolyVec], v: &[Poly]) -> PolyVec {
        m.iter().map(|row| self.inner(row, v)).collect()
    }

    /// Pre-transform a matrix for repeated CRT-NTT matrix-vector
    /// products.  Returns `None` when the ring is using schoolbook
    /// multiplication.
    pub fn mat_to_ntt(&self, m: &[PolyVec]) -> Option<NttPolyMat> {
        self.backend.as_ref().map(|bk| bk.mat_to_ntt(m))
    }

    /// Matrix-vector product using a pre-transformed matrix.  Returns
    /// `None` when the ring has no CRT backend.
    pub fn mat_vec_ntt(&self, m: &NttPolyMat, v: &[Poly]) -> Option<PolyVec> {
        self.backend.as_ref().map(|bk| bk.mat_vec_ntt(m, v))
    }

    pub fn vec_concat(&self, parts: &[&[Poly]]) -> PolyVec {
        let mut out = Vec::new();
        for p in parts {
            out.extend_from_slice(p);
        }
        out
    }

    // ---- centred decomposition (Bai-Galbraith) --------------------------

    /// Split each centred coefficient `c` as `hi * 2^K + lo` with
    /// `lo ∈ (-2^{K-1}, 2^{K-1}]`.  Returns `(hi, lo)` both in unsigned
    /// canonical form.
    ///
    /// The split matches `lotrs-py/ring.py::Ring.centered_decompose`
    /// exactly — in particular, when `c ≡ 2^{K-1} (mod 2^K)` the low
    /// part is `+2^{K-1}` (not `-2^{K-1}`), which affects the value of
    /// `hi` and therefore the Fiat-Shamir hash.
    pub fn centered_decompose(&self, a: &[u64], k_bits: u32) -> (Poly, Poly) {
        let q = self.q as i64;
        let h = self.half_q as i64;
        let step = 1i64 << k_bits;
        let half_step = 1i64 << (k_bits - 1);

        let mut hi = vec![0u64; self.d];
        let mut lo = vec![0u64; self.d];

        for i in 0..self.d {
            let c = a[i] as i64;
            let c_centred = if c > h { c - q } else { c };
            // Python:  low = c % power; if low > half_power: low -= power
            let mut low = c_centred.rem_euclid(step);
            if low > half_step {
                low -= step;
            }
            let high = (c_centred - low) / step;
            lo[i] = low.rem_euclid(q) as u64;
            hi[i] = high.rem_euclid(q) as u64;
        }
        (hi, lo)
    }
}

// -------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{PRODUCTION_PARAMS, TEST_PARAMS};

    fn lcg(mut s: u64) -> impl FnMut(u64) -> u64 {
        move |bound: u64| {
            s = s
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (s >> 16) % bound
        }
    }

    #[test]
    fn commutativity_d32() {
        let R = Ring::new(TEST_PARAMS.q, 32);
        let mut r = lcg(1);
        let a: Poly = (0..32).map(|_| r(R.q)).collect();
        let b: Poly = (0..32).map(|_| r(R.q)).collect();
        assert_eq!(R.mul(&a, &b), R.mul(&b, &a));
    }

    #[test]
    fn multiplicative_identity() {
        let R = Ring::new(TEST_PARAMS.q, 32);
        let mut r = lcg(2);
        let a: Poly = (0..32).map(|_| r(R.q)).collect();
        assert_eq!(R.mul(&a, &R.one()), a);
    }

    #[test]
    fn x_to_the_d_is_minus_one() {
        let R = Ring::new(TEST_PARAMS.q, 32);
        let mut x = vec![0u64; 32];
        x[1] = 1;
        let mut prod = R.one();
        for _ in 0..32 {
            prod = R.mul(&prod, &x);
        }
        let mut expected = vec![0u64; 32];
        expected[0] = R.q - 1;
        assert_eq!(prod, expected);
    }

    #[test]
    fn ntt_agrees_with_schoolbook_d128() {
        let R = Ring::new(PRODUCTION_PARAMS.q, 128);
        let mut r = lcg(3);
        for _ in 0..5 {
            let a: Poly = (0..128).map(|_| r(R.q)).collect();
            let b: Poly = (0..128).map(|_| r(R.q)).collect();
            assert_eq!(R.mul(&a, &b), R.mul_schoolbook(&a, &b));
        }
    }

    #[test]
    fn ntt_mat_vec_agrees_with_plain_mat_vec() {
        let R = Ring::new(PRODUCTION_PARAMS.q, 128);
        let mut r = lcg(5);
        let mat: PolyMat = (0..3)
            .map(|_| (0..4).map(|_| (0..128).map(|_| r(R.q)).collect()).collect())
            .collect();
        let v: PolyVec = (0..4).map(|_| (0..128).map(|_| r(R.q)).collect()).collect();
        let mat_ntt = R.mat_to_ntt(&mat).expect("production ring has CRT backend");
        assert_eq!(R.mat_vec_ntt(&mat_ntt, &v).unwrap(), R.mat_vec(&mat, &v));
    }

    #[test]
    fn centered_decompose_roundtrip() {
        let R = Ring::new(PRODUCTION_PARAMS.q, 128);
        let mut r = lcg(4);
        for k_bits in [4u32, 8, 12] {
            let a: Poly = (0..128).map(|_| r(R.q)).collect();
            let (hi, lo) = R.centered_decompose(&a, k_bits);
            let reconstructed = R.add(&R.scale(1i64 << k_bits, &hi), &lo);
            assert_eq!(reconstructed, a);
            assert!(R.inf_norm(&lo) <= 1u64 << (k_bits - 1));
        }
    }
}
