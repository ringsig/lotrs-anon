//! CRT-NTT backend — port of `lotrs-py/aux_ntt.py`.
//!
//! Fast negacyclic multiplication in `R_q = Z_q[X]/(X^d + 1)` via two
//! auxiliary NTT-friendly primes near `2^48`.  Both primes are `≡ 1 mod
//! 512`, so each supports a standard radix-2 negacyclic NTT for `d ∈
//! {128, 256}`.  The scheme modulus `q` does **not** need to be
//! NTT-friendly — it only needs to satisfy the paper-facing conditions
//! (primality, `q ≡ 5 mod 8`, size).
//!
//! Pipeline for one multiplication:
//!
//! 1. Centre-lift each coefficient of `a, b` from `[0, q)` to `(-q/2, q/2]`.
//! 2. Reduce mod each auxiliary prime.
//! 3. Forward negacyclic NTT under each prime.
//! 4. Pointwise multiply.
//! 5. Inverse negacyclic NTT.
//! 6. Garner-CRT to the centred integer in `(-P/2, P/2]`.
//! 7. Reduce mod `q`.
//!
//! The auxiliary primes have the form `2^48 - c` with `c ∈ {16383, 19967}`.
//! The Rust path uses `u128` products and the Mersenne-style identity
//! `x mod p ≡ hi·c + lo (mod p)` for 96-bit `x = hi·2^48 + lo`; final
//! reductions in the transform hot path are mask-style conditional
//! subtracts.  The backend also exposes pre-transformed matrix-vector
//! products so public matrices are transformed once and reused.

// ---- constants -----------------------------------------------------------

/// First auxiliary prime: `2^48 - 16383`.
pub const AUX_P1: u64 = 281_474_976_694_273;
/// Second auxiliary prime: `2^48 - 19967`.
pub const AUX_P2: u64 = 281_474_976_690_689;

pub const AUX_MERSENNE_C1: u64 = 16_383;
pub const AUX_MERSENNE_C2: u64 = 19_967;

/// Supported ring dimensions for the NTT fast path.
pub const SUPPORTED_D: &[usize] = &[128, 256];

// Product of the two auxiliary primes as a `u128`.
const AUX_P: u128 = AUX_P1 as u128 * AUX_P2 as u128;

// ---- compile-time specialisation on the auxiliary prime -----------------
//
// The butterfly inner loop is the single hottest piece of code in the
// crate.  Making `P` (and the Mersenne constant `C = 2^48 − P`) compile-
// time constants lets the compiler specialise `reduce_pseudo_mersenne`
// into a pair of constant multiplies, fold all `% p` into conditional
// subtracts, and unroll / vectorise the butterfly.
//
// Concrete implementors: [`P1Marker`] and [`P2Marker`].

/// Compile-time parameters for a specialised auxiliary NTT.
pub trait AuxPrime {
    const P: u64;
    const C: u64;
}

/// Marker type representing `AUX_P1`.
#[derive(Clone, Copy)]
pub enum P1Marker {}
impl AuxPrime for P1Marker {
    const P: u64 = AUX_P1;
    const C: u64 = AUX_MERSENNE_C1;
}

/// Marker type representing `AUX_P2`.
#[derive(Clone, Copy)]
pub enum P2Marker {}
impl AuxPrime for P2Marker {
    const P: u64 = AUX_P2;
    const C: u64 = AUX_MERSENNE_C2;
}

/// Monomorphic `mul_mod` — equivalent to `mul_mod(a, b, T::P)` but with
/// the Mersenne constants folded at compile time.
#[inline(always)]
fn mul_mod_aux<T: AuxPrime>(a: u64, b: u64) -> u64 {
    debug_assert!(a < T::P && b < T::P);
    reduce_pseudo_mersenne(a as u128 * b as u128, T::C, T::P)
}

#[inline(always)]
fn add_mod_aux<T: AuxPrime>(a: u64, b: u64) -> u64 {
    let s = a + b;
    let mask = 0u64.wrapping_sub((s >= T::P) as u64);
    s.wrapping_sub(T::P & mask)
}

#[inline(always)]
fn sub_mod_aux<T: AuxPrime>(a: u64, b: u64) -> u64 {
    let d = a.wrapping_sub(b);
    let mask = 0u64.wrapping_sub((a < b) as u64);
    d.wrapping_add(T::P & mask)
}

// ---- primitive root search ----------------------------------------------

/// Return the smallest `x ∈ [2, p)` with `x^d ≡ -1 (mod p)`, i.e. exact
/// multiplicative order `2d`.  Requires `2d | p - 1`.
pub fn find_primitive_2d_root(p: u64, d: usize) -> Option<u64> {
    let two_d = (2 * d) as u64;
    if (p - 1) % two_d != 0 {
        return None;
    }
    let exponent = (p - 1) / two_d;
    for x in 2..p {
        let psi = pow_mod(x, exponent, p);
        if pow_mod(psi, d as u64, p) == p - 1 {
            return Some(psi);
        }
    }
    None
}

/// `base^exp mod p` using 128-bit multiplies.  Constant-time in `exp`
/// only up to the bit-length; not meant for secret exponents.
#[inline]
pub fn pow_mod(mut base: u64, mut exp: u64, p: u64) -> u64 {
    let mut acc: u64 = 1;
    base %= p;
    while exp > 0 {
        if exp & 1 == 1 {
            acc = mul_mod(acc, base, p);
        }
        exp >>= 1;
        if exp > 0 {
            base = mul_mod(base, base, p);
        }
    }
    acc
}

/// Modular inverse via Fermat (p is prime).
#[inline]
pub fn inv_mod(a: u64, p: u64) -> u64 {
    pow_mod(a, p - 2, p)
}

/// `(a * b) mod p` via `u128`.  Both `a` and `b` assumed `< p < 2^63`.
#[inline]
pub fn mul_mod(a: u64, b: u64, p: u64) -> u64 {
    debug_assert!(a < p && b < p);
    match p {
        AUX_P1 => reduce_pseudo_mersenne(a as u128 * b as u128, AUX_MERSENNE_C1, AUX_P1),
        AUX_P2 => reduce_pseudo_mersenne(a as u128 * b as u128, AUX_MERSENNE_C2, AUX_P2),
        _ => ((a as u128 * b as u128) % p as u128) as u64,
    }
}

#[inline]
fn add_mod(a: u64, b: u64, p: u64) -> u64 {
    let s = a + b;
    let mask = 0u64.wrapping_sub((s >= p) as u64);
    s.wrapping_sub(p & mask)
}

#[inline]
fn sub_mod(a: u64, b: u64, p: u64) -> u64 {
    let d = a.wrapping_sub(b);
    let mask = 0u64.wrapping_sub((a < b) as u64);
    d.wrapping_add(p & mask)
}

#[inline]
fn reduce_pseudo_mersenne(x: u128, c: u64, p: u64) -> u64 {
    const MASK48: u128 = (1u128 << 48) - 1;
    let c128 = c as u128;
    let mut r = (x & MASK48) + c128 * (x >> 48);
    r = (r & MASK48) + c128 * (r >> 48);
    let mut out = r as u64;
    let mask = 0u64.wrapping_sub((out >= p) as u64);
    out = out.wrapping_sub(p & mask);
    out
}

// ---- bit-reversed twiddle table (FIPS 204 convention) ------------------

#[inline]
fn bitrev(mut k: usize, bits: u32) -> usize {
    let mut r = 0usize;
    for _ in 0..bits {
        r = (r << 1) | (k & 1);
        k >>= 1;
    }
    r
}

/// `zetas[k] = psi^bitrev(k, log2 d) mod p` for `k ∈ [0, d)`.
pub fn make_zeta_table(p: u64, d: usize, psi: u64) -> Vec<u64> {
    let bits = d.trailing_zeros();
    (0..d)
        .map(|k| pow_mod(psi, bitrev(k, bits) as u64, p))
        .collect()
}

// ---- radix-2 negacyclic NTT (FIPS 204 Algs. 41 / 42) -------------------

/// Forward negacyclic NTT.  Input in natural order; output in bit-reversed
/// order; values in `[0, p)`.
pub fn ntt(f: &mut [u64], zetas: &[u64], p: u64) {
    let d = f.len();
    let mut m = 0usize;
    let mut le = d / 2;
    while le >= 1 {
        let mut st = 0usize;
        while st < d {
            m += 1;
            let z = zetas[m];
            for j in st..st + le {
                let t = mul_mod(z, f[j + le], p);
                let fj = f[j];
                f[j + le] = sub_mod(fj, t, p);
                f[j] = add_mod(fj, t, p);
            }
            st += 2 * le;
        }
        le /= 2;
    }
}

/// Inverse negacyclic NTT.  Input bit-reversed; output natural order.
/// Includes the trailing `1/d` scale.
pub fn intt(f: &mut [u64], zetas: &[u64], p: u64) {
    let d = f.len();
    let mut m = d;
    let mut le = 1usize;
    while le < d {
        let mut st = 0usize;
        while st < d {
            m -= 1;
            let z = (p - zetas[m]) % p; // -zetas[m] mod p
            for j in st..st + le {
                let t = f[j];
                let u = f[j + le];
                f[j] = add_mod(t, u, p);
                f[j + le] = mul_mod(z, sub_mod(t, u, p), p);
            }
            st += 2 * le;
        }
        le *= 2;
    }
    let inv_d = inv_mod(d as u64, p);
    for x in f.iter_mut() {
        *x = mul_mod(*x, inv_d, p);
    }
}

// ---- monomorphic NTT butterflies ----------------------------------------

/// Forward negacyclic NTT, specialised on `T`.  Equivalent to
/// `ntt(f, zetas, T::P)` but with every `mul_mod` / `add_mod` /
/// `sub_mod` folded at compile time.
pub fn ntt_aux<T: AuxPrime>(f: &mut [u64], zetas: &[u64]) {
    let d = f.len();
    let mut m = 0usize;
    let mut le = d / 2;
    while le >= 1 {
        let mut st = 0usize;
        while st < d {
            m += 1;
            let z = zetas[m];
            for j in st..st + le {
                let t = mul_mod_aux::<T>(z, f[j + le]);
                let fj = f[j];
                f[j + le] = sub_mod_aux::<T>(fj, t);
                f[j] = add_mod_aux::<T>(fj, t);
            }
            st += 2 * le;
        }
        le /= 2;
    }
}

/// Inverse negacyclic NTT, specialised on `T`.  Includes the trailing
/// `1/d` scale.
pub fn intt_aux<T: AuxPrime>(f: &mut [u64], zetas: &[u64]) {
    let d = f.len();
    let mut m = d;
    let mut le = 1usize;
    while le < d {
        let mut st = 0usize;
        while st < d {
            m -= 1;
            // -zetas[m] mod P.  zetas[m] is always in (0, T::P) for a
            // valid twiddle table, but we guard the 0 case defensively.
            let z = if zetas[m] == 0 { 0 } else { T::P - zetas[m] };
            for j in st..st + le {
                let t = f[j];
                let u = f[j + le];
                f[j] = add_mod_aux::<T>(t, u);
                f[j + le] = mul_mod_aux::<T>(z, sub_mod_aux::<T>(t, u));
            }
            st += 2 * le;
        }
        le *= 2;
    }
    let inv_d = inv_mod(d as u64, T::P);
    for x in f.iter_mut() {
        *x = mul_mod_aux::<T>(*x, inv_d);
    }
}

// ---- per-prime context --------------------------------------------------

#[derive(Clone)]
struct AuxContext<T: AuxPrime> {
    zetas: Vec<u64>,
    _phantom: core::marker::PhantomData<fn() -> T>,
}

impl<T: AuxPrime> AuxContext<T> {
    fn new(d: usize) -> Self {
        let psi = find_primitive_2d_root(T::P, d).expect("aux prime does not support this d");
        Self {
            zetas: make_zeta_table(T::P, d, psi),
            _phantom: core::marker::PhantomData,
        }
    }

    fn forward(&self, input: &[u64]) -> Vec<u64> {
        let mut out = input.to_vec();
        ntt_aux::<T>(&mut out, &self.zetas);
        out
    }

    fn inverse_in_place(&self, f: &mut [u64]) {
        intt_aux::<T>(f, &self.zetas);
    }
}

/// One polynomial represented in the two auxiliary NTT domains.
#[derive(Clone)]
pub struct CrtNttPoly {
    c1: Vec<u64>,
    c2: Vec<u64>,
}

/// Matrix of polynomials represented in the two auxiliary NTT domains.
pub type CrtNttMat = Vec<Vec<CrtNttPoly>>;

// ---- CRT backend --------------------------------------------------------

/// Exact negacyclic multiplication in `R_q = Z_q[X]/(X^d + 1)`.
#[derive(Clone)]
pub struct CrtBackend {
    q: u64,
    d: usize,
    half_q: u64,
    ctx1: AuxContext<P1Marker>,
    ctx2: AuxContext<P2Marker>,
    p1_inv_mod_p2: u64,
}

impl CrtBackend {
    /// Build a backend for `R_q` at dimension `d`.
    ///
    /// Returns `None` if
    ///
    /// * `d` is not in [`SUPPORTED_D`], or
    /// * the CRT-reconstruction bound  `2 · d · ((q-1)/2)^2  <  P`
    ///   (where `P = p1 · p2 ≈ 2^96`) is violated — in that case CRT
    ///   cannot recover the exact integer product coefficient, so the
    ///   backend would silently produce wrong results for worst-case
    ///   operands.  The caller should fall back to schoolbook.
    pub fn new(q: u64, d: usize) -> Option<Self> {
        if !SUPPORTED_D.contains(&d) {
            return None;
        }
        if AUX_P1 <= q || AUX_P2 <= q {
            return None;
        }
        // Worst-case |c_k| for a negacyclic product of centred operands
        // in [-(q-1)/2, (q-1)/2]:
        //     |c_k|  <=  d * ((q-1)/2)^2
        // For exact CRT reconstruction we need  2 * max|c_k|  <  P.
        let half_q: u128 = ((q - 1) / 2) as u128;
        let half_q_sq = half_q * half_q; // always fits in u128
        match (d as u128).checked_mul(half_q_sq) {
            Some(bound) if bound < AUX_P / 2 => {}
            _ => return None,
        }
        Some(Self {
            q,
            d,
            half_q: q / 2,
            ctx1: AuxContext::<P1Marker>::new(d),
            ctx2: AuxContext::<P2Marker>::new(d),
            p1_inv_mod_p2: inv_mod(AUX_P1, AUX_P2),
        })
    }

    /// `a * b` mod `X^d + 1`, mod `q`.  Inputs and output in `[0, q)`.
    pub fn mul(&self, a: &[u64], b: &[u64]) -> Vec<u64> {
        assert_eq!(a.len(), self.d);
        assert_eq!(b.len(), self.d);

        // centre-lift and reduce mod each auxiliary prime
        let (a1, a2) = self.split(a);
        let (b1, b2) = self.split(b);

        let mut c1 =
            mul_pointwise_aux::<P1Marker>(&self.ctx1.forward(&a1), &self.ctx1.forward(&b1));
        self.ctx1.inverse_in_place(&mut c1);

        let mut c2 =
            mul_pointwise_aux::<P2Marker>(&self.ctx2.forward(&a2), &self.ctx2.forward(&b2));
        self.ctx2.inverse_in_place(&mut c2);

        self.crt_combine(&c1, &c2)
    }

    /// Transform a canonical `[0,q)` polynomial into the auxiliary NTT
    /// domains.  This is useful when the operand is reused across many
    /// multiplications, e.g. public matrices in `mat_vec`.
    pub fn to_ntt(&self, a: &[u64]) -> CrtNttPoly {
        assert_eq!(a.len(), self.d);
        let (a1, a2) = self.split(a);
        CrtNttPoly {
            c1: self.ctx1.forward(&a1),
            c2: self.ctx2.forward(&a2),
        }
    }

    /// Transform a matrix of canonical polynomials into the auxiliary NTT
    /// domains.
    pub fn mat_to_ntt(&self, m: &[Vec<Vec<u64>>]) -> CrtNttMat {
        m.iter()
            .map(|row| row.iter().map(|p| self.to_ntt(p)).collect())
            .collect()
    }

    /// Matrix-vector product where the matrix is already in NTT form.
    ///
    /// The vector is transformed once per input polynomial and then reused
    /// across every output row.  This avoids the two forward NTTs per
    /// matrix coefficient that `mul` would otherwise perform.
    pub fn mat_vec_ntt(&self, m: &CrtNttMat, v: &[Vec<u64>]) -> Vec<Vec<u64>> {
        assert!(!m.is_empty());
        assert_eq!(m[0].len(), v.len());

        let v_ntt: Vec<CrtNttPoly> = v.iter().map(|p| self.to_ntt(p)).collect();
        let mut out = Vec::with_capacity(m.len());
        for row in m {
            assert_eq!(row.len(), v_ntt.len());
            let mut acc1 = vec![0u64; self.d];
            let mut acc2 = vec![0u64; self.d];
            for (a, b) in row.iter().zip(v_ntt.iter()) {
                for i in 0..self.d {
                    let p1 = mul_mod_aux::<P1Marker>(a.c1[i], b.c1[i]);
                    acc1[i] = add_mod_aux::<P1Marker>(acc1[i], p1);
                    let p2 = mul_mod_aux::<P2Marker>(a.c2[i], b.c2[i]);
                    acc2[i] = add_mod_aux::<P2Marker>(acc2[i], p2);
                }
            }
            self.ctx1.inverse_in_place(&mut acc1);
            self.ctx2.inverse_in_place(&mut acc2);
            out.push(self.crt_combine(&acc1, &acc2));
        }
        out
    }

    fn split(&self, a: &[u64]) -> (Vec<u64>, Vec<u64>) {
        let (q, hq) = (self.q, self.half_q);
        debug_assert!(AUX_P1 > q && AUX_P2 > q);
        let mut a1 = Vec::with_capacity(self.d);
        let mut a2 = Vec::with_capacity(self.d);
        for &c in a {
            // Centre-lift from [0,q) to [-(q-1)/2,(q-1)/2] and reduce
            // into each auxiliary field.  Avoid the coefficient-dependent
            // branch `if c > q/2 { c-q } else { c }`; since q < p1,p2
            // for every supported parameter set, the negative case is
            // represented as `c + p - q`.
            let neg_mask = 0u64.wrapping_sub((c > hq) as u64);
            a1.push(c.wrapping_add((AUX_P1 - q) & neg_mask));
            a2.push(c.wrapping_add((AUX_P2 - q) & neg_mask));
        }
        (a1, a2)
    }

    fn crt_combine(&self, c1: &[u64], c2: &[u64]) -> Vec<u64> {
        let q_i128 = self.q as i128;
        let p1 = AUX_P1;
        let p2 = AUX_P2;
        let inv = self.p1_inv_mod_p2;
        let half_p: i128 = (AUX_P / 2) as i128;
        let big_p: i128 = AUX_P as i128;

        let mut out = Vec::with_capacity(self.d);
        for i in 0..self.d {
            let r1 = c1[i];
            let r2 = c2[i];
            // t = (r2 - r1) * inv mod p2
            let diff: u64 = sub_mod_aux::<P2Marker>(r2, r1 % p2);
            let t = mul_mod_aux::<P2Marker>(diff, inv);
            // x = r1 + t * p1  in  [0, P)
            let x: u128 = r1 as u128 + t as u128 * p1 as u128;
            let mut xi: i128 = x as i128;
            let center_mask = 0i128.wrapping_sub((xi > half_p) as i128);
            xi = xi.wrapping_sub(big_p & center_mask);
            let mut y = xi % q_i128;
            let neg_mask = 0i128.wrapping_sub((y < 0) as i128);
            y = y.wrapping_add(q_i128 & neg_mask);
            out.push(y as u64);
        }
        out
    }
}

#[inline]
fn mul_pointwise_aux<T: AuxPrime>(a: &[u64], b: &[u64]) -> Vec<u64> {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b.iter())
        .map(|(&x, &y)| mul_mod_aux::<T>(x, y))
        .collect()
}

// -------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn schoolbook(a: &[u64], b: &[u64], q: u64) -> Vec<u64> {
        let d = a.len();
        let mut c = vec![0i128; d];
        for i in 0..d {
            for j in 0..d {
                let prod = (a[i] as i128) * (b[j] as i128);
                let k = i + j;
                if k < d {
                    c[k] += prod;
                } else {
                    c[k - d] -= prod;
                }
            }
        }
        c.into_iter()
            .map(|v| {
                let m = v.rem_euclid(q as i128);
                m as u64
            })
            .collect()
    }

    fn rng(seed: u64) -> impl FnMut(u64) -> u64 {
        // Linear-congruential — fine for test fixtures, not crypto.
        let mut s = seed;
        move |bound: u64| {
            s = s
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (s >> 16) % bound
        }
    }

    #[test]
    fn aux_primes_support_d256() {
        for &p in &[AUX_P1, AUX_P2] {
            assert_eq!((p - 1) % 512, 0);
        }
    }

    #[test]
    fn primitive_root_has_exact_order_2d() {
        for &p in &[AUX_P1, AUX_P2] {
            for &d in &[128, 256] {
                let psi = find_primitive_2d_root(p, d).unwrap();
                assert_eq!(pow_mod(psi, d as u64, p), p - 1);
                assert_eq!(pow_mod(psi, 2 * d as u64, p), 1);
            }
        }
    }

    #[test]
    fn forward_inverse_roundtrip() {
        for &p in &[AUX_P1, AUX_P2] {
            let d = 128usize;
            let psi = find_primitive_2d_root(p, d).unwrap();
            let zetas = make_zeta_table(p, d, psi);
            let mut r = rng(0xABCDEF01);
            let original: Vec<u64> = (0..d).map(|_| r(p)).collect();
            let mut f = original.clone();
            ntt(&mut f, &zetas, p);
            intt(&mut f, &zetas, p);
            assert_eq!(f, original);
        }
    }

    #[test]
    fn crt_mul_matches_schoolbook() {
        use crate::params::PRODUCTION_PARAMS as PP;
        for &q in &[PP.q, PP.q_hat] {
            for &d in &[128usize] {
                let backend = CrtBackend::new(q, d).unwrap();
                let mut r = rng(0x1234_5678);
                for _ in 0..5 {
                    let a: Vec<u64> = (0..d).map(|_| r(q)).collect();
                    let b: Vec<u64> = (0..d).map(|_| r(q)).collect();
                    let want = schoolbook(&a, &b, q);
                    let got = backend.mul(&a, &b);
                    assert_eq!(got, want, "d={}, q={}", d, q);
                }
            }
        }
    }

    #[test]
    fn crt_backend_rejects_q_beyond_reconstruction_bound() {
        // q near 2^48 at d = 256 blows the CRT bound: 256 * (2^47)^2 = 2^102
        // vastly exceeds P/2 ≈ 2^95.  Construction must refuse.
        assert!(CrtBackend::new(281_474_976_710_597, 256).is_none());
        // Same q at d = 128 also fails: 128 * (2^47)^2 = 2^101 > 2^95.
        assert!(CrtBackend::new(281_474_976_710_597, 128).is_none());
    }

    #[test]
    fn crt_backend_accepts_production_moduli() {
        use crate::params::PRODUCTION_PARAMS as PP;
        assert!(CrtBackend::new(PP.q, 128).is_some());
        assert!(CrtBackend::new(PP.q_hat, 128).is_some());
        assert!(CrtBackend::new(PP.q, 256).is_some());
        assert!(CrtBackend::new(PP.q_hat, 256).is_some());
    }

    #[test]
    fn crt_mul_edge_cases() {
        use crate::params::PRODUCTION_PARAMS as PP;
        let d = 128usize;
        let backend = CrtBackend::new(PP.q, d).unwrap();
        let mut x = vec![0u64; d];
        x[1] = 1; // X
        let mut x127 = vec![0u64; d];
        x127[127] = 1; // X^127
        let prod = backend.mul(&x, &x127);
        // expected: -1 at position 0 (i.e. q-1), zeros elsewhere
        let mut expected = vec![0u64; d];
        expected[0] = PP.q - 1;
        assert_eq!(prod, expected);
    }
}
