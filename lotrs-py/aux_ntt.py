"""
aux_ntt.py -- Exact negacyclic multiplication in R_q = Z_q[X]/(X^d + 1)
via CRT across two auxiliary NTT-friendly primes.

The scheme's moduli q, q_hat only need to satisfy the paper-facing
conditions (primality, q ≡ 5 mod 8, Lemma 1 size constraints).  They do
NOT need to support a native NTT.  Fast multiplication is delegated to
two 48-bit auxiliary primes p1, p2 that both support a standard radix-2
negacyclic NTT at d = 128 and d = 256:

    p1 = 281474976694273  ==  2^48 - 16383
    p2 = 281474976690689  ==  2^48 - 19967

Both primes satisfy  p ≡ 1 (mod 512), so a primitive 2d-th root of unity
exists in Z_p for every d ∈ {128, 256}.  Their product

    P = p1 * p2 ≈ 2^96

is comfortably larger than 2 * max|c_k| for every combination of
(d, q) and (d, q_hat) used in the scheme, even with inner-product
accumulation up to 256 terms — so CRT reconstruction recovers the
*exact* integer product coefficient, not a residue.

Pipeline for one multiplication in R_q of two operands in [0, q):

    1. Centre-lift coefficients from [0, q)  to  (-q/2, q/2].
    2. Reduce each coefficient mod p1 and mod p2.
    3. Forward negacyclic NTT of length d in each auxiliary prime.
    4. Pointwise multiply in each.
    5. Inverse negacyclic NTT.
    6. CRT-reconstruct the exact integer coefficient (signed).
    7. Reduce mod q (or mod q_hat for the binary-proof ring).

Bit budget (headroom against P/2 = 2^95):

                           worst |c_k|     inner <= 256
    d=128, q   < 2^38      < 2^81          < 2^89
    d=128, q̂ < 2^34      < 2^73          < 2^81
    d=256, q   < 2^38      < 2^82          < 2^90
    d=256, q̂ < 2^34      < 2^74          < 2^82

Optimisation note (Mersenne-style reduction)
--------------------------------------------
Because each prime is very close to 2^48, modular reduction in C / asm
can skip the general divison: for a 96-bit intermediate  x = hi·2^48+lo
with  hi, lo ∈ [0, 2^48),

    x mod p1  ≡  (hi · 16383 + lo)  mod p1
    x mod p2  ≡  (hi · 19967 + lo)  mod p2

since 2^48 ≡ 16383 (mod p1) and 2^48 ≡ 19967 (mod p2).  The reduced
value needs one additional trim to land in [0, p), but the multiply
constant is only 14-15 bits, so the whole reduction is two int64 fused
multiply-adds on typical hardware.

In pure Python that trick has no measurable benefit — CPython's `%` on
modest big-ints is already efficient — so we use the direct `% p` form
here and leave the Mersenne form for a C / Rust / Go port.
"""


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

AUX_P1 = 281474976694273        # 2^48 - 16383
AUX_P2 = 281474976690689        # 2^48 - 19967
AUX_P  = AUX_P1 * AUX_P2
AUX_HALF_P = AUX_P // 2

# Mersenne-style reduction constants — see module docstring.  Unused in the
# Python reference path; kept as named constants for C / Rust ports.
AUX_MERSENNE_C1 = 16383
AUX_MERSENNE_C2 = 19967

SUPPORTED_D = (128, 256)


# ---------------------------------------------------------------------------
# Primitive 2d-th root of unity
# ---------------------------------------------------------------------------

def find_primitive_2d_root(p: int, d: int) -> int:
    """Smallest x in [2, p) whose d-th power is -1 mod p.

    Equivalently, x has multiplicative order exactly 2d.  Requires
    2d | (p - 1).
    """
    if (p - 1) % (2 * d) != 0:
        raise ValueError(f"2d = {2*d} does not divide p - 1")
    exponent = (p - 1) // (2 * d)
    for x in range(2, p):
        psi = pow(x, exponent, p)
        if pow(psi, d, p) == p - 1:
            return psi
    raise ValueError(f"no primitive {2*d}-th root of unity mod {p}")


# ---------------------------------------------------------------------------
# Bit-reversed twiddle table (FIPS 204 convention)
# ---------------------------------------------------------------------------

def _bitrev(k: int, bits: int) -> int:
    r = 0
    for _ in range(bits):
        r = (r << 1) | (k & 1)
        k >>= 1
    return r


def make_zeta_table(p: int, d: int, psi: int) -> list[int]:
    """zetas[k] = psi^bitrev(k, log2 d) mod p, for k in [0, d)."""
    bits = d.bit_length() - 1
    return [pow(psi, _bitrev(k, bits), p) for k in range(d)]


# ---------------------------------------------------------------------------
# Radix-2 negacyclic NTT  (Cooley-Tukey, FIPS 204 Algs. 41 / 42)
# ---------------------------------------------------------------------------

def ntt(f: list[int], zetas: list[int], p: int) -> list[int]:
    """Forward negacyclic NTT in R_p = Z_p[X]/(X^d + 1).

    Input  in natural order (length d, values in [0, p)).
    Output in bit-reversed order.
    """
    f = f.copy()
    d = len(f)
    m = 0
    le = d // 2
    while le >= 1:
        st = 0
        while st < d:
            m += 1
            z = zetas[m]
            for j in range(st, st + le):
                t = z * f[j + le] % p
                f[j + le] = (f[j] - t) % p
                f[j] = (f[j] + t) % p
            st += 2 * le
        le //= 2
    return f


def intt(f: list[int], zetas: list[int], p: int) -> list[int]:
    """Inverse negacyclic NTT.  Input bit-reversed; output natural order.

    Includes the final scaling by d^-1 mod p.
    """
    f = f.copy()
    d = len(f)
    m = d
    le = 1
    while le < d:
        st = 0
        while st < d:
            m -= 1
            z = (-zetas[m]) % p
            for j in range(st, st + le):
                t = f[j]
                f[j] = (t + f[j + le]) % p
                f[j + le] = z * (t - f[j + le]) % p
            st += 2 * le
        le *= 2
    inv_d = pow(d, -1, p)
    return [inv_d * x % p for x in f]


# ---------------------------------------------------------------------------
# Per-prime context
# ---------------------------------------------------------------------------

class _AuxContext:
    """Precomputed constants for one auxiliary prime at one dimension d."""

    __slots__ = ("p", "d", "psi", "zetas")

    def __init__(self, p: int, d: int):
        self.p = p
        self.d = d
        self.psi = find_primitive_2d_root(p, d)
        self.zetas = make_zeta_table(p, d, self.psi)

    def forward(self, f: list[int]) -> list[int]:
        return ntt(f, self.zetas, self.p)

    def inverse(self, F: list[int]) -> list[int]:
        return intt(F, self.zetas, self.p)


# ---------------------------------------------------------------------------
# CRT backend  (public API for ring.py)
# ---------------------------------------------------------------------------

class CRTBackend:
    """Exact negacyclic multiplication in R_q = Z_q[X]/(X^d + 1).

    Uses the two auxiliary primes AUX_P1, AUX_P2 and radix-2 negacyclic
    NTTs of length d.  Inputs and outputs are in [0, q).
    """

    __slots__ = ("q", "d", "half_q", "ctx1", "ctx2", "_p1_inv_mod_p2")

    def __init__(self, q: int, d: int):
        if d not in SUPPORTED_D:
            raise ValueError(f"unsupported d={d}; must be in {SUPPORTED_D}")
        self.q = q
        self.d = d
        self.half_q = q // 2
        self.ctx1 = _AuxContext(AUX_P1, d)
        self.ctx2 = _AuxContext(AUX_P2, d)
        self._p1_inv_mod_p2 = pow(AUX_P1, -1, AUX_P2)

    def mul(self, a: list[int], b: list[int]) -> list[int]:
        """Negacyclic product a * b mod X^d + 1, reduced mod q."""
        a1, a2 = self._split(a)
        b1, b2 = self._split(b)

        A1 = self.ctx1.forward(a1)
        B1 = self.ctx1.forward(b1)
        c1 = self.ctx1.inverse([A1[i] * B1[i] % AUX_P1 for i in range(self.d)])

        A2 = self.ctx2.forward(a2)
        B2 = self.ctx2.forward(b2)
        c2 = self.ctx2.inverse([A2[i] * B2[i] % AUX_P2 for i in range(self.d)])

        return self._crt_combine(c1, c2)

    # ---- helpers ---------------------------------------------------------

    def _split(self, a: list[int]) -> tuple[list[int], list[int]]:
        """Centre-lift from [0, q) and reduce mod the two auxiliary primes."""
        q, hq = self.q, self.half_q
        a1 = [0] * self.d
        a2 = [0] * self.d
        for i, c in enumerate(a):
            cc = c - q if c > hq else c          # centred integer
            a1[i] = cc % AUX_P1
            a2[i] = cc % AUX_P2
        return a1, a2

    def _crt_combine(self, c1: list[int], c2: list[int]) -> list[int]:
        """Garner CRT: reconstruct centred integer, then reduce mod q."""
        q, P = self.q, AUX_P
        P1, P2 = AUX_P1, AUX_P2
        inv = self._p1_inv_mod_p2
        half_P = AUX_HALF_P
        out = [0] * self.d
        for i in range(self.d):
            r1, r2 = c1[i], c2[i]
            t = (r2 - r1) * inv % P2              # t in [0, P2)
            x = r1 + t * P1                       # x in [0, P)
            if x > half_P:
                x -= P                            # centre-lift
            out[i] = x % q
        return out


# ---------------------------------------------------------------------------
# Module self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import random

    def _schoolbook(a, b, q, d):
        c = [0] * d
        for i in range(d):
            for j in range(d):
                k = i + j
                if k < d:
                    c[k] = (c[k] + a[i] * b[j]) % q
                else:
                    c[k - d] = (c[k - d] - a[i] * b[j]) % q
        return c

    # prime / root sanity
    for p in (AUX_P1, AUX_P2):
        assert (p - 1) % 512 == 0                  # supports d=128 and d=256
    print(f"p1 = {AUX_P1}  (2^48 - {AUX_MERSENNE_C1})")
    print(f"p2 = {AUX_P2}  (2^48 - {AUX_MERSENNE_C2})")
    print(f"P  = {AUX_P}  ({AUX_P.bit_length()} bits)")

    rng = random.Random(1)
    for d in SUPPORTED_D:
        for q_scheme in (274877906837, 8589934237):
            backend = CRTBackend(q_scheme, d)
            for _ in range(3):
                a = [rng.randrange(q_scheme) for _ in range(d)]
                b = [rng.randrange(q_scheme) for _ in range(d)]
                want = _schoolbook(a, b, q_scheme, d)
                got = backend.mul(a, b)
                assert want == got, f"mismatch at d={d}, q={q_scheme}"
    print("aux_ntt.py: all self-tests passed")
