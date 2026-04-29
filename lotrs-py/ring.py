"""
ring.py -- Polynomial ring  R_q = Z_q[X] / (X^d + 1).

Coefficients are Python ints in [0, q) (unsigned canonical form).
Signed / centred form [-q//2, q//2] is used where noted.

Two multiplication paths:
  * Schoolbook negacyclic convolution  (always available; correctness
    reference for tests).
  * CRT-NTT backend via two 48-bit auxiliary primes  (see aux_ntt.py).
    Enabled when d in {128, 256}.  The scheme modulus q need not be
    NTT-friendly; fast multiplication is carried by the auxiliary primes
    and reconstruction is exact because P = p1 * p2 ≈ 2^96 exceeds every
    possible integer product coefficient for the sizes of q, q_hat, d
    used in the scheme.

Vectors of ring elements are plain Python lists; matrices are lists of
rows (each row itself a list of ring elements).
"""

from aux_ntt import CRTBackend, SUPPORTED_D


class Ring:
    """Arithmetic in  R_q = Z_q[X] / (X^d + 1)."""

    def __init__(self, q: int, d: int):
        self.q = q
        self.d = d
        self.half_q = q // 2
        self._ntt_ctx = None
        if d in SUPPORTED_D:
            self._ntt_ctx = CRTBackend(q, d)

    # ---- element creation ------------------------------------------------

    def zero(self):
        """Additive identity."""
        return [0] * self.d

    def one(self):
        """Multiplicative identity (constant 1)."""
        r = [0] * self.d
        r[0] = 1
        return r

    def const(self, c):
        """Constant polynomial  c mod q."""
        r = [0] * self.d
        r[0] = int(c) % self.q
        return r

    # ---- coefficient form conversions ------------------------------------

    def reduce(self, a):
        """Reduce every coefficient into [0, q)."""
        q = self.q
        return [int(c) % q for c in a]

    def centered(self, a):
        """Unsigned [0,q) --> signed [-q//2, q//2]."""
        h, q = self.half_q, self.q
        return [c - q if c > h else c for c in a]

    def from_centered(self, a):
        """Signed --> unsigned [0, q)."""
        q = self.q
        return [c % q for c in a]

    # ---- element-wise arithmetic -----------------------------------------

    def add(self, a, b):
        q = self.q
        return [(a[i] + b[i]) % q for i in range(self.d)]

    def sub(self, a, b):
        q = self.q
        return [(a[i] - b[i]) % q for i in range(self.d)]

    def neg(self, a):
        q = self.q
        return [(q - c) % q for c in a]

    def scale(self, c, a):
        """Integer scalar  c  times polynomial  a."""
        q = self.q
        c = int(c) % q
        return [(c * ai) % q for ai in a]

    def mul(self, a, b):
        """Negacyclic convolution  a * b  in R_q."""
        if self._ntt_ctx is not None:
            return self._ntt_ctx.mul(a, b)
        return self._mul_schoolbook(a, b)

    def _mul_schoolbook(self, a, b):
        """Schoolbook negacyclic convolution (reference implementation)."""
        d, q = self.d, self.q
        c = [0] * d
        for i in range(d):
            ai = a[i]
            if ai == 0:
                continue
            for j in range(d):
                k = i + j
                if k < d:
                    c[k] = (c[k] + ai * b[j]) % q
                else:
                    c[k - d] = (c[k - d] - ai * b[j]) % q
        return c

    # ---- norms (operate on unsigned-canonical coefficients) ---------------

    def inf_norm(self, a):
        """||a||_inf of the centred representation."""
        h, q = self.half_q, self.q
        return max((q - c if c > h else c) for c in a)

    def l2_norm_sq(self, a):
        """||a||_2^2 of the centred representation."""
        h, q = self.half_q, self.q
        return sum((c - q) ** 2 if c > h else c * c for c in a)

    def l1_norm(self, a):
        """||a||_1 of the centred representation."""
        h, q = self.half_q, self.q
        return sum((q - c if c > h else c) for c in a)

    # ---- vector operations -----------------------------------------------
    # A "vector" is  list[poly]  where  poly = list[int].

    def vec_zero(self, n):
        return [self.zero() for _ in range(n)]

    def vec_add(self, u, v):
        return [self.add(u[i], v[i]) for i in range(len(u))]

    def vec_sub(self, u, v):
        return [self.sub(u[i], v[i]) for i in range(len(u))]

    def vec_neg(self, v):
        return [self.neg(vi) for vi in v]

    def vec_scale(self, c_poly, v):
        """Ring-element scalar  c_poly  times vector  v."""
        return [self.mul(c_poly, vi) for vi in v]

    def vec_scale_int(self, c, v):
        """Integer scalar  c  times vector  v."""
        return [self.scale(c, vi) for vi in v]

    def inner(self, u, v):
        """<u, v> = sum_i u_i * v_i  in R_q."""
        acc = self.zero()
        for i in range(len(u)):
            acc = self.add(acc, self.mul(u[i], v[i]))
        return acc

    def vec_inf_norm(self, v):
        if not v:
            return 0
        return max(self.inf_norm(vi) for vi in v)

    def vec_l2_norm_sq(self, v):
        return sum(self.l2_norm_sq(vi) for vi in v)

    # ---- matrix operations -----------------------------------------------
    # A "matrix" is  list[list[poly]]  (row-major).

    def mat_vec(self, M, v):
        """Matrix-vector product  M * v."""
        return [self.inner(row, v) for row in M]

    # ---- helpers ---------------------------------------------------------

    def vec_concat(self, *vecs):
        out = []
        for v in vecs:
            out.extend(v)
        return out

    def centered_decompose(self, a, K):
        """
        Coefficient-wise centred decomposition.

        For each coefficient c (in centred form), write
            c = high * 2^K  +  low     with  |low| <= 2^{K-1}.

        Returns (high_poly, low_poly) both in unsigned canonical form.
        """
        q, d = self.q, self.d
        h = self.half_q
        power = 1 << K
        half_power = power >> 1
        hi = [0] * d
        lo = [0] * d
        for i in range(d):
            c = a[i] - q if a[i] > h else a[i]      # signed
            # Python's divmod rounds toward -inf, so use explicit formula
            low = c % power
            if low > half_power:
                low -= power
            high = (c - low) // power
            hi[i] = high % q
            lo[i] = low % q
        return hi, lo


# --------------------------------------------------------------------------
if __name__ == "__main__":
    R = Ring(7681, 32)

    # basic identities
    a = R.one()
    z = R.zero()
    assert R.mul(a, a) == a, "1*1 == 1"
    assert R.add(z, a) == a, "0+1 == 1"

    # negacyclic property:  X^d == -1
    x_d = [0] * 32
    x_d[0] = 7680            # -1 mod 7681
    x_1 = [0] * 32
    x_1[1] = 1               # X
    prod = x_1
    for _ in range(31):       # X^32
        prod = R.mul(prod, x_1)
    assert prod == x_d, "X^d == -1 in R_q"

    # commutativity / associativity spot-check
    import random
    rng = random.Random(42)
    def rand_poly():
        return [rng.randrange(7681) for _ in range(32)]

    a, b, c = rand_poly(), rand_poly(), rand_poly()
    assert R.mul(a, b) == R.mul(b, a), "commutativity"
    assert R.mul(R.mul(a, b), c) == R.mul(a, R.mul(b, c)), "associativity"
    assert R.mul(a, R.add(b, c)) == R.add(R.mul(a, b), R.mul(a, c)), \
        "distributivity"

    # centred decomposition round-trip
    p = rand_poly()
    hi, lo = R.centered_decompose(p, 5)
    recon = R.add(R.scale(32, hi), lo)
    assert recon == p, "centred decomposition round-trip"

    print("ring.py: all self-tests passed")
