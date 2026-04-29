"""
test_ring.py -- Unit tests for polynomial ring arithmetic.
"""

import random
from ring import Ring
from params import BENCH_PARAMS, PRODUCTION_PARAMS
from aux_ntt import (AUX_P1, AUX_P2, AUX_P, SUPPORTED_D,
                     CRTBackend, find_primitive_2d_root,
                     make_zeta_table, ntt, intt)


def _rand_poly(rng, q, d):
    return [rng.randrange(q) for _ in range(d)]


# ---- basic identities ----------------------------------------------------

def test_additive_identity():
    R = Ring(7681, 32)
    a = _rand_poly(random.Random(1), 7681, 32)
    assert R.add(a, R.zero()) == a


def test_multiplicative_identity():
    R = Ring(7681, 32)
    a = _rand_poly(random.Random(2), 7681, 32)
    assert R.mul(a, R.one()) == a


def test_additive_inverse():
    R = Ring(7681, 32)
    a = _rand_poly(random.Random(3), 7681, 32)
    assert R.add(a, R.neg(a)) == R.zero()


# ---- ring axioms ---------------------------------------------------------

def test_commutativity():
    R = Ring(7681, 64)
    rng = random.Random(10)
    a, b = _rand_poly(rng, 7681, 64), _rand_poly(rng, 7681, 64)
    assert R.mul(a, b) == R.mul(b, a)


def test_associativity():
    R = Ring(7681, 32)
    rng = random.Random(20)
    a = _rand_poly(rng, 7681, 32)
    b = _rand_poly(rng, 7681, 32)
    c = _rand_poly(rng, 7681, 32)
    assert R.mul(R.mul(a, b), c) == R.mul(a, R.mul(b, c))


def test_distributivity():
    R = Ring(7681, 32)
    rng = random.Random(30)
    a = _rand_poly(rng, 7681, 32)
    b = _rand_poly(rng, 7681, 32)
    c = _rand_poly(rng, 7681, 32)
    lhs = R.mul(a, R.add(b, c))
    rhs = R.add(R.mul(a, b), R.mul(a, c))
    assert lhs == rhs


# ---- negacyclic property  X^d = -1 --------------------------------------

def test_negacyclic():
    for d in (16, 32, 64):
        R = Ring(7681, d)
        x1 = [0] * d
        x1[1] = 1
        prod = R.one()
        for _ in range(d):
            prod = R.mul(prod, x1)
        expected = [0] * d
        expected[0] = R.q - 1             # -1 mod q
        assert prod == expected, f"X^{d} != -1 for d={d}"


# ---- norms ---------------------------------------------------------------

def test_inf_norm():
    R = Ring(101, 8)
    # 50 = q//2, centred: 50 -> 50,  51 -> -50,  100 -> -1
    a = [0, 50, 51, 100, 1, 0, 0, 0]
    assert R.inf_norm(a) == 50


def test_l1_norm():
    R = Ring(101, 4)
    a = [1, 100, 0, 50]        # centred: 1, -1, 0, 50
    assert R.l1_norm(a) == 52


# ---- centred decomposition round-trip ------------------------------------

def test_centered_decompose_roundtrip():
    R = Ring(8380417, 32)
    rng = random.Random(40)
    for K in (4, 8, 12):
        a = _rand_poly(rng, 8380417, 32)
        hi, lo = R.centered_decompose(a, K)
        recon = R.add(R.scale(1 << K, hi), lo)
        assert recon == a, f"round-trip failed for K={K}"


def test_centered_decompose_low_bound():
    R = Ring(8380417, 32)
    rng = random.Random(50)
    a = _rand_poly(rng, 8380417, 32)
    for K in (4, 8):
        _, lo = R.centered_decompose(a, K)
        assert R.inf_norm(lo) <= (1 << (K - 1)), \
            f"|low| exceeded 2^(K-1) for K={K}"


# ---- vector / matrix operations ------------------------------------------

def test_mat_vec():
    R = Ring(7681, 16)
    rng = random.Random(60)
    k, l = 3, 4
    M = [[_rand_poly(rng, 7681, 16) for _ in range(l)]
         for _ in range(k)]
    v = [_rand_poly(rng, 7681, 16) for _ in range(l)]

    result = R.mat_vec(M, v)
    assert len(result) == k

    # compare with manual inner product for first row
    manual = R.zero()
    for j in range(l):
        manual = R.add(manual, R.mul(M[0][j], v[j]))
    assert result[0] == manual


# ---- scale ---------------------------------------------------------------

def test_scale():
    R = Ring(7681, 16)
    a = R.one()
    assert R.scale(5, a) == R.const(5)
    assert R.scale(-1, a) == R.neg(a)


# ---- CRT-NTT path --------------------------------------------------------

_NTT_MODULI = [
    ("q",     PRODUCTION_PARAMS.q),
    ("q_hat", PRODUCTION_PARAMS.q_hat),
]


def test_ntt_enabled_when_d_supported():
    for d in SUPPORTED_D:
        assert Ring(PRODUCTION_PARAMS.q, d)._ntt_ctx is not None
    # unsupported dimensions fall back to schoolbook
    assert Ring(PRODUCTION_PARAMS.q, 32)._ntt_ctx is None
    assert Ring(PRODUCTION_PARAMS.q, 64)._ntt_ctx is None


def test_aux_primes_ntt_friendly():
    """Both auxiliary primes support radix-2 negacyclic NTT at every supported d."""
    for p in (AUX_P1, AUX_P2):
        for d in SUPPORTED_D:
            # primitive 2d-th root: psi^d ≡ -1 (mod p)
            psi = find_primitive_2d_root(p, d)
            assert pow(psi, d, p) == p - 1
            assert pow(psi, 2 * d, p) == 1


def test_aux_crt_bound_comfortable():
    """P = p1 * p2 must dominate twice the worst-case |c_k| for every
    (d, q) combination actually used in the scheme."""
    for q in (PRODUCTION_PARAMS.q, PRODUCTION_PARAMS.q_hat):
        for d in SUPPORTED_D:
            max_ck = d * ((q - 1) // 2) ** 2
            assert 2 * max_ck < AUX_P


def test_ntt_roundtrip():
    rng = random.Random(111)
    for d in SUPPORTED_D:
        for p in (AUX_P1, AUX_P2):
            psi = find_primitive_2d_root(p, d)
            zetas = make_zeta_table(p, d, psi)
            for _ in range(3):
                f = [rng.randrange(p) for _ in range(d)]
                assert intt(ntt(f, zetas, p), zetas, p) == f


def test_ntt_mul_matches_schoolbook():
    rng = random.Random(222)
    for _name, q in _NTT_MODULI:
        R = Ring(q, 128)
        for _ in range(20):
            a = _rand_poly(rng, q, 128)
            b = _rand_poly(rng, q, 128)
            assert R.mul(a, b) == R._mul_schoolbook(a, b)


def test_ntt_mul_edge_cases():
    for _name, q in _NTT_MODULI:
        R = Ring(q, 128)
        zero = R.zero()
        one = R.one()
        assert R.mul(zero, one) == zero
        assert R.mul(one, one) == one

        # X * X^127 = X^128 = -1
        x      = [0] * 128; x[1]       = 1
        x127   = [0] * 128; x127[127]  = 1
        neg_one = [0] * 128; neg_one[0] = q - 1
        assert R.mul(x, x127) == neg_one

        # monomial pair straddling the wrap: X^65 * X^64 = X^129 = -X
        x65 = [0] * 128; x65[65] = 1
        x64 = [0] * 128; x64[64] = 1
        neg_x = [0] * 128; neg_x[1] = q - 1
        assert R.mul(x65, x64) == neg_x

        # extreme magnitudes — stresses CRT reconstruction near the bound
        a = [q - 1] * 128
        b = [q - 1] * 128
        assert R.mul(a, b) == R._mul_schoolbook(a, b)


def test_ntt_inner_matches_schoolbook():
    """mat_vec routes through inner -> mul; guard that chain at d=128."""
    rng = random.Random(333)
    q = PRODUCTION_PARAMS.q
    R_ntt = Ring(q, 128)
    R_ref = Ring(q, 128)
    R_ref._ntt_ctx = None                    # force schoolbook on reference

    k, l = 3, 4
    M = [[_rand_poly(rng, q, 128) for _ in range(l)] for _ in range(k)]
    v = [_rand_poly(rng, q, 128) for _ in range(l)]

    assert R_ntt.mat_vec(M, v) == R_ref.mat_vec(M, v)


def test_ntt_d256_matches_schoolbook():
    """Direct check at the (currently unused) d=256 path — ensures the
    backend stays future-proof for double-size parameters."""
    rng = random.Random(444)
    q = PRODUCTION_PARAMS.q
    backend = CRTBackend(q, 256)

    def ref(a, b):
        d = 256
        c = [0] * d
        for i in range(d):
            for j in range(d):
                k = i + j
                if k < d:
                    c[k] = (c[k] + a[i] * b[j]) % q
                else:
                    c[k - d] = (c[k - d] - a[i] * b[j]) % q
        return c

    for _ in range(2):
        a = [rng.randrange(q) for _ in range(256)]
        b = [rng.randrange(q) for _ in range(256)]
        assert backend.mul(a, b) == ref(a, b)


def test_find_primitive_2d_root_rejects_wrong_modulus():
    # p1 - 1 has 2-adic valuation 14, so 2d > 2^14 fails.  d = 16384 -> 2d = 2^15.
    try:
        find_primitive_2d_root(AUX_P1, 16384)
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError when 2d does not divide p-1")


# ---- run all -------------------------------------------------------------

if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"  {t.__name__}: ok")
    print(f"test_ring.py: {len(tests)} tests passed")
