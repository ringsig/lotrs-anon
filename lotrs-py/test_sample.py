"""
test_sample.py -- Unit tests for XOF-based samplers.
"""

import math
from Crypto.Hash import SHAKE256

from sample import (build_cdt, make_xof,
                    derive_subseed,
                    xof_sample_uniform, xof_sample_gaussian,
                    build_cdt_sampler, build_facct_sampler,
                    xof_sample_gaussian_facct,
                    _prepare_large_gaussian_sampler,
                    _facct_exp_neg_threshold,
                    xof_sample_short, xof_sample_challenge,
                    xof_sample_ternary, xof_sample_bounded,
                    needs_large_gaussian_sampler)


# ---- determinism ---------------------------------------------------------

def test_uniform_determinism():
    a = xof_sample_uniform(make_xof(b"\x00" * 32, b"u"), 7681, 64)
    b = xof_sample_uniform(make_xof(b"\x00" * 32, b"u"), 7681, 64)
    assert a == b


def test_gaussian_determinism():
    cdt = build_cdt(10.0, lam=128)
    a = xof_sample_gaussian(make_xof(b"\x01" * 32, b"g"), cdt, 128, 64)
    b = xof_sample_gaussian(make_xof(b"\x01" * 32, b"g"), cdt, 128, 64)
    assert a == b


def test_challenge_determinism():
    a = xof_sample_challenge(make_xof(b"\x02" * 32, b"c"), 10, 64)
    b = xof_sample_challenge(make_xof(b"\x02" * 32, b"c"), 10, 64)
    assert a == b


# ---- domain separation ---------------------------------------------------

def test_different_tags_differ():
    a = xof_sample_uniform(make_xof(b"\x00" * 32, b"A"), 7681, 16)
    b = xof_sample_uniform(make_xof(b"\x00" * 32, b"B"), 7681, 16)
    assert a != b


def test_different_seeds_differ():
    a = xof_sample_uniform(make_xof(b"\x00" * 32, b"X"), 7681, 16)
    b = xof_sample_uniform(make_xof(b"\x01" * 32, b"X"), 7681, 16)
    assert a != b


def test_buffered_xof_matches_raw_shake_stream():
    raw = SHAKE256.new(b"\x11" * 32)
    raw.update(b"buf")
    xof = make_xof(b"\x11" * 32, b"buf")
    chunks = [1, 7, 8, 31, 2, 64, 5, 96]
    for n in chunks:
        assert xof.read(n) == raw.read(n)


def test_derive_subseed_matches_single_block_expand():
    want = make_xof(b"\x12" * 32, b"sub", 7).read(32)
    got = derive_subseed(b"\x12" * 32, b"sub", 7)
    assert got == want


# ---- uniform sampler bounds ----------------------------------------------

def test_uniform_range():
    for q in (7681, 8380417, 101):
        poly = xof_sample_uniform(make_xof(b"\x03" * 32, b"U"), q, 64)
        assert all(0 <= c < q for c in poly), f"out of [0,{q})"


# ---- Gaussian sampler ----------------------------------------------------

def test_gaussian_tail():
    sigma = 10.0
    cdt = build_cdt(sigma, lam=128)
    poly = xof_sample_gaussian(
        make_xof(b"\x04" * 32, b"G"), cdt, 128, 1024)
    assert max(abs(c) for c in poly) <= 14 * sigma


def test_gaussian_symmetry():
    """Mean of a large sample should be near zero."""
    sigma = 10.0
    cdt = build_cdt(sigma, lam=128)
    poly = xof_sample_gaussian(
        make_xof(b"\x05" * 32, b"G"), cdt, 128, 4096)
    mean = sum(poly) / len(poly)
    assert abs(mean) < 3 * sigma / (len(poly) ** 0.5), \
        f"mean {mean:.3f} too far from zero"


def test_needs_large_gaussian_sampler_threshold():
    # Decision surface around MAX_CDT_ENTRIES = 1_000_000.
    assert needs_large_gaussian_sampler(1.0e7)
    assert not needs_large_gaussian_sampler(1.0e4)


def test_build_facct_sampler_produces_facct_kind():
    assert build_facct_sampler(1.0e7, lam=128).kind == "facct"
    # Also usable at moderate sigma if the caller deliberately picks
    # it, e.g. for tests that compare CDT and FACCT on the same sigma.
    assert build_facct_sampler(100.0, lam=128).kind == "facct"


def test_build_cdt_sampler_produces_cdt_kind():
    assert build_cdt_sampler(10.0, lam=128).kind == "cdt"


def test_facct_sampler_determinism():
    sampler = build_facct_sampler(1.0e7, lam=128)
    a = xof_sample_gaussian_facct(
        make_xof(b"\x0a" * 32, b"LG"), sampler.data, 32)
    b = xof_sample_gaussian_facct(
        make_xof(b"\x0a" * 32, b"LG"), sampler.data, 32)
    assert a == b


def test_facct_sampler_has_expected_scale():
    sigma = 1.0e6
    sampler = build_facct_sampler(sigma, lam=128)
    poly = xof_sample_gaussian_facct(
        make_xof(b"\x0b" * 32, b"LG"), sampler.data, 2048)
    mean_sq = sum(c * c for c in poly) / len(poly)
    # Very loose sanity check: RMS should be in the right ballpark.
    assert 0.5 * sigma < math.sqrt(mean_sq) < 1.5 * sigma


def test_facct_sampler_is_roughly_centered():
    sigma = 1.0e6
    sampler = build_facct_sampler(sigma, lam=128)
    poly = xof_sample_gaussian_facct(
        make_xof(b"\x0c" * 32, b"LG"), sampler.data, 4096)
    mean = sum(poly) / len(poly)
    assert abs(mean) < 3 * sigma / (len(poly) ** 0.5)


def test_facct_exp_remainder_threshold_is_accurate():
    params = _prepare_large_gaussian_sampler(1.0e6)
    u_bits = params["u_bits"]
    prob_scale = params["prob_scale"]
    for frac in range(0, 257):
        r = math.log(2.0) * frac / 256.0
        r_q = int(round(r * (1 << u_bits)))
        got = _facct_exp_neg_threshold(
            r_q, u_bits, prob_scale, params["exp_poly"]) / prob_scale
        want = math.exp(-r)
        assert abs(got - want) < 3e-15


def test_explicit_builders_reject_nonpositive_sigma():
    # Both builders share the same sigma guard.
    for build in (build_cdt_sampler, build_facct_sampler):
        for bad in (0.0, -1.0, float("nan")):
            try:
                build(bad, lam=128)
            except ValueError:
                continue
            raise AssertionError(
                f"expected ValueError for sigma={bad!r} in {build.__name__}")


def test_facct_matches_cdt_distribution_window():
    sigma = 100.0
    n = 32768
    cdt = build_cdt(sigma, lam=128)
    facct = build_facct_sampler(sigma, lam=128)
    a = xof_sample_gaussian(
        make_xof(b"\x0d" * 32, b"CDT"), cdt, 128, n)
    b = xof_sample_gaussian_facct(
        make_xof(b"\x0e" * 32, b"FACCT"), facct.data, n)
    for k in range(-3, 4):
        pa = sum(1 for x in a if x == k) / n
        pb = sum(1 for x in b if x == k) / n
        assert abs(pa - pb) < 0.01
    ea = sum(x * x for x in a) / n
    eb = sum(x * x for x in b) / n
    assert abs(math.sqrt(ea) - math.sqrt(eb)) < 5.0


# ---- short sampler -------------------------------------------------------

def test_short_bounds():
    for eta in (1, 2, 4):
        poly = xof_sample_short(
            make_xof(b"\x06" * 32, b"S"), eta, 256)
        assert all(-eta <= c <= eta for c in poly)


# ---- challenge sampler ---------------------------------------------------

def test_challenge_weight():
    for w in (4, 10, 20):
        ch = xof_sample_challenge(
            make_xof(b"\x07" * 32, b"C", w), w, 128)
        assert sum(abs(c) for c in ch) == w, "L1 != w"
        assert max(abs(c) for c in ch) == 1, "Linf != 1"
        assert sum(1 for c in ch if c != 0) == w, "wrong support size"


def test_challenge_signs():
    """Both +1 and -1 should appear over many samples."""
    positives, negatives = 0, 0
    for seed_byte in range(20):
        ch = xof_sample_challenge(
            make_xof(bytes([seed_byte]) * 32, b"CS"), 10, 64)
        positives += sum(1 for c in ch if c == 1)
        negatives += sum(1 for c in ch if c == -1)
    assert positives > 0 and negatives > 0


# ---- ternary sampler -----------------------------------------------------

def test_ternary_range():
    poly = xof_sample_ternary(make_xof(b"\x08" * 32, b"T"), 256)
    assert all(c in (-1, 0, 1) for c in poly)


# ---- bounded sampler -----------------------------------------------------

def test_bounded_range():
    for bound in (1, 5, 10):
        poly = xof_sample_bounded(
            make_xof(b"\x09" * 32, b"B", bound), bound, 128)
        assert all(-bound <= c <= bound for c in poly)


# ---- CDT construction ----------------------------------------------------

def test_cdt_monotone():
    cdt = build_cdt(10.0, lam=128)
    for i in range(1, len(cdt)):
        assert cdt[i] >= cdt[i - 1], "CDT not monotone"


def test_cdt_terminates():
    cdt = build_cdt(10.0, lam=128)
    assert cdt[-1] == (1 << 128)


# ---- run all -------------------------------------------------------------

if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"  {t.__name__}: ok")
    print(f"test_sample.py: {len(tests)} tests passed")
