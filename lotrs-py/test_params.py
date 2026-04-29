"""
test_params.py -- Unit tests for LoTRS parameter sets.
"""

import math
from params import LoTRSParams, TEST_PARAMS


# ---- parameter consistency ---------------------------------------------------

def test_N_equals_beta_kappa():
    p = TEST_PARAMS
    assert p.N == p.beta ** p.kappa

def test_d_is_power_of_two():
    p = TEST_PARAMS
    assert p.d > 0 and (p.d & (p.d - 1)) == 0

def test_check_passes():
    """The built-in consistency check should not raise."""
    TEST_PARAMS.check()

def test_gaussian_widths_positive():
    p = TEST_PARAMS
    assert p.sigma_0 > 0
    assert p.sigma_0_prime > 0
    assert p.sigma_a > 0
    assert p.sigma_b > 0

def test_gaussian_widths_within_q():
    p = TEST_PARAMS
    for sigma in (p.sigma_0, p.sigma_0_prime, p.sigma_a, p.sigma_b):
        assert 12 * sigma < p.q, \
            f"sigma={sigma:.0f} too large for q={p.q}"

def test_mu_phi_greater_than_one():
    p = TEST_PARAMS
    assert p.mu_phi > 1.0

def test_mu_is_monotone_decreasing():
    """mu(phi) should decrease as phi increases."""
    p = TEST_PARAMS
    assert p.mu(10.0) > p.mu(20.0) > p.mu(50.0)

def test_B_a_B_b_positive():
    p = TEST_PARAMS
    assert p.B_a > 0
    assert p.B_b > 0

def test_B_g0_B_g1_positive():
    p = TEST_PARAMS
    assert p.B_g0 > 0
    assert p.B_g1 > 0

def test_B_eta_w_nonnegative():
    p = TEST_PARAMS
    assert p.B_eta_w >= 0

def test_G_cols_matches_formula():
    p = TEST_PARAMS
    expected = p.n_hat + p.k_hat + 2 * p.kappa * p.beta
    assert p.G_cols == expected

def test_B_0_positive():
    p = TEST_PARAMS
    assert p.B_0 > 0
    assert p.B_0_prime > 0

def test_K_values_positive():
    p = TEST_PARAMS
    assert p.K_A > 0
    assert p.K_B > 0
    assert p.K_w > 0
    assert p.K_A > p.K_B, "K_A should exceed K_B"


# ---- custom parameter set ----------------------------------------------------

def test_custom_params():
    """Create a minimal custom parameter set and verify it."""
    p = LoTRSParams(
        name="custom",
        d=16, q=4194389, q_hat=7000061,   # primes, 5 mod 8
        kappa=1, beta=2, T=2,
        k=2, l=2, l_prime=2,
        n_hat=2, k_hat=2,
        w=2, eta=1,
        phi=12.0, phi_a=12.0, phi_b=12.0,
        K_A=10, K_B=3, K_w=4,
    )
    p.check()
    assert p.N == 2
    assert p.sigma_0 > 0


def test_eta_prime_defaults_to_eta():
    p = TEST_PARAMS
    assert p.eta_p == p.eta


def test_prime_check_rejects_wrong_residue():
    """q must be 5 mod 8 for Lemma 1."""
    try:
        p = LoTRSParams(
            name="bad_q",
            d=16, q=8380417, q_hat=7000061,  # 8380417 is 1 mod 8
            kappa=1, beta=2, T=2,
            k=2, l=2, l_prime=2, n_hat=2, k_hat=2,
            w=2, eta=1,
            phi=12.0, phi_a=12.0, phi_b=12.0,
            K_A=10, K_B=3, K_w=4,
        )
        p.check()
        assert False, "should reject q with wrong residue"
    except AssertionError as e:
        assert "5 mod 8" in str(e)


def test_security_checks_pass_production():
    from params import PRODUCTION_PARAMS, BENCH_PARAMS
    PRODUCTION_PARAMS.check_security()
    BENCH_PARAMS.check_security()


def test_mask_sampler_declared_per_parameter_set():
    # TEST uses CDT for the mask widths (sigma_0 ≈ 2172 is small);
    # BENCH_4OF32 / BENCH / PRODUCTION declare FACCT because sigma_0 is
    # in the 10^6-10^7 regime and a CDT would need 10^7-10^8 entries.
    from params import PRODUCTION_PARAMS, BENCH_PARAMS, BENCH_4OF32
    assert TEST_PARAMS.mask_sampler == "cdt"
    assert BENCH_4OF32.mask_sampler == "facct"
    assert BENCH_PARAMS.mask_sampler == "facct"
    assert PRODUCTION_PARAMS.mask_sampler == "facct"


def test_check_rejects_inconsistent_mask_sampler_declaration():
    # Declaring `cdt` when sigma_0 is clearly in the FACCT regime must
    # fail `check()` — we don't want the scheme to silently pick one
    # backend while the parameter set claims another.
    from params import PRODUCTION_PARAMS
    from dataclasses import replace
    bad = replace(PRODUCTION_PARAMS, mask_sampler="cdt")
    try:
        bad.check()
    except AssertionError as e:
        assert "cdt" in str(e) and "facct" in str(e)
    else:
        raise AssertionError("check() must reject cdt declaration at sigma_0 ~ 2.6e7")


def test_check_rejects_unknown_mask_sampler():
    from dataclasses import replace
    bad = replace(TEST_PARAMS, mask_sampler="bogus")
    try:
        bad.check()
    except AssertionError as e:
        assert "bogus" in str(e)
    else:
        raise AssertionError("check() must reject unknown mask_sampler")


def test_B_f1_tail_formula():
    """B_f1 = 1 + sqrt(2*ln(2*M/eps)) * sigma_a  with  M = kappa(beta-1)d.

    Matches the paper's Appendix (Choice of B_{f_1}, B_{f_0}, ...).
    """
    p = TEST_PARAMS
    import math
    M = p.kappa * (p.beta - 1) * p.d
    tau = math.sqrt(2.0 * math.log(2.0 * M / (p.eps_tot / 4.0)))
    expected = 1.0 + tau * p.sigma_a
    assert abs(p.B_f1 - expected) < 1e-9


def test_B_f0_tail_formula():
    """B_f0 = 1 + sqrt(2*ln(2*M/eps)) * sqrt(beta-1) * sigma_a,
    M = kappa*d."""
    p = TEST_PARAMS
    import math
    M = p.kappa * p.d
    tau = math.sqrt(2.0 * math.log(2.0 * M / (p.eps_tot / 4.0)))
    expected = 1.0 + tau * math.sqrt(p.beta - 1) * p.sigma_a
    assert abs(p.B_f0 - expected) < 1e-9


def test_B_g0_matches_paper_formula():
    """B_g0 = w*B_f0 + d*B_f0^2  per paper Appendix (bound for g_0)."""
    p = TEST_PARAMS
    expected = p.w * p.B_f0 + p.d * p.B_f0 * p.B_f0
    assert abs(p.B_g0 - expected) < 1e-9


def test_B_g1_matches_paper_formula():
    """B_g1 = w*B_f1 + d*B_f1^2  per paper Appendix (bound for g_1)."""
    p = TEST_PARAMS
    expected = p.w * p.B_f1 + p.d * p.B_f1 * p.B_f1
    assert abs(p.B_g1 - expected) < 1e-9


# ---- run all -----------------------------------------------------------------

if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"  {t.__name__}: ok")
    print(f"test_params.py: {len(tests)} tests passed")
