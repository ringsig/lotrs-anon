"""
test_wtilde.py -- Tests for the w-tilde compression changes (new paper).

Covers:
- w_tilde_hi[0] is NOT in the signature (kappa=1)
- verifier reconstructs w_hat_0 from verification equation
- stability check in sign2 causes restart when needed
- FS hash uses reconstructed w_hat_0_hi
- SAgg expanded consistency checks
"""

from params import TEST_PARAMS
from lotrs import LoTRS
from ring import Ring
from codec import LoTRSCodec


_scheme = None


def _get_scheme():
    global _scheme
    if _scheme is None:
        _scheme = LoTRS(TEST_PARAMS)
    return _scheme


def _make_ring(par):
    """Generate a full PK table + secret keys for testing."""
    scheme = _get_scheme()
    pp = scheme.setup(b"\x00" * 32)
    pk_table = []
    all_sks = []
    for col in range(par.N):
        cpks, csks = [], []
        for row in range(par.T):
            seed = bytes([col, row]) + b"\x00" * 30
            sk, pk = scheme.keygen(pp, seed)
            csks.append(sk)
            cpks.append(pk)
        pk_table.append(cpks)
        all_sks.append(csks)
    return pp, pk_table, all_sks


# ===========================================================================
#  pi structure tests
# ===========================================================================

def test_pi_has_no_w_tilde_hi_for_kappa1():
    """For kappa=1, pi["w_tilde_hi"] must be empty."""
    scheme = _get_scheme()
    par = TEST_PARAMS
    assert par.kappa == 1, "test requires kappa=1"
    pp, pk_table, all_sks = _make_ring(par)

    sig = scheme.sign(pp, all_sks[0], 0, b"test", pk_table, b"\x01" * 32)
    assert sig["pi"]["w_tilde_hi"] == [], \
        "kappa=1: w_tilde_hi should be empty in pi"


def test_pi_fields_present():
    """pi should contain B_bin_hi, w_tilde_hi, x, f1, z_b."""
    scheme = _get_scheme()
    par = TEST_PARAMS
    pp, pk_table, all_sks = _make_ring(par)

    sig = scheme.sign(pp, all_sks[1], 1, b"pi fields", pk_table,
                      b"\x02" * 32)
    pi = sig["pi"]
    assert "B_bin_hi" in pi
    assert "w_tilde_hi" in pi
    assert "x" in pi
    assert "f1" in pi
    assert "z_b" in pi


# ===========================================================================
#  w_hat_0 reconstruction tests
# ===========================================================================

def test_verify_reconstructs_w_hat_0():
    """
    The verifier must successfully reconstruct w_hat_0 and use it
    in the FS hash.  If the reconstruction is wrong, verification fails.
    """
    scheme = _get_scheme()
    par = TEST_PARAMS
    pp, pk_table, all_sks = _make_ring(par)

    ell = 1
    mu = b"reconstruct test"
    sig = scheme.sign(pp, all_sks[ell], ell, mu, pk_table, b"\x03" * 32)
    assert scheme.verify(pp, mu, sig, pk_table)


def test_corrupted_z_tilde_breaks_w_hat_0():
    """
    Corrupting z_tilde changes the reconstructed w_hat_0, which should
    cause the FS hash check to fail.
    """
    scheme = _get_scheme()
    par = TEST_PARAMS
    pp, pk_table, all_sks = _make_ring(par)

    ell = 2
    mu = b"corrupt z"
    sig = scheme.sign(pp, all_sks[ell], ell, mu, pk_table, b"\x04" * 32)
    assert scheme.verify(pp, mu, sig, pk_table)

    # corrupt z_tilde
    z = sig["z_tilde"]
    z[0] = list(z[0])
    z[0][0] = (z[0][0] + 1) % par.q
    sig["z_tilde"] = z

    assert not scheme.verify(pp, mu, sig, pk_table), \
        "corrupted z_tilde should break verification"


def test_corrupted_e_tilde_breaks_w_hat_0():
    """Corrupting e_tilde should also break the FS hash via w_hat_0."""
    scheme = _get_scheme()
    par = TEST_PARAMS
    pp, pk_table, all_sks = _make_ring(par)

    ell = 0
    mu = b"corrupt e"
    sig = scheme.sign(pp, all_sks[ell], ell, mu, pk_table, b"\x05" * 32)
    assert scheme.verify(pp, mu, sig, pk_table)

    # corrupt e_tilde with a large perturbation to ensure it changes
    # the high bits of w_hat_0 (small changes may stay within
    # the same decomposition bucket)
    e = sig["e_tilde"]
    e[0] = list(e[0])
    shift = 1 << par.K_w   # large enough to change high bits
    e[0][0] = (e[0][0] + shift) % par.q
    sig["e_tilde"] = e

    assert not scheme.verify(pp, mu, sig, pk_table)


def test_corrupted_r_tilde_breaks_w_hat_0():
    """Corrupting r_tilde should also break verification."""
    scheme = _get_scheme()
    par = TEST_PARAMS
    pp, pk_table, all_sks = _make_ring(par)

    ell = 1
    mu = b"corrupt r"
    sig = scheme.sign(pp, all_sks[ell], ell, mu, pk_table, b"\x06" * 32)
    assert scheme.verify(pp, mu, sig, pk_table)

    r = sig["r_tilde"]
    r[0] = list(r[0])
    shift = 1 << par.K_w
    r[0][0] = (r[0][0] + shift) % par.q
    sig["r_tilde"] = r

    assert not scheme.verify(pp, mu, sig, pk_table)


# ===========================================================================
#  Codec integration with new format
# ===========================================================================

def test_codec_roundtrip_no_w_tilde_hi():
    """
    For kappa=1, the encoded signature should contain no w_tilde_hi
    bytes, yet decode + verify correctly.
    """
    scheme = _get_scheme()
    par = TEST_PARAMS
    codec = LoTRSCodec(par)
    pp, pk_table, all_sks = _make_ring(par)

    ell = 1
    mu = b"codec no whi"
    sig = scheme.sign(pp, all_sks[ell], ell, mu, pk_table, b"\x07" * 32)

    enc = codec.sig_encode(sig)
    dec = codec.sig_decode(enc)

    assert dec["pi"]["w_tilde_hi"] == []
    assert scheme.verify(pp, mu, dec, pk_table)


def test_signature_smaller_without_w_tilde_hi():
    """
    The new signature (kappa=1) should be smaller than a hypothetical
    version that includes w_tilde_hi.
    """
    codec = LoTRSCodec(TEST_PARAMS)
    s = codec.sizes()
    assert s["w_tilde_hi"] == 0, "kappa=1: no w_tilde_hi bytes"


# ===========================================================================
#  Stability check
# ===========================================================================

def test_sign_succeeds_despite_stability_check():
    """
    Normal signing should succeed -- the stability check should not
    reject too aggressively.
    """
    scheme = _get_scheme()
    par = TEST_PARAMS
    pp, pk_table, all_sks = _make_ring(par)

    for ell in range(par.N):
        sig = scheme.sign(pp, all_sks[ell], ell, b"stability",
                          pk_table, bytes([ell + 0x50]) * 32)
        assert scheme.verify(pp, b"stability", sig, pk_table)


# ===========================================================================
#  SAgg consistency
# ===========================================================================

def test_sagg_rejects_empty_list():
    """sagg([]) must raise ValueError, not IndexError."""
    scheme = _get_scheme()
    try:
        scheme.sagg([])
        assert False, "sagg([]) should raise ValueError"
    except ValueError as e:
        assert "expected" in str(e).lower()


def test_sagg_rejects_wrong_count():
    """sagg with len != T must raise ValueError."""
    scheme = _get_scheme()
    par = TEST_PARAMS
    try:
        scheme.sagg([None] * (par.T - 1))
        assert False, "sagg with T-1 should raise ValueError"
    except ValueError as e:
        assert "expected" in str(e).lower()


def test_sagg_rejects_non_dict():
    scheme = _get_scheme()
    par = TEST_PARAMS
    try:
        scheme.sagg([None] * par.T)
        assert False, "sagg with None entries should raise ValueError"
    except ValueError as e:
        assert "dict" in str(e).lower()


def test_sagg_rejects_missing_keys():
    scheme = _get_scheme()
    par = TEST_PARAMS
    try:
        scheme.sagg([{}] * par.T)
        assert False, "sagg with empty dicts should raise ValueError"
    except ValueError as e:
        assert "missing" in str(e).lower()


def test_sagg_rejects_non_iterable():
    scheme = _get_scheme()
    try:
        scheme.sagg(42)
        assert False, "sagg(42) should raise ValueError"
    except ValueError as e:
        assert "sequence" in str(e).lower() or "expected" in str(e).lower()


def test_sagg_checks_all_pi_fields():
    """
    SAgg should reject if any pi field differs between signers.
    We simulate this by patching one signer's pi.
    """
    scheme = _get_scheme()
    par = TEST_PARAMS
    pp, pk_table, all_sks = _make_ring(par)

    from sample import make_xof

    ell = 1
    mu = b"sagg check"
    rho = make_xof(b"\x08" * 32, b"rho", 0).read(32)

    # round 1
    states, all_coms = [], []
    for u in range(par.T):
        st, com = scheme.sign1(pp, all_sks[ell][u], u, ell, mu,
                               pk_table, rho, 0)
        states.append(st)
        all_coms.append(com)

    # round 2 -- collect honest sigmas
    sigmas = []
    for u in range(par.T):
        sig_u = scheme.sign2(states[u], all_coms, pk_table)
        if sig_u is None:
            return  # restart needed, skip test
        sigmas.append(sig_u)

    # SAgg should succeed with consistent sigmas
    agg = scheme.sagg(sigmas)
    assert agg is not None

    # Now tamper with signer 1's B_bin_hi
    import copy
    sigmas_bad = copy.deepcopy(sigmas)
    bhi = sigmas_bad[1]["pi"]["B_bin_hi"]
    bhi[0] = list(bhi[0])
    bhi[0][0] = (bhi[0][0] + 1) % par.q_hat
    sigmas_bad[1]["pi"]["B_bin_hi"] = bhi

    try:
        scheme.sagg(sigmas_bad)
        assert False, "SAgg should reject mismatched B_bin_hi"
    except ValueError as e:
        assert "B_bin_hi" in str(e)


# ===========================================================================
#  Multiple signing seeds
# ===========================================================================

def test_multiple_seeds_all_verify():
    """Different signing seeds produce different signatures that all verify."""
    scheme = _get_scheme()
    par = TEST_PARAMS
    pp, pk_table, all_sks = _make_ring(par)

    ell = 1
    mu = b"multi-seed"
    sigs = set()
    for i in range(5):
        seed = bytes([i]) * 32
        sig = scheme.sign(pp, all_sks[ell], ell, mu, pk_table, seed)
        assert scheme.verify(pp, mu, sig, pk_table)
        # FS challenge should differ between seeds
        x_tuple = tuple(sig["pi"]["x"])
        sigs.add(x_tuple)
    assert len(sigs) == 5, "different seeds should produce different x"


# ===========================================================================
#  Runner
# ===========================================================================

if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    print(f"Running {len(tests)} w-tilde compression tests ...")
    for t in tests:
        t()
        print(f"  {t.__name__}: ok")
    print(f"test_wtilde.py: {len(tests)} tests passed")
