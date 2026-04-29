"""
test_e2e.py -- End-to-end integration tests for LoTRS.

Tests the complete pipeline: keygen -> sign -> serialize -> deserialize
-> verify, plus various negative / tampering scenarios.
"""

import copy
import math

from params import TEST_PARAMS
from lotrs import LoTRS
from ring import Ring
from codec import LoTRSCodec


_scheme = None
_codec = None


def _get():
    global _scheme, _codec
    if _scheme is None:
        _scheme = LoTRS(TEST_PARAMS)
        _codec = LoTRSCodec(TEST_PARAMS)
    return _scheme, _codec


def _make_ring(par):
    scheme, _ = _get()
    pp = scheme.setup(b"\x00" * 32)
    pk_table, all_sks = [], []
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
#  Full pipeline: sign -> encode -> decode -> verify
# ===========================================================================

def test_full_pipeline_all_columns():
    """Complete pipeline for every hidden column."""
    scheme, codec = _get()
    par = TEST_PARAMS
    pp, pk_table, all_sks = _make_ring(par)

    for ell in range(par.N):
        mu = f"pipeline ell={ell}".encode()
        sig = scheme.sign(pp, all_sks[ell], ell, mu, pk_table,
                          bytes([ell + 0x60]) * 32)

        # encode -> decode -> verify
        enc = codec.sig_encode(sig)
        dec = codec.sig_decode(enc)
        assert scheme.verify(pp, mu, dec, pk_table), \
            f"pipeline failed for ell={ell}"


def test_full_pipeline_various_messages():
    """Different messages produce valid signatures."""
    scheme, codec = _get()
    par = TEST_PARAMS
    pp, pk_table, all_sks = _make_ring(par)

    ell = 1
    for i in range(5):
        mu = f"message #{i}".encode()
        sig = scheme.sign(pp, all_sks[ell], ell, mu, pk_table,
                          bytes([i + 0x70]) * 32)
        enc = codec.sig_encode(sig)
        dec = codec.sig_decode(enc)
        assert scheme.verify(pp, mu, dec, pk_table)


def test_empty_message():
    """Signing and verifying an empty message."""
    scheme, codec = _get()
    par = TEST_PARAMS
    pp, pk_table, all_sks = _make_ring(par)

    sig = scheme.sign(pp, all_sks[0], 0, b"", pk_table, b"\xA1" * 32)
    enc = codec.sig_encode(sig)
    dec = codec.sig_decode(enc)
    assert scheme.verify(pp, b"", dec, pk_table)


def test_long_message():
    """Signing and verifying a long message."""
    scheme, codec = _get()
    par = TEST_PARAMS
    pp, pk_table, all_sks = _make_ring(par)

    mu = b"A" * 10000
    sig = scheme.sign(pp, all_sks[0], 0, mu, pk_table, b"\xA2" * 32)
    assert scheme.verify(pp, mu, sig, pk_table)


# ===========================================================================
#  Negative tests: wrong message, wrong keys
# ===========================================================================

def test_wrong_message_all_columns():
    """Wrong message must fail for every hidden column."""
    scheme, _ = _get()
    par = TEST_PARAMS
    pp, pk_table, all_sks = _make_ring(par)

    for ell in range(par.N):
        sig = scheme.sign(pp, all_sks[ell], ell, b"correct", pk_table,
                          bytes([ell + 0x80]) * 32)
        assert not scheme.verify(pp, b"incorrect", sig, pk_table)


def test_wrong_pk_table():
    """Verification with a different PK table must fail."""
    scheme, _ = _get()
    par = TEST_PARAMS
    pp, pk_table, all_sks = _make_ring(par)

    ell = 0
    sig = scheme.sign(pp, all_sks[ell], ell, b"pk test", pk_table,
                      b"\xB1" * 32)

    # create a different PK table
    pk_table2 = []
    for col in range(par.N):
        cpks = []
        for row in range(par.T):
            seed = bytes([col + 0x40, row + 0x80]) + b"\x00" * 30
            _, pk = scheme.keygen(pp, seed)
            cpks.append(pk)
        pk_table2.append(cpks)

    assert not scheme.verify(pp, b"pk test", sig, pk_table2)


# ===========================================================================
#  Tampering tests
# ===========================================================================

def test_tamper_B_bin_hi():
    """Corrupting B_bin_hi should break verification."""
    scheme, _ = _get()
    par = TEST_PARAMS
    pp, pk_table, all_sks = _make_ring(par)

    sig = scheme.sign(pp, all_sks[1], 1, b"tamper B", pk_table,
                      b"\xC1" * 32)
    assert scheme.verify(pp, b"tamper B", sig, pk_table)

    bhi = sig["pi"]["B_bin_hi"]
    bhi[0] = list(bhi[0])
    bhi[0][0] = (bhi[0][0] + 1) % par.q_hat
    assert not scheme.verify(pp, b"tamper B", sig, pk_table)


def test_tamper_f1():
    """Corrupting f1 should break verification."""
    scheme, _ = _get()
    par = TEST_PARAMS
    pp, pk_table, all_sks = _make_ring(par)

    sig = scheme.sign(pp, all_sks[0], 0, b"tamper f1", pk_table,
                      b"\xC2" * 32)
    assert scheme.verify(pp, b"tamper f1", sig, pk_table)

    f1 = sig["pi"]["f1"]
    f1[0] = list(f1[0])
    f1[0][0] = (f1[0][0] + 1) % par.q
    assert not scheme.verify(pp, b"tamper f1", sig, pk_table)


def test_tamper_z_b():
    """Corrupting z_b should break verification."""
    scheme, _ = _get()
    par = TEST_PARAMS
    pp, pk_table, all_sks = _make_ring(par)

    sig = scheme.sign(pp, all_sks[2], 2, b"tamper zb", pk_table,
                      b"\xC3" * 32)
    assert scheme.verify(pp, b"tamper zb", sig, pk_table)

    zb = sig["pi"]["z_b"]
    zb[0] = list(zb[0])
    zb[0][0] = (zb[0][0] + 1) % par.q_hat
    assert not scheme.verify(pp, b"tamper zb", sig, pk_table)


def test_tamper_x():
    """Corrupting the challenge x should break verification."""
    scheme, _ = _get()
    par = TEST_PARAMS
    pp, pk_table, all_sks = _make_ring(par)

    sig = scheme.sign(pp, all_sks[1], 1, b"tamper x", pk_table,
                      b"\xC4" * 32)
    assert scheme.verify(pp, b"tamper x", sig, pk_table)

    x = sig["pi"]["x"]
    # find a nonzero position and flip its sign
    for i in range(par.d):
        c = x[i] if x[i] <= par.q // 2 else x[i] - par.q
        if c != 0:
            x[i] = (-c) % par.q
            break
    assert not scheme.verify(pp, b"tamper x", sig, pk_table)


# ===========================================================================
#  Signature determinism
# ===========================================================================

def test_deterministic_signatures():
    """Same inputs produce identical signatures."""
    scheme, codec = _get()
    par = TEST_PARAMS
    pp, pk_table, all_sks = _make_ring(par)

    ell, mu, seed = 1, b"determinism", b"\xD1" * 32
    sig1 = scheme.sign(pp, all_sks[ell], ell, mu, pk_table, seed)
    sig2 = scheme.sign(pp, all_sks[ell], ell, mu, pk_table, seed)

    enc1 = codec.sig_encode(sig1)
    enc2 = codec.sig_encode(sig2)
    assert enc1 == enc2, "signatures should be byte-identical"


# ===========================================================================
#  Signature structure
# ===========================================================================

def test_signature_contains_required_fields():
    """Aggregated signature must have the expected top-level fields."""
    scheme, _ = _get()
    par = TEST_PARAMS
    pp, pk_table, all_sks = _make_ring(par)

    sig = scheme.sign(pp, all_sks[0], 0, b"fields", pk_table,
                      b"\xE1" * 32)
    assert "pi" in sig
    assert "z_tilde" in sig
    assert "r_tilde" in sig
    assert "e_tilde" in sig
    assert len(sig["z_tilde"]) == par.l
    assert len(sig["r_tilde"]) == par.l_prime
    assert len(sig["e_tilde"]) == par.k


def test_signature_norms_bounded():
    """Aggregated response norms should be within expected bounds."""
    scheme, _ = _get()
    par = TEST_PARAMS
    Rq = scheme.Rq
    pp, pk_table, all_sks = _make_ring(par)

    sig = scheme.sign(pp, all_sks[1], 1, b"norms", pk_table,
                      b"\xE2" * 32)

    B_z = par.sigma_0 * par.tail_t * math.sqrt(par.T * par.d * par.l)
    B_r = par.sigma_0_prime * par.tail_t * math.sqrt(
        par.T * par.d * par.l_prime)

    z_norm = math.sqrt(Rq.vec_l2_norm_sq(sig["z_tilde"]))
    r_norm = math.sqrt(Rq.vec_l2_norm_sq(sig["r_tilde"]))

    assert z_norm <= B_z, f"z_tilde norm {z_norm:.0f} > {B_z:.0f}"
    assert r_norm <= B_r, f"r_tilde norm {r_norm:.0f} > {B_r:.0f}"


# ===========================================================================
#  Runner
# ===========================================================================

if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    print(f"Running {len(tests)} end-to-end tests ...")
    for t in tests:
        t()
        print(f"  {t.__name__}: ok")
    print(f"test_e2e.py: {len(tests)} tests passed")
