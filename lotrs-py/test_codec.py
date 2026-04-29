"""
test_codec.py -- Tests for serialization, deserialization, and test vectors.
"""

import json
import math
import random

from params import TEST_PARAMS
from lotrs import LoTRS
from codec import (LoTRSCodec, BitWriter, BitReader,
                   _pack_fixed, _unpack_fixed,
                   _pack_rice, _unpack_rice,
                   _pack_challenge, _unpack_challenge,
                   optimal_rice_k)
from vectors import generate, verify_vectors


_scheme = None
_codec = None


def _get():
    global _scheme, _codec
    if _scheme is None:
        _scheme = LoTRS(TEST_PARAMS)
        _codec = LoTRSCodec(TEST_PARAMS)
    return _scheme, _codec


# ===========================================================================
#  Bitstream
# ===========================================================================

def test_bitwriter_reader_roundtrip():
    w = BitWriter()
    w.write_bits(0b10110, 5)
    w.write_bits(0xFF, 8)
    w.write_bits(0, 3)
    data = w.to_bytes()
    r = BitReader(data)
    assert r.read_bits(5) == 0b10110
    assert r.read_bits(8) == 0xFF
    assert r.read_bits(3) == 0


def test_bitreader_rejects_truncated():
    r = BitReader(b"\x00")
    r.read_bits(8)
    try:
        r.read_bits(1)
        assert False, "should have raised"
    except ValueError:
        pass


def test_bitreader_rejects_nonzero_padding():
    # 5 data bits + 3 padding bits.  Set a padding bit to 1.
    data = bytes([0b11111_111])        # all ones -> padding bits nonzero
    r = BitReader(data)
    r.read_bits(5)
    try:
        r.check_padding()
        assert False, "should reject nonzero padding"
    except ValueError:
        pass


# ===========================================================================
#  Fixed-width packing
# ===========================================================================

def test_fixed_roundtrip():
    coeffs = [0, 1, 100, 255, 128]
    enc = _pack_fixed(coeffs, 8)
    dec, _ = _unpack_fixed(enc, 5, 8)
    assert dec == coeffs


def test_fixed_signed_roundtrip():
    coeffs = [-50, 0, 50, -1, 1]
    offset = 50
    enc = _pack_fixed(coeffs, 7, offset)
    dec, _ = _unpack_fixed(enc, 5, 7, offset)
    assert dec == coeffs


def test_fixed_rejects_overflow():
    try:
        _pack_fixed([256], 8)
        assert False, "should reject"
    except ValueError:
        pass


# ===========================================================================
#  Rice coding
# ===========================================================================

def test_rice_roundtrip_zeros():
    coeffs = [0] * 32
    enc = _pack_rice(coeffs, 4, 100)
    dec, _ = _unpack_rice(enc, 32, 4, 100)
    assert dec == coeffs


def test_rice_roundtrip_random():
    rng = random.Random(99)
    for sigma in (5, 50, 500):
        rk = optimal_rice_k(sigma)
        bound = int(6 * sigma)
        coeffs = [max(-bound, min(bound, int(rng.gauss(0, sigma))))
                  for _ in range(64)]
        enc = _pack_rice(coeffs, rk, bound)
        dec, _ = _unpack_rice(enc, 64, rk, bound)
        assert dec == coeffs, f"Rice roundtrip failed for sigma={sigma}"


def test_rice_rejects_oob():
    try:
        _pack_rice([200], 4, 100)
        assert False, "should reject"
    except ValueError:
        pass


def test_rice_rejects_tampered():
    """Flip a bit and expect decode to reject or differ."""
    coeffs = [10, -20, 30, 0, -5]
    enc = bytearray(_pack_rice(coeffs, 3, 100))
    enc[1] ^= 0x80                         # flip a bit
    try:
        dec, _ = _unpack_rice(bytes(enc), 5, 3, 100)
        # if it decodes, it must differ
        assert dec != coeffs
    except ValueError:
        pass                               # rejection is also acceptable


def test_rice_compression_ratio():
    """Rice should beat fixed-width for Gaussian data."""
    rng = random.Random(42)
    sigma = 100
    d = 256
    rk = optimal_rice_k(sigma)
    bound = int(6 * sigma)
    coeffs = [max(-bound, min(bound, int(rng.gauss(0, sigma))))
              for _ in range(d)]
    rice_bytes = len(_pack_rice(coeffs, rk, bound))
    fixed_bits = (2 * bound + 1).bit_length()
    fixed_bytes = (d * fixed_bits + 7) // 8
    assert rice_bytes < fixed_bytes, \
        f"Rice ({rice_bytes}) should be smaller than fixed ({fixed_bytes})"


# ===========================================================================
#  Challenge encoding
# ===========================================================================

def test_challenge_roundtrip():
    ch = [0] * 64
    ch[3] = 1
    ch[7] = -1
    ch[31] = 1
    ch[60] = -1
    enc = _pack_challenge(ch, 4, 64)
    dec, _ = _unpack_challenge(enc, 4, 64)
    assert dec == ch


def test_challenge_rejects_unsorted():
    """Manually craft non-ascending positions."""
    w = BitWriter()
    pos_bits = 6                           # for d=64
    # write positions out of order: 10, 5
    w.write_bits(10, pos_bits)
    w.write_bits(5, pos_bits)
    w.write_bits(0, 1)                     # signs
    w.write_bits(0, 1)
    data = w.to_bytes()
    try:
        _unpack_challenge(data, 2, 64)
        assert False, "should reject non-ascending positions"
    except ValueError:
        pass


def test_challenge_rejects_out_of_range():
    """For non-power-of-2 d, positions can exceed d."""
    d = 50                                 # not a power of 2
    pos_bits = max(1, (d - 1).bit_length())  # 6 bits, max val 63
    wr = BitWriter()
    wr.write_bits(55, pos_bits)            # 55 >= d=50
    wr.write_bits(0, 1)                    # sign
    data = wr.to_bytes()
    try:
        _unpack_challenge(data, 1, d)
        assert False, "should reject position >= d"
    except ValueError:
        pass


# ===========================================================================
#  pp / sk / pk encoding
# ===========================================================================

def test_pp_roundtrip():
    _, codec = _get()
    pp = b"\x00" * 32
    assert codec.pp_decode(codec.pp_encode(pp)) == pp


def test_pp_rejects_bad_length():
    _, codec = _get()
    try:
        codec.pp_decode(b"\x00" * 31)
        assert False
    except ValueError:
        pass


def test_sk_roundtrip():
    _, codec = _get()
    seed = b"\xAB" * 32
    assert codec.sk_decode(codec.sk_encode(seed)) == seed


def test_pk_roundtrip():
    scheme, codec = _get()
    pp = scheme.setup(b"\x00" * 32)
    _, pk = scheme.keygen(pp, b"\x01" * 32)
    enc = codec.pk_encode(pk)
    dec = codec.pk_decode(enc)
    assert dec == pk


def test_pk_rejects_trailing():
    scheme, codec = _get()
    pp = scheme.setup(b"\x00" * 32)
    _, pk = scheme.keygen(pp, b"\x01" * 32)
    enc = codec.pk_encode(pk) + b"\x00"
    try:
        codec.pk_decode(enc)
        assert False, "should reject trailing bytes"
    except ValueError:
        pass


# ===========================================================================
#  Full signature round-trip
# ===========================================================================

def test_sig_roundtrip():
    """Encode -> decode -> verify."""
    scheme, codec = _get()
    par = TEST_PARAMS
    pp = scheme.setup(b"\x00" * 32)

    pk_table, all_sks = [], []
    for col in range(par.N):
        cpks, csks = [], []
        for row in range(par.T):
            sk, pk = scheme.keygen(pp, bytes([col, row]) + b"\x00" * 30)
            csks.append(sk); cpks.append(pk)
        pk_table.append(cpks); all_sks.append(csks)

    ell = 1
    mu = b"roundtrip"
    sig = scheme.sign(pp, all_sks[ell], ell, mu, pk_table, b"\xBB" * 32)

    enc = codec.sig_encode(sig)
    dec = codec.sig_decode(enc)
    assert scheme.verify(pp, mu, dec, pk_table), \
        "decoded signature must verify"


def test_sig_rejects_trailing():
    scheme, codec = _get()
    par = TEST_PARAMS
    pp = scheme.setup(b"\x00" * 32)

    pk_table, all_sks = [], []
    for col in range(par.N):
        cpks, csks = [], []
        for row in range(par.T):
            sk, pk = scheme.keygen(pp, bytes([col, row]) + b"\x00" * 30)
            csks.append(sk); cpks.append(pk)
        pk_table.append(cpks); all_sks.append(csks)

    ell = 0
    sig = scheme.sign(pp, all_sks[ell], ell, b"t", pk_table, b"\xCC" * 32)
    enc = codec.sig_encode(sig) + b"\xFF"
    try:
        codec.sig_decode(enc)
        assert False, "should reject trailing bytes"
    except ValueError:
        pass


def test_sig_rejects_truncated():
    scheme, codec = _get()
    par = TEST_PARAMS
    pp = scheme.setup(b"\x00" * 32)

    pk_table, all_sks = [], []
    for col in range(par.N):
        cpks, csks = [], []
        for row in range(par.T):
            sk, pk = scheme.keygen(pp, bytes([col, row]) + b"\x00" * 30)
            csks.append(sk); cpks.append(pk)
        pk_table.append(cpks); all_sks.append(csks)

    ell = 0
    sig = scheme.sign(pp, all_sks[ell], ell, b"t", pk_table, b"\xDD" * 32)
    enc = codec.sig_encode(sig)
    try:
        codec.sig_decode(enc[:len(enc) // 2])
        assert False, "should reject truncated"
    except ValueError:
        pass


def test_sig_encode_determinism():
    scheme, codec = _get()
    par = TEST_PARAMS
    pp = scheme.setup(b"\x00" * 32)

    pk_table, all_sks = [], []
    for col in range(par.N):
        cpks, csks = [], []
        for row in range(par.T):
            sk, pk = scheme.keygen(pp, bytes([col, row]) + b"\x00" * 30)
            csks.append(sk); cpks.append(pk)
        pk_table.append(cpks); all_sks.append(csks)

    ell = 2
    seed = b"\x55" * 32
    sig1 = scheme.sign(pp, all_sks[ell], ell, b"det", pk_table, seed)
    sig2 = scheme.sign(pp, all_sks[ell], ell, b"det", pk_table, seed)
    assert codec.sig_encode(sig1) == codec.sig_encode(sig2)


# ===========================================================================
#  Test vector generation and verification
# ===========================================================================

def test_vector_generation():
    vec = generate()
    assert vec["verification"] is True
    assert vec["signature"]["byte_length"] > 0
    assert len(vec["keygen"]) == TEST_PARAMS.N * TEST_PARAMS.T


def test_vector_self_verify():
    vec = generate()
    ok, errs = verify_vectors(vec)
    assert ok, f"self-verify failed: {errs}"


# ===========================================================================
#  Size reporting
# ===========================================================================

def test_sizes_nonnegative():
    _, codec = _get()
    s = codec.sizes()
    for name, val in s.items():
        assert val >= 0, f"{name} size must be non-negative"
    # w_tilde_hi is 0 for kappa=1 (new compression)
    if TEST_PARAMS.kappa == 1:
        assert s["w_tilde_hi"] == 0
    assert codec.total_size_estimate() == sum(s.values())


# ===========================================================================
#  Runner
# ===========================================================================

if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    print(f"Running {len(tests)} tests ...")
    for t in tests:
        t()
        print(f"  {t.__name__}: ok")
    print(f"test_codec.py: {len(tests)} tests passed")
