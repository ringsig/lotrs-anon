"""
vectors.py -- Deterministic test-vector generation for LoTRS.

Runs a fully pinned keygen -> sign -> serialize -> deserialize -> verify
workflow and emits a JSON blob with hex-encoded binary for every object
and selected intermediate values.  Useful for cross-implementation
comparison.

Fixed seeds
-----------
  pp seed          : 00 01 02 ... 1f
  signer (col, row): bytes([col ^ 0x40, row ^ 0x80]) + b'\\x00' * 30
  signing seed     : aa aa ... aa  (32 bytes)

Usage
-----
  python vectors.py                      # print JSON to stdout
  python vectors.py --out vectors.json   # write to file
  python vectors.py --verify vectors.json  # re-derive and compare
"""

import json
import sys
import math

from params import TEST_PARAMS
from lotrs import LoTRS
from codec import LoTRSCodec


# ---- fixed seeds ---------------------------------------------------------

PP_SEED = bytes(range(32))                            # 00..1f

def signer_seed(col, row):
    return bytes([col ^ 0x40, row ^ 0x80]) + b"\x00" * 30

SIGNING_SEED = b"\xAA" * 32


# ---- helpers -------------------------------------------------------------

def _hex(data):
    if isinstance(data, (bytes, bytearray)):
        return data.hex()
    return data


def _poly_hex(poly, ring):
    """Hex string of centred coefficients (compact, for JSON)."""
    return [int(c) for c in ring.centered(poly)]


# ---- generation ----------------------------------------------------------

def generate(par=None):
    """
    Generate a complete test-vector dict.

    Returns a JSON-serialisable dict.
    """
    if par is None:
        par = TEST_PARAMS

    scheme = LoTRS(par)
    codec = LoTRSCodec(par)
    Rq = scheme.Rq

    pp = scheme.setup(PP_SEED)
    pp_bytes = codec.pp_encode(pp)

    # ---- keygen ----------------------------------------------------------
    N, T = par.N, par.T
    pk_table = []
    all_sks = []
    keygen_records = []

    for col in range(N):
        col_pks, col_sks = [], []
        for row in range(T):
            seed = signer_seed(col, row)
            sk, pk = scheme.keygen(pp, seed)
            col_sks.append(sk)
            col_pks.append(pk)
            keygen_records.append({
                "col": col,
                "row": row,
                "seed": _hex(seed),
                "sk_bytes": _hex(codec.sk_encode(seed)),
                "pk_bytes": _hex(codec.pk_encode(pk)),
            })
        pk_table.append(col_pks)
        all_sks.append(col_sks)

    # ---- sign ------------------------------------------------------------
    ell = 1
    mu = b"test vector"
    sks = all_sks[ell]

    sig = scheme.sign(pp, sks, ell, mu, pk_table, SIGNING_SEED)

    # ---- serialize -------------------------------------------------------
    sig_bytes = codec.sig_encode(sig)

    # ---- deserialize and verify ------------------------------------------
    sig_decoded = codec.sig_decode(sig_bytes)
    ok = scheme.verify(pp, mu, sig_decoded, pk_table)

    # ---- intermediate values ---------------------------------------------
    x_centered = Rq.centered(sig["pi"]["x"])

    vec = {
        "schema_version": 1,
        "params": par.name,
        "d": par.d,
        "N": N,
        "T": T,
        "kappa": par.kappa,
        "beta": par.beta,
        "q": par.q,
        "q_hat": par.q_hat,
        "pp_seed": _hex(PP_SEED),
        "pp_bytes": _hex(pp_bytes),
        "signing_seed": _hex(SIGNING_SEED),
        "ell": ell,
        "message": mu.decode(),
        "keygen": keygen_records,
        "signature": {
            "bytes": _hex(sig_bytes),
            "byte_length": len(sig_bytes),
            "x": x_centered,
            "z_tilde_inf_norm": Rq.vec_inf_norm(sig["z_tilde"]),
            "r_tilde_inf_norm": Rq.vec_inf_norm(sig["r_tilde"]),
            "e_tilde_inf_norm": Rq.vec_inf_norm(sig["e_tilde"]),
        },
        "verification": ok,
        "sizes": codec.sizes(),
        "encoding": {
            "rice_k_f1": codec.rice_f1,
            "rice_k_zb": codec.rice_zb,
            "rice_k_zt": codec.rice_zt,
            "rice_k_rt": codec.rice_rt,
            "rice_k_et": codec.rice_et,
            "dx_bbin": codec.dx_bbin,
            "dx_whi": codec.dx_whi,
            "dx_pk": codec.dx_pk,
        },
    }
    return vec


# ---- verification --------------------------------------------------------

def verify_vectors(vec):
    """
    Re-derive everything from seeds and compare against the stored
    test vector.  Returns (ok: bool, errors: list[str]).
    """
    par = TEST_PARAMS
    if vec["params"] != par.name:
        return False, [f"params mismatch: {vec['params']} vs {par.name}"]

    errors = []
    scheme = LoTRS(par)
    codec = LoTRSCodec(par)

    pp = scheme.setup(bytes.fromhex(vec["pp_seed"]))

    # ---- pp round-trip ---------------------------------------------------
    pp_bytes = codec.pp_encode(pp)
    if "pp_bytes" in vec and pp_bytes.hex() != vec["pp_bytes"]:
        errors.append("pp_bytes mismatch")
    pp_roundtrip = codec.pp_encode(codec.pp_decode(pp_bytes))
    if pp_roundtrip != pp_bytes:
        errors.append("pp decode/encode is not canonical")

    # ---- re-derive keys --------------------------------------------------
    N, T = par.N, par.T
    pk_table = []
    all_sks = []
    for col in range(N):
        col_pks, col_sks = [], []
        for row in range(T):
            rec = vec["keygen"][col * T + row]
            seed = bytes.fromhex(rec["seed"])
            sk, pk = scheme.keygen(pp, seed)
            col_sks.append(sk)
            col_pks.append(pk)

            actual_pk_hex = codec.pk_encode(pk).hex()
            if actual_pk_hex != rec["pk_bytes"]:
                errors.append(
                    f"pk_bytes mismatch at col={col} row={row}")
            pk_roundtrip = codec.pk_encode(
                codec.pk_decode(bytes.fromhex(rec["pk_bytes"])))
            if pk_roundtrip.hex() != rec["pk_bytes"]:
                errors.append(
                    f"pk decode/encode not canonical at col={col} row={row}")

            if "sk_bytes" in rec:
                actual_sk_hex = codec.sk_encode(seed).hex()
                if actual_sk_hex != rec["sk_bytes"]:
                    errors.append(
                        f"sk_bytes mismatch at col={col} row={row}")
                sk_roundtrip = codec.sk_encode(
                    codec.sk_decode(bytes.fromhex(rec["sk_bytes"])))
                if sk_roundtrip.hex() != rec["sk_bytes"]:
                    errors.append(
                        f"sk decode/encode not canonical at col={col} row={row}")
        pk_table.append(col_pks)
        all_sks.append(col_sks)

    # ---- re-derive signature ---------------------------------------------
    ell = vec["ell"]
    mu = vec["message"].encode()
    signing_seed = bytes.fromhex(vec["signing_seed"])
    sks = all_sks[ell]

    sig = scheme.sign(pp, sks, ell, mu, pk_table, signing_seed)
    sig_bytes = codec.sig_encode(sig)

    expected_hex = vec["signature"]["bytes"]
    actual_hex = sig_bytes.hex()
    if actual_hex != expected_hex:
        errors.append("signature bytes mismatch")
        # find first divergence byte for diagnostics
        for i in range(min(len(actual_hex), len(expected_hex)) // 2):
            a = actual_hex[2 * i:2 * i + 2]
            e = expected_hex[2 * i:2 * i + 2]
            if a != e:
                errors.append(f"  first diff at byte {i}: "
                              f"got 0x{a} expected 0x{e}")
                break

    if len(sig_bytes) != vec["signature"]["byte_length"]:
        errors.append(
            f"sig length {len(sig_bytes)} vs "
            f"expected {vec['signature']['byte_length']}")

    # ---- decode and verify -----------------------------------------------
    try:
        sig_dec = codec.sig_decode(sig_bytes)
    except ValueError as exc:
        errors.append(f"sig_decode failed: {exc}")
        return len(errors) == 0, errors

    sig_roundtrip = codec.sig_encode(sig_dec)
    if sig_roundtrip != sig_bytes:
        errors.append("sig decode/encode is not canonical")

    ok = scheme.verify(pp, mu, sig_dec, pk_table)
    if not ok:
        errors.append("verification failed on re-derived signature")

    return len(errors) == 0, errors


# ---- CLI -----------------------------------------------------------------

def main():
    if len(sys.argv) >= 3 and sys.argv[1] == "--verify":
        path = sys.argv[2]
        with open(path) as f:
            vec = json.load(f)
        print(f"Verifying test vectors from {path} ...")
        ok, errs = verify_vectors(vec)
        if ok:
            print("ALL CHECKS PASSED")
        else:
            print(f"FAILED ({len(errs)} errors):")
            for e in errs:
                print(f"  {e}")
        sys.exit(0 if ok else 1)

    print("Generating test vectors ...", file=sys.stderr)
    vec = generate()

    if len(sys.argv) >= 3 and sys.argv[1] == "--out":
        path = sys.argv[2]
        with open(path, "w") as f:
            json.dump(vec, f, indent=2)
        print(f"Written to {path}", file=sys.stderr)
    else:
        json.dump(vec, sys.stdout, indent=2)
        print(file=sys.stdout)


if __name__ == "__main__":
    main()
