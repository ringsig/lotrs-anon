#!/usr/bin/env python3
"""
gen_sampler_kat.py -- Emit cross-language KAT vectors for the Gaussian samplers.

Produces `tests/sampler_kat.json` in the Rust crate, with entries that
pin the byte-for-byte behaviour of each backend:

- FACCT at moderate sigma (forced path)
- FACCT at large sigma
- CDT at small sigma

Every entry carries:

    { "name":       str,
      "backend":    "facct" | "cdt",
      "sigma":      float,
      "d":          int,
      "seed":       hex,
      "tags":       [ hex-encoded bytes, ... ],
      "samples":    [ i64, ... ]                  # length d
    }

The Rust integration test reconstructs the same XOF (seed || tags),
runs the matching backend at the same sigma, and asserts the sample
list matches exactly.

Usage
-----
    python scripts/gen_sampler_kat.py  > tests/sampler_kat.json
"""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LOTRS_PY = os.path.normpath(os.path.join(HERE, "..", "..", "lotrs-py"))
if LOTRS_PY not in sys.path:
    sys.path.insert(0, LOTRS_PY)

from params import TEST_PARAMS, BENCH_4OF32, BENCH_PARAMS, PRODUCTION_PARAMS  # noqa
from sample import (                                                # noqa
    make_xof,
    build_cdt_sampler,
    build_facct_sampler,
    xof_sample_gaussian,
    xof_sample_gaussian_facct,
)


def run_facct(name, sigma, d, seed, tag):
    sampler = build_facct_sampler(sigma, lam=128)
    xof = make_xof(seed, tag)
    samples = xof_sample_gaussian_facct(xof, sampler.data, d)
    return {
        "name": name, "backend": "facct",
        "sigma": float(sigma), "d": d,
        "seed": seed.hex(), "tags": [tag.hex()],
        "samples": samples,
    }


def run_cdt(name, sigma, d, seed, tag, cdt_name=None):
    """Build a CDT sampler for this sigma.  If `cdt_name` is provided
    it names a `pub static` shipped in the Rust `cdt` module, letting
    the Rust KAT fetch the exact table instead of inferring one from
    sigma; otherwise the entry is informational / Python-only."""
    sampler = build_cdt_sampler(sigma, lam=128)
    xof = make_xof(seed, tag)
    samples = xof_sample_gaussian(xof, sampler.data, sampler.lam, d)
    entry = {
        "name": name, "backend": "cdt",
        "sigma": float(sigma), "d": d,
        "seed": seed.hex(), "tags": [tag.hex()],
        "samples": samples,
    }
    if cdt_name is not None:
        entry["cdt_name"] = cdt_name
    return entry


def main():
    vectors = [
        # FACCT at moderate sigma where CDT would also be feasible
        # (two separate seeds so the sampler sweep is non-trivial).
        run_facct("facct_moderate_seed0", sigma=100.0, d=256,
                  seed=bytes([0x10] * 32), tag=b"KAT-facct-0"),
        run_facct("facct_moderate_seed1", sigma=250.5, d=128,
                  seed=bytes(range(32)),  tag=b"KAT-facct-1"),
        # FACCT at the BENCH / PRODUCTION sigma_0 regime.
        run_facct("facct_large_bench",    sigma=8.5e6, d=128,
                  seed=b"\xaa" * 32,      tag=b"KAT-facct-bench"),
        # CDT at small sigma.  The `cdt_name` tells the Rust side
        # which shipped `pub static` to fetch; without it the CDT
        # fixtures would only be covered indirectly via the full-
        # signature interop tests.
        run_cdt("cdt_test_sigma_a",
                sigma=TEST_PARAMS.sigma_a, d=128,
                seed=b"\x0c" * 32, tag=b"KAT-cdt-test-a",
                cdt_name="CDT_SIGMA_A_TEST"),
        run_cdt("cdt_test_sigma_b",
                sigma=TEST_PARAMS.sigma_b, d=64,
                seed=b"\x0d" * 32, tag=b"KAT-cdt-test-b",
                cdt_name="CDT_SIGMA_B_TEST"),
        run_cdt("cdt_bench_4of32_sigma_a",
                sigma=BENCH_4OF32.sigma_a, d=64,
                seed=b"\x1a" * 32, tag=b"KAT-cdt-bench-4of32-a",
                cdt_name="CDT_SIGMA_A_BENCH_4OF32"),
        run_cdt("cdt_bench_sigma_a",
                sigma=BENCH_PARAMS.sigma_a, d=64,
                seed=b"\x0e" * 32, tag=b"KAT-cdt-bench-a",
                cdt_name="CDT_SIGMA_A_BENCH"),
        run_cdt("cdt_bench_sigma_b",
                sigma=BENCH_PARAMS.sigma_b, d=64,
                seed=b"\x0f" * 32, tag=b"KAT-cdt-bench-b",
                cdt_name="CDT_SIGMA_B_BENCH_PRODUCTION"),
        run_cdt("cdt_production_sigma_a",
                sigma=PRODUCTION_PARAMS.sigma_a, d=64,
                seed=b"\x11" * 32, tag=b"KAT-cdt-prod-a",
                cdt_name="CDT_SIGMA_A_PRODUCTION"),
    ]
    json.dump({"kat": vectors}, sys.stdout, indent=2)
    print()


if __name__ == "__main__":
    main()
