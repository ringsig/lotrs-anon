# LoTRS — Anonymous Artifact

This repository is the anonymous artifact for the paper **LoTRS:
Practical Post-Quantum Structured Threshold Ring Signatures from
Lattices**.  It contains two complete reference implementations
(Python and Rust), the parameter-estimation scripts that reproduce
the concrete parameter set, the FACCT-style sampler specification
used by both implementations, and the bench harness that produced
the timings reported in the paper.

```
README.md                  this file
LICENSE                    MIT
Makefile                   thin wrapper (top-level `clean`)
lotrs-facct-sampler.md     specification of the large-sigma sampler
lotrs-py/                  Python reference implementation
lotrs-rs/                  Rust performance implementation
estimator/                 SageMath parameter-estimation scripts
```

The Python implementation is the **golden reference** — it is short,
readable, and emits a deterministic test-vector blob that the Rust
implementation reproduces byte-for-byte.

## Quick reproducibility check

The fastest end-to-end smoke test:

```bash
# 1. Python reference: run all 137 unit / e2e tests + verify shipped test vectors.
cd lotrs-py
pip install -r requirements.txt              # numpy, pycryptodome, mpmath
for t in test_ring test_sample test_params \
         test_lotrs test_codec test_wtilde test_e2e; do
    python "$t.py"
done
python vectors.py --verify vectors.json
cd ..

# 2. Rust mirror: run all unit tests + 12 interop tests against the same vectors.json.
cd lotrs-rs
cargo test --release                          # finishes in ≈ 1 minute
```

The Rust interop tests load `../lotrs-py/vectors.json` and compare
`pp`, `sk`, `pk`, and full signatures byte-for-byte against the Python
reference.  All tests must pass.

## Reproducing the benchmarks

The numbers reported in the paper come from a single back-to-back
run of `examples/bench` on the Rust implementation.  To reproduce:

```bash
cd lotrs-rs
mkdir -p bench-out
cargo run --release --example bench --      \
    --grid "32,1:4:8:16:32;100,1:5:10:25:50;1,2:4:8:16" \
    > bench-out/grid-d128.md 2>&1
python3 scripts/plot_bench.py bench-out/grid-d128.md
```

Wall-clock is **≈ 15 minutes** on a single modern x86_64 core
(release build, single-threaded).  The plot script writes
`bench-out/{data.csv, summary.md, sign-vs-T.pdf, sigsize-vs-T.pdf,
breakdown-vs-T.pdf, rs-alone-vs-N.pdf, dualms-alone-vs-T.pdf}`.

For a faster smoke test that produces the headline parameter rows
without the full grid:

```bash
cargo run --release --example bench -- --skip-test --with-prod
```

The reference output captured by the authors lives under
`lotrs-rs/bench-out/` — `summary.md` is the human-readable digest,
`data.csv` is the machine-readable form.

The `bench-out/` files committed to the repository are the precise
outputs that produced the figures and tables in the paper.

See [`lotrs-rs/README.md`](lotrs-rs/README.md) § Benchmarks for the
hardware / methodology footnotes (sample counts per cell, the
empirical-vs-heuristic attempt-count gap, and column definitions).

## Reproducing the parameter selection

The concrete parameters in Table 3 of the paper are the output of
the SageMath scripts in `estimator/`:

```bash
cd estimator
make                # runs lotrs_estimate.py — prints the headline output
make reference      # writes a fresh comparison log to compare against
                    # the checked-in LoTRS-Estimate-Output-N100T50.txt
```

Headline values for the paper parameter point (`N = 100`, `T = 50`):

```
Signature size                            ≈ 35.14 KB
Single public key size                    ≈  7.13 KB
Ring PK size                              ≈ 35,625 KB
Number of repetitions for rejection samp. ≈ 12.34
```

Requires SageMath with Python support available as `sage`.  The
estimator vendors local copies of the lattice estimator and the
ASIS/MSIS estimator so no network access is needed at run time —
see [`estimator/README.md`](estimator/README.md) for provenance and
the full list of helper scripts.

## Component overviews

### `lotrs-py/`

Python reference implementation of the concrete `kappa = 1` LoTRS
instantiation.  Covers parameters / derived bounds, ring arithmetic
(with auxiliary-prime CRT-NTT), CDT and FACCT-style Gaussian
samplers, signature encoding/decoding, the full signing and
verification pipeline, unit tests, and the deterministic
test-vector emitter.  Pure Python; no compiled extensions.

### `lotrs-rs/`

Rust implementation of the same concrete LoTRS path, aligned with
the Python reference at the public wire-format boundary.  Covers
setup, key generation, aggregation, signing, and verification;
serialization compatible with the Python reference; a CRT-NTT ring
backend with pseudo-Mersenne reduction; CDT and FACCT-style
samplers; the Python/Rust interop tests; and the bench harness.

Both implementations precompute a 256-bit SHAKE-256 digest of the
canonical PK-table serialization once per signing or verification
call and feed *that* digest into `H_agg`, `H_com`, and the
Fiat-Shamir hash, instead of re-hashing the multi-MiB ring
public-key table at every call site.  Binding to the full PK is
preserved by SHAKE-256 collision-resistance.

### `estimator/`

SageMath scripts that reproduce the concrete parameter selection,
size calculations, and post-quantum security estimates.  Vendors
local copies of the lattice estimator and the ASIS/MSIS estimator
so the artifact can be evaluated without fetching dependencies.
See `estimator/README.md` for requirements and provenance notes for
the bundled external estimator components.

### `lotrs-facct-sampler.md`

Specification of the FACCT-style integer Gaussian sampler used for
the large masking distributions (`sigma_0`, `sigma_0_prime`) in the
benchmark and production parameter sets.  Covers the truncated
target distribution, the SHAKE-256/XOF-driven randomness, the
uniform proposal, the fixed-point Bernoulli-exp acceptance test,
the dispatch rule used by the implementations, and the validation
requirements that the Python/Rust cross-language KAT enforces.

## Artifact scope

The artifact focuses on the concrete `kappa = 1` implementation
path and the parameter set used for the paper's reported sizes and
benchmarks.  The Python implementation is the readable reference,
the Rust implementation is the performance-oriented implementation,
and the estimator reproduces the parameter and size calculations.

## License

MIT — see [`LICENSE`](LICENSE).
