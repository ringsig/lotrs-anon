# lotrs-rs

Rust implementation of **LoTRS** — practical post-quantum threshold ring
signatures from lattices.  Test-vector compatible with the Python reference
in [`../lotrs-py`](../lotrs-py/).

## Goals

- **Byte-for-byte interop with the Python reference.** Every seed-derived
  randomness stream, every serialized object, every Fiat–Shamir hash input
  must match `lotrs-py` exactly.  Integration tests load
  `../lotrs-py/vectors.json` and compare Rust-derived artefacts byte-wise.
- **Performance-oriented.** Uses the same CRT-NTT backend as `aux_ntt.py`,
  but with 64-bit native arithmetic, cache-friendly layouts, and
  pre-transformed CRT-NTT matrices for the large matrix-vector products.
  The auxiliary primes are both near `2^48`, so the Rust backend uses
  `u128` products plus pseudo-Mersenne reduction (`hi·c + lo mod p`,
  constants `c ∈ {16383, 19967}`) in the transform hot path.
- **More constant-time, but not CT-audited.** The arithmetic hot paths
  avoid coefficient-dependent branches for modular add/sub/neg,
  pseudo-Mersenne final reduction, CRT input splitting, and CRT output
  centering.  The signer as a whole is still **not** constant-time:
  Gaussian / challenge / uniform rejection loops, CDT binary search,
  FACCT acceptance, and restart logic remain data-dependent.  See
  *Constant-time considerations* below.
- **Robust verification.** `verify()` returns `false` for every malformed
  input class (wrong-length bytes, out-of-range coefficients, nonzero
  padding bits, trailing bytes, oversized Rice unary runs, ...) rather
  than panicking or surfacing error details that might leak signer state.
- **Fail-closed on unsupported parameter sets.** `LoTRS::try_new`
  refuses any parameter set it can't serve — currently `kappa != 1`
  or a sigma whose matching small-sigma CDT isn't shipped.

## Paper alignment

The public API mirrors the protocol figures of the accompanying paper:

| Paper            | Rust                                  |
|------------------|---------------------------------------|
| Fig. 4 `Setup`   | `LoTRS::setup`                        |
| Fig. 4 `KGen`    | `LoTRS::keygen`                       |
| Fig. 4 `KAgg`    | `LoTRS::kagg`                         |
| Fig. 4 `Sign₁`   | `LoTRS::sign1`                        |
| Fig. 4 `SAgg`    | `LoTRS::sagg`                         |
| Fig. 5 `Sign₂`   | `LoTRS::sign2`                        |
| Fig. 5 `Sign_bin`| `LoTRS::sign_bin`                     |
| Fig. 6 `Vf`      | `LoTRS::verify`                       |
| Fig. 1 `Rej` / `RejOp` | `sample::rej` / `sample::rej_op`|
| Table 2 bounds   | methods on `LoTRSParams`              |
| Table 3 params   | `PRODUCTION_PARAMS`                   |

Implementation notes for this artifact:

- The concrete artifact implements the `kappa = 1` parameter sets used by the
  paper.
- The ring backend uses exact CRT-NTT multiplication, matching the Python
  reference arithmetic while using Rust-native layouts.
- Large masking widths use the FACCT-style sampler described in
  [`../lotrs-facct-sampler.md`](../lotrs-facct-sampler.md).
- The structured public-key table is hashed once to a 256-bit digest and that
  digest is used in the `H_agg`, `H_com`, and Fiat-Shamir inputs.

## Crate layout

```
src/
  lib.rs         Public API (re-exports).
  params.rs      Parameter sets (TEST / BENCH_4OF32 / BENCH / PRODUCTION).
  aux_ntt.rs     Exact negacyclic multiplication via CRT over two
                 48-bit auxiliary primes  (mirrors ../lotrs-py/aux_ntt.py).
  ring.rs        Polynomial ring R_q = Z_q[X]/(X^d + 1).
  sample.rs      SHAKE256 XOF + uniform / short / ternary / challenge /
                 Gaussian samplers (matches ../lotrs-py/sample.py),
                 plus FACCT-style large-sigma sampler.
  codec.rs       Canonical serialization of pp / sk / pk / signature.
  cdt.rs         Precomputed CDT tables for every small-sigma width
                 used across TEST / BENCH_4OF32 / BENCH / PRODUCTION.  Ported from
                 ../lotrs-py via scripts/gen_cdt.py.
  lotrs.rs       Scheme implementation: Setup / KeyGen / Sign / Verify.
tests/
  interop.rs     Loads ../lotrs-py/vectors.json and checks byte-exact
                 reproduction of every wire object; Rust-vs-Python
                 signature bytes; tampering rejection; BENCH sign/verify.
  sampler_kat.rs Cross-language Gaussian-sampler KAT vs
                 tests/sampler_kat.json (emitted by the Python side).
scripts/
  gen_cdt.py         Regenerates src/cdt.rs from the Python reference.
  gen_sampler_kat.py Regenerates tests/sampler_kat.json.
```

## Status

| Pass | Deliverable                                             | Status  |
|------|---------------------------------------------------------|---------|
| 1    | `params` + `aux_ntt` + `ring` + `sample` + `codec`      | done ✅  |
| 2    | `lotrs::setup` / `keygen`, interop tests pp / sk / pk   | done ✅  |
| 3    | `lotrs::verify` + signature codec                       | done ✅  |
| 4    | CDT tables (TEST + `sigma_a` for `BENCH_4OF32` / `BENCH` / `PRODUCTION`) + generator | done ✅  |
| 5    | `lotrs::sign` (two-round, rejection sampling)           | done ✅  |
| 6    | FACCT-style large-sigma sampler + cross-language KAT    | done ✅  |
| 7    | Signing at `BENCH` / `PRODUCTION`                       | done ✅  |
| 8    | Constant-time hardening pass                            | partial: arithmetic hot paths hardened; signer not CT-audited |
| 9    | Benchmarks (`examples/bench.rs`)                        | partial: benchmark tool exists; full grid refresh deferred |

### Byte-for-byte interop vs `../lotrs-py/vectors.json`

Every externally-visible artefact of the Python reference is reproduced
bit-identically by the Rust code, under identical seeds:

* `pp_bytes`, `sk_bytes`, `pk_bytes` (all key pairs)
* full signature bytes: `LoTRS::sign(pp, sks, ell, mu, pk_table, seed)`
* verification:   a Python-generated signature verifies `true` in Rust;
  a Rust-generated signature verifies `true` in Rust (and would in Python)
* negative paths: tampered message / flipped bit / wrong pk-table shape
  / truncated / trailing bytes  all return `false` without panics

### Constant-time considerations (status)

The current implementation meets the **panic-free on any malformed input**
contract for `verify()` and the public codec entry points.  It has also
received a first arithmetic hardening pass:

* `Ring::{add,sub,neg}` and the in-place variants use mask-style
  conditional reductions rather than coefficient-dependent branches.
* Auxiliary-prime add/sub and pseudo-Mersenne final reduction use the
  same mask pattern.
* CRT input splitting no longer branches on `c > q/2`; it maps the
  negative centered case to `c + p_i - q` under a mask.
* CRT output centering and final non-negative reduction avoid
  branch-on-secret fixups.

This is **not** a full constant-time signer claim.  A deeper audit is
still outstanding; the concentrated risk areas are all on the
**signing** side:

1. **CDT binary search** in `xof_sample_gaussian` — branches on a
   secret-derived 128-bit value.  Leak channel: the sampled magnitude
   (which is secret), not the seed itself.  Mitigation for a hardened
   build: constant-time scan (linear or bit-sliced tournament) over the
   full table.  Cost scales with table length — tolerable for
   `sigma_a`, `sigma_b` (≤ ~20k entries); impractical at
   `sigma_0 ≈ 26M`, which is another reason a different sampler is
   needed at `BENCH` / `PRODUCTION` (see below).
2. **FACCT / rejection sampling** — FACCT proposals, Bernoulli
   accept/reject tests, and `rej` / `rej_op` all have data-dependent
   loop counts or branches.  Some accept / reject outcomes are already
   externally visible as signing restarts, but this is still not a
   constant-time sampler.
3. **Centred decomposition** in `Ring::centered_decompose` — branches
   on the signed coefficient value.  Inputs here are public (the
   commitments `w_tilde`, the proof commitments `B_bin`), so no secret
   data flows through the branch.
4. **Secret-key polynomials** (`sk_u`) are consumed through
   `Ring::vec_scale` and ring multiplication.  The CRT-NTT fast path
   now avoids the obvious coefficient-dependent arithmetic branches,
   but this is not a substitute for a compiler / microarchitecture CT
   audit.  The Python-matching schoolbook fallback at `d = 32` is for
   `TEST_PARAMS`, not a production path.

Bottom line: every public API is panic-free and returns `false` on
garbage input; arithmetic is now more CT-friendly; a full CT-audited
signer remains a future pass.

### Gaussian sampling

Two backends, matching the Python reference one-for-one:

* **CDT** (`xof_sample_gaussian`) — used for every width where the
  table fits in memory.  `src/cdt.rs` ships the eight CDTs needed by
  the four supported parameter sets: the four TEST widths
  (`sigma_0`, `sigma_0_prime`, `sigma_a`, `sigma_b`) plus
  `sigma_a` at `BENCH_4OF32` / `BENCH_PARAMS` / `PRODUCTION_PARAMS`
  and a single shared `sigma_b` for all three 4of32 / 16of32 / 50of100
  sets (they use the same `phi_b` and `B_b`).
* **FACCT-style** (`xof_sample_gaussian_facct` +
  `prepare_facct`) — used for the two mask widths `sigma_0`,
  `sigma_0_prime` at `BENCH_4OF32` / `BENCH_PARAMS` / `PRODUCTION_PARAMS`,
  where the CDT would need 10⁷–10⁸ entries.  Spec at
  [`../lotrs-facct-sampler.md`](../lotrs-facct-sampler.md);
  integer-only runtime, degree-20 Q(64) polynomial evaluated by
  Horner.

The dispatch is decided once in `LoTRS::try_new`, from the explicit
`mask_sampler` field on `LoTRSParams`: `sigma_a` / `sigma_b` always use
CDT; `sigma_0` / `sigma_0_prime` use the backend declared on the
parameter set (`Cdt` for `TEST_PARAMS`, `Facct` for `BENCH_4OF32` /
`BENCH_PARAMS` / `PRODUCTION_PARAMS`). The call sites in `sign1` /
`sign_bin` match on the pre-resolved sampler; there is no hot-path
dispatch.

**Cross-language sampler KAT.**  `scripts/gen_sampler_kat.py` emits
`tests/sampler_kat.json` from the Python reference (fixed seeds,
fixed sigmas, expected `i64[]`).  The integration test
`cross_language_sampler_kat_matches_python` asserts Rust reproduces
every entry byte-for-byte.  Covers:

* forced-FACCT at moderate sigma (CDT and FACCT would both work
  there, so we can compare)
* FACCT at `sigma ≈ 8.5 × 10⁶` (the `BENCH_PARAMS` `sigma_0` regime)
* CDT at small sigma

**BENCH / PRODUCTION end-to-end** signing + verification works
(see the `#[ignore]`d `bench_and_production_signing_round_trip`
test; ~1 min per signature in release mode).

Regenerate the CDT tables or KAT with

```bash
python scripts/gen_cdt.py            > src/cdt.rs
python scripts/gen_sampler_kat.py    > tests/sampler_kat.json
```

### Verified against `../lotrs-py/vectors.json`

* `pp_bytes`, `sk_bytes`, `pk_bytes` for every key pair — byte-for-byte.
* Signature bytes round-trip (decode then re-encode → identical).
* A Python-generated signature verifies `true` in Rust.
* Decoder refuses truncated, trailing-byte, and tampered signatures.
* Verifier returns `false` for a wrong message and for a flipped sig bit.

### CDT tables

At `TEST_PARAMS` all four Gaussian widths fit in CDT form and are
shipped as `const` arrays in `src/cdt.rs`.  At `BENCH_4OF32` /
`BENCH_PARAMS` / `PRODUCTION_PARAMS`, only the small binary-proof
widths `sigma_a` / `sigma_b` stay on CDT; the two large mask widths
`sigma_0` / `sigma_0_prime` use the FACCT-style sampler described
above.

Regenerate the shipped CDT tables with

```bash
python scripts/gen_cdt.py > src/cdt.rs
```

## Build

```bash
cd lotrs-rs
cargo build --release
cargo test
cargo test --test interop           # byte-exact vs lotrs-py/vectors.json
```

The interop test requires `../lotrs-py/vectors.json` to exist — generate
it by running `python vectors.py --out vectors.json` in `../lotrs-py/`.

## Benchmarks

`examples/bench.rs` measures the four high-level primitives
(`keygen` / `kagg` / `sign` / `verify`) and the wire sizes.  Two modes:

```bash
# Regression presets (TEST at d=32 + d=128 4-of-32, 16-of-32, optionally 50-of-100):
cargo run --release --example bench                  # TEST, 4-of-32, 16-of-32
cargo run --release --example bench -- --with-prod   # + 50-of-100

# Grid of (N, T) pairs on the d=128 lattice (full grid for the
# captured data and plots in bench-out/, including the standalone
# RS (T=1) and DualMS (N=1) sweeps):
cargo run --release --example bench -- \
    --grid "32,1:4:8:16:32;100,1:5:10:25:50;1,2:4:8:16"
```

The reproducible-grid command above completes in ≈ **15 minutes** on
the hardware below; it produces the table shown under
[Grid results](#grid-results) and corresponds to the captured data in
[`bench-out/`](bench-out/).

### Plotting bench output

`scripts/plot_bench.py` parses the bench binary's Markdown report and
emits CSV + paper-ready plots. Capture the bench output to a file and
hand it to the script:

```bash
mkdir -p bench-out
cargo run --release --example bench -- \
    --grid "32,1:4:8:16:32;100,1:5:10:25:50;1,2:4:8:16" \
    > bench-out/grid-d128.md 2>&1
python3 scripts/plot_bench.py bench-out/grid-d128.md
# → bench-out/{data.csv, summary.md,
#              sign-vs-T.pdf, sigsize-vs-T.pdf, breakdown-vs-T.pdf,
#              rs-alone-vs-N.pdf, dualms-alone-vs-T.pdf}
```

The grid has three regimes mixed into one run:

* **threshold** rows (`N≥2, T≥2`): the LoTRS combined protocol.
* **RS-alone** rows (`T=1, N≥2`): plain ring signature, one signer.
* **DualMS-alone** rows (`N=1, T≥2`): plain multi-sig, no ring hiding.

Outputs land alongside the input file:

| file | content |
|---|---|
| `data.csv`              | one row per `(N, T)` cell with every timing column + every byte size |
| `sign-vs-T.pdf`         | log-log Sign / Verify time vs threshold `T`, one line per `N` (threshold rows only) |
| `breakdown-vs-T.pdf`    | Sign decomposition: total / DualMS / RS sub-times vs `T` |
| `sigsize-vs-T.pdf`      | signature size (KiB) vs `T` (threshold rows only) |
| `rs-alone-vs-N.pdf`     | RS-alone (`T=1`): Sign / Verify / sig size vs ring size `N` |
| `dualms-alone-vs-T.pdf` | DualMS-alone (`N=1`): Sign$_{\text{DualMS}}$ / Verify$_{\text{DualMS}}$ vs signer count `T` |
| `summary.md`            | Markdown tables, one per regime (paper-friendly subsets of the columns) |

Requires `pip install matplotlib`. Running with no path argument reads
the report from stdin and writes outputs to the current working directory.

If `examples/bench.rs` gains or removes columns in any of its three
output tables (`Primitive timings` / `Sign / Verify breakdown` /
`Key / signature sizes`), the regexes at the top of `plot_bench.py`
must be updated to match.

### Parameters used for the numbers below

All rows use the `d=128` lattice from `estimator/lotrs_estimate.py`
(2026-04-19 revision):

```
κ = 1,  k = 12,  l = 5,  l' = 6,  n̂ = 10,  k̂ = 8,  w = 31,  η = 1
q     = 274877906837   (largest prime ≤ 2^38 with q ≡ 5 mod 8)
q_hat = 8589934237     (largest prime ≤ 2^33 with q_hat ≡ 5 mod 8)
φ_a = 50,  φ_b = 4,  φ = 11.75·T    (rejection-sampling slack)
K_A = 20,  K_B = 5,  K_w = 5,  ε_tot = 0.01
mask_sampler = FACCT,  tail_t = 1.2
```

`N = β^κ = β` (κ = 1), so only `N`, `T`, and `φ = 11.75·T` vary across
rows; everything else is fixed.

### Methodology

**Sign** and **Verify** are arithmetic means over multiple signing
seeds — the rejection-sampling attempt count per signing call is
geometrically distributed.  The estimator's heuristic predicts
`μ_total ≈ 12.3` accepted-attempt mean; the implementation also
restarts on the additional `w̃₀`-stability check used by the κ=1
commitment-compression path, so empirical means are in the 25–40
range across the grid.  Per-attempt sign cost is the more
parameter-stable indicator.  Sample counts per cell: **12** at T≤5,
**10** at T≤10, **8** at T≤25, **5** at T≥32.  The residual
cell-to-cell sampling noise is ≈ ±10% at these counts.

**KeyGen** and **KAgg** are single-run (deterministic given `pp` / PK).

### Hardware / toolchain

* **CPU**: AMD Ryzen AI 9 HX 370 (Zen 5, 12 cores / 24 threads,
  5.1 GHz max boost)
* **OS**: Debian forky/sid (`Linux 6.19.13-1 amd64`)
* **Toolchain**: `rustc 1.95.0`, single-threaded

A production deployment would parallelise the per-signer loops and
amortise the CRT-NTT transforms; the numbers here are a
single-threaded baseline.

### Grid results

Single grid covering both ring sizes (N = 32, 100) with the same
`d = 128` lattice, plus standalone RS (T=1) and DualMS (N=1) sweeps.
Captured in
[`bench-out/grid-d128.md`](bench-out/grid-d128.md) (raw report),
[`bench-out/data.csv`](bench-out/data.csv),
[`bench-out/summary.md`](bench-out/summary.md), and the five paper
plots: [`sign-vs-T.pdf`](bench-out/sign-vs-T.pdf),
[`sigsize-vs-T.pdf`](bench-out/sigsize-vs-T.pdf),
[`breakdown-vs-T.pdf`](bench-out/breakdown-vs-T.pdf),
[`rs-alone-vs-N.pdf`](bench-out/rs-alone-vs-N.pdf),
[`dualms-alone-vs-T.pdf`](bench-out/dualms-alone-vs-T.pdf).

|   N |  T | Sign       | Verify   | KAgg     |  signature  | single pk |  ring PK   |
|  -: | -: | -:         | -:       | -:       |   -:        |   -:      |    -:      |
|  32 |  4 |   1.38 s   |   60 ms  |   35 ms  |  22.53 KiB  | 7.12 KiB  | 0.89 MiB   |
|  32 |  8 |   4.29 s   |  100 ms  |   71 ms  |  22.99 KiB  | 7.12 KiB  | 1.78 MiB   |
|  32 | 16 |   6.01 s   |  177 ms  |  144 ms  |  23.62 KiB  | 7.12 KiB  | 3.56 MiB   |
|  32 | 32 |   6.05 s   |  334 ms  |  288 ms  |  24.07 KiB  | 7.12 KiB  | 7.12 MiB   |
| 100 |  5 |   4.60 s   |  208 ms  |  134 ms  |  33.74 KiB  | 7.12 KiB  | 3.48 MiB   |
| 100 | 10 |  12.40 s   |  364 ms  |  276 ms  |  34.36 KiB  | 7.12 KiB  | 6.96 MiB   |
| 100 | 25 |  29.24 s   |  828 ms  |  702 ms  |  35.09 KiB  | 7.12 KiB  | 17.40 MiB  |
| 100 | 50 | **50.39 s**| **1.61 s** | **1.41 s** | **35.53 KiB** | 7.12 KiB | 34.79 MiB  |

The `(N, T) = (100, 50)` row is Table 3 of the paper (= the
`PRODUCTION_PARAMS` preset). The four `(32, *)` rows include the
named `BENCH_4OF32` (T=4) and `BENCH_PARAMS` (T=16) presets used by
the regression suite. Total grid wall-clock ≈ 15 minutes on the
hardware listed above.

`KeyGen` sits at ~1.6 ms across every row (independent of `N`, `T`).
`KAgg` scales as `N·T` (deterministic, single-run): 5000 pk products
at `(100, 50)` → 1.41 s. `pk` is 7.125 KiB in every row; the
`N·T·pk` ring-PK table scales accordingly. Sign time grows
near-linearly in `T` at fixed `N` and grows by roughly the `N` ratio
at fixed `T` (compare `(32, T)` vs `(100, T)` rows).

Implementation note: the structured PK table is pre-hashed once as a
256-bit SHAKE256 digest before the protocol-specific hash calls. The
digest, not the full table serialization, is then fed to `H_agg`,
`H_com`, and the Fiat-Shamir hash. This keeps the transcript bound to
the full PK table while avoiding repeated hashing of the multi-MiB
`N·T·pk` input.

Signature size grows only ~6 % from T=5 to T=50 at fixed `N` — the
T-dependent component is small versus the binary-proof / `B_bin`
terms, which scale with `N`.

### Smoke-test presets

For quick regression runs without the full grid, the named
parameter sets `BENCH_4OF32` (T=4, N=32), `BENCH_PARAMS` (T=16, N=32),
and `PRODUCTION_PARAMS` (T=50, N=100) reproduce three rows of the
grid above:

```bash
cargo run --release --example bench -- --skip-test            # 4-of-32 + 16-of-32 (< 1 min)
cargo run --release --example bench -- --skip-test --with-prod # + 50-of-100 (~ 5 min)
```

The `BENCH_PARAMS` (16-of-32) row is the canonical small-bench point.

Observations:

* **KeyGen** is independent of `N` and `T` — it's one `mat_vec`
  `A·s` in `R_q` at fixed `(k, l)`.  The ~1 ms baseline sits where
  you'd expect for `k·(l+k) = 12·17 = 204` ring multiplications at
  `d=128`.
* **KAgg** computes `pk_hash = H(PK)` once, then runs `T` small
  `hash_agg(pk_hash, u)` calls plus `N·T` ring mul-adds. The dominant
  deterministic cost is therefore the `N·T` arithmetic, not repeated
  hashing of the full PK table.
* **Sign** is dominated by rejection-sampling attempts.  The
  estimator's `μ_total = μ·μ_a·μ_b·μ_BG·μ_fg ≈ 12.3` is the
  multiplicative restart factor of the four samplers and is roughly
  constant across `T` because `φ = 11.75·T` holds `μ^T` constant; the
  implementation also restarts on the additional `w̃₀`-stability
  check, so the empirical mean attempts in this grid are 25–40 (see
  the `att.` column of `data.csv`).  Per-attempt cost scales
  ~linearly in `T·N` for the signer / kagg / sign_bin loops.
* **Verify** tracks `kagg` closely — the bulk of the verifier cost
  is the same `N·T` product.
* **Signature size** grows as `√T` at fixed `N` (Rice-coded
  `z̃`, `r̃`, `ẽ` have widths `∝ √T·σ₀`).  Moving from `N=25` to
  `N=100` at fixed `T` adds ~55% because `f_1` carries
  `N·(β−1)·d` = `N·(N−1)·d` coefficients (quadratic in `N`).
* **pk** is the same 7.125 KiB everywhere — a `k`-poly vector in
  `R_q` at `d=128, k=12, log q = 38`, giving
  `k·d·⌈log₂ q⌉ / 8 = 12·128·38/8 = 7.125 KiB`.
* Sizes agree with the `estimator/lotrs_estimate.py` per-component
  totals to within codec overhead.  At `(N, T) = (100, 50)` the
  estimator predicts 35.14 KB and the codec emits 35.53 KiB; the
  small surplus is the Golomb-Rice header / per-stream parameter
  overhead, which the estimator ignores.
