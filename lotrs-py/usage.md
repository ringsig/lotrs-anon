# LoTRS Python reference — usage

Reference implementation and test-vector generator for the LoTRS lattice-based structured threshold ring signature scheme (kappa = 1 case). Run from the `lotrs-py/` directory.

## Prerequisites

```bash
pip install -r requirements.txt   # numpy, pycryptodome, mpmath
```

Python 3.10+. No compiled extensions.

## Running tests

```bash
python test_ring.py       # ring arithmetic + CRT-NTT (22 tests)
python test_sample.py     # XOF samplers, CDT + FACCT, rejection sampling (26 tests)
python test_params.py     # parameter consistency (24 tests)
python test_lotrs.py      # scheme unit + end-to-end tests (10 tests)
python test_codec.py      # serialization, Rice coding, test vectors (26 tests)
python test_wtilde.py     # w-tilde compression tests (16 tests)
python test_e2e.py        # full pipeline and tampering tests (13 tests)
```

137 tests total across 7 files. First run takes a few seconds to build the small-sigma CDTs; subsequent scheme operations reuse the cached tables.  Rust mirror at `../lotrs-rs/` adds 48 unit + 12 interop + 1 sampler KAT — see [`lotrs-rs/README.md`](../lotrs-rs/README.md).

### Module self-tests

Each module has an `if __name__ == "__main__"` block with basic smoke tests:

```bash
python ring.py            # negacyclic property, commutativity, decomposition
python aux_ntt.py         # auxiliary-prime NTT, CRT reconstruction vs schoolbook (d=128 and d=256)
python sample.py          # determinism, Gaussian tail, challenge weight
python params.py          # prints derived parameter values for TEST, BENCH, and PRODUCTION
python lotrs.py           # full keygen → sign → verify cycle
python codec.py           # Rice round-trip, challenge round-trip, size report
```

## Programmatic usage

### End-to-end signing

```python
from params import TEST_PARAMS
from lotrs import LoTRS

scheme = LoTRS(TEST_PARAMS)              # builds CDT tables (~2s)
pp = scheme.setup(b"\x00" * 32)          # public parameters

# Key generation: N columns, T rows
pk_table = []
all_sks = []
for col in range(TEST_PARAMS.N):         # N = 4
    col_pks, col_sks = [], []
    for row in range(TEST_PARAMS.T):     # T = 2
        seed = bytes([col, row]) + b"\x00" * 30
        sk, pk = scheme.keygen(pp, seed)
        col_sks.append(sk)
        col_pks.append(pk)
    pk_table.append(col_pks)
    all_sks.append(col_sks)

# Sign with the signers from column ell
ell = 1
sks = all_sks[ell]
mu = b"approve transaction"

sig = scheme.sign(pp, sks, ell, mu, pk_table, b"\xAA" * 32)

# Verify
assert scheme.verify(pp, mu, sig, pk_table)
```

### Step-by-step signing (two-round protocol)

For applications that need to inspect intermediate state or simulate the interactive protocol across a network:

```python
from sample import make_xof

rho = make_xof(signing_seed, b"rho", attempt).read(32)

# Round 1: each signer produces commitments
states, all_coms = [], []
for u in range(T):
    st, com = scheme.sign1(pp, sks[u], u, ell, mu, pk_table, rho, attempt)
    states.append(st)
    all_coms.append(com)

# Round 2: each signer produces a partial signature
sigmas = []
for u in range(T):
    sig_u = scheme.sign2(states[u], all_coms, pk_table)
    if sig_u is None:
        break                             # restart with next attempt
    sigmas.append(sig_u)

# Aggregate
sigma = scheme.sagg(sigmas)
```

### Serialization

```python
from codec import LoTRSCodec

codec = LoTRSCodec(TEST_PARAMS)

# Encode / decode public key
pk_bytes = codec.pk_encode(pk)            # fixed-width at log_q bits/coeff
pk_back = codec.pk_decode(pk_bytes)       # raises ValueError on malformation

# Encode / decode signature (Rice coding for Gaussian components)
sig_bytes = codec.sig_encode(sig)
sig_back = codec.sig_decode(sig_bytes)    # raises ValueError on malformation
assert scheme.verify(pp, mu, sig_back, pk_table)

# Size breakdown
codec.print_sizes()
```

Decoders reject: out-of-range coefficients, nonzero padding bits, trailing bytes, truncated data, and excessive Rice unary runs (DoS guard).

### Test vector generation

```bash
python vectors.py --out vectors.json          # generate from fixed seeds
python vectors.py --verify vectors.json       # re-derive and compare byte-for-byte
```

All randomness derives from deterministic seeds. The `sign()` convenience method handles the restart loop deterministically: attempt `i` derives `rho_i = SHAKE256(signing_seed || "rho" || i)`.

Fixed seeds used by `vectors.py`:
- pp seed: `00 01 02 ... 1f`
- signer (col, row): `bytes([col ^ 0x40, row ^ 0x80]) + b'\x00' * 30`
- signing seed: `aa aa ... aa`

The JSON blob contains hex-encoded binary for pp, pk, and the full signature, plus intermediate values (challenge x, norms) for cross-implementation debugging.

## Parameter sets

### TEST_PARAMS (correctness testing)

Small, fast, not secure. Uses distinct primes q != q_hat to exercise both rings.

| Parameter | Value | Notes |
|-----------|-------|-------|
| d | 32 | Ring dimension |
| q | 4194389 | Prime, 5 mod 8, ~2^22 |
| q_hat | 7000061 | Distinct prime, 5 mod 8, ~2^23 |
| kappa | 1 | Only supported value |
| beta | 4 | Column index base |
| N | 4 | Ring size (columns in PK) |
| T | 2 | Threshold (signers) |
| k, l, l' | 2, 2, 3 | Matrix dimensions |
| n_hat, k_hat | 2, 3 | Binary proof dimensions |
| w | 4 | Challenge weight |
| eta | 1 | Ternary secret keys |
| phi, phi_a, phi_b | 12.0, 12.0, 12.0 | Rejection sampling slack |
| K_A, K_B, K_w | 13, 4, 5 | Bit-dropping |
| tail_t | 2.0 | Gaussian tail factor (looser than 1.2 for small dims) |
| `mask_sampler` | `"cdt"` (σ₀ ≈ 2.2 × 10³) | Small enough for a CDT |

Expected ~33 signing attempts per signature (μ·μ_a·μ_b·μ_BG·μ_fg at the estimator-formula aggregate, `μ(12)² ≈ 7.4` for the T=2 z_u factor alone).

### BENCH_4OF32 (4-of-32 benchmark-only variant)

The three d=128 sets all share the same lattice (`k=12, l=5, l'=6, n̂=10, k̂=8`, `q = 274877906837`, `q_hat = 8589934237`, `phi_a=50, phi_b=4`, `K_A=20, K_B=5, K_w=5`, `mask_sampler="facct"`). Only `beta`, `T`, and `phi = 11.75·T` differ. All are tracked against `estimator/lotrs_estimate.py`.

Probes the smaller-T regime against the same MSIS lattice as `BENCH_PARAMS`.

| Parameter | Value |
|-----------|-------|
| beta, N, T | 32, 32, 4 |
| phi | 47.0 (= 11.75 · T) |
| Signature size | ~22 KB |
| Expected attempts | ~12.3 |
| σ₀ | ≈ 2.1 × 10⁶ |

Not independently re-run by `estimator/lotrs_estimate.py` for this specific T — security argument is qualitative (smaller T against the same MSIS lattice). Use `BENCH_PARAMS` or `PRODUCTION_PARAMS` for any security-sensitive claim.

### BENCH_PARAMS (16-of-32 benchmark set)

Recommended benchmark point.

| Parameter | Value |
|-----------|-------|
| beta, N, T | 32, 32, 16 |
| phi | 188.0 (= 11.75 · T) |
| Signature size | ~23 KB |
| Expected attempts | ~12.3 |
| σ₀ | ≈ 8.5 × 10⁶ |
| Security (BKZ cost_pq) | 92 (binary ASIS) / 86 (DualMS ASIS) |

### PRODUCTION_PARAMS (Table 3 of the paper)

Full 50-of-100 set. Matches `estimator/LoTRS-Estimate-Output-N100T50.txt`. At d=128 the CRT-NTT path (`aux_ntt.py`) is active — signing and verification are tractable in pure Python.

| Parameter | Value |
|-----------|-------|
| beta, N, T | 100, 100, 50 |
| phi | 587.5 (= 11.75 · T) |
| Signature size | ~35 KB |
| Expected attempts | ~12.3 |
| σ₀ | ≈ 2.6 × 10⁷ |
| Security (BKZ cost_pq) | 92 (binary ASIS) / 86 (DualMS ASIS) |

## Benchmarks

Wall-clock timings of the high-level primitives are produced from the
Rust implementation, not from Python.  The authoritative numbers live
in [`../lotrs-rs/README.md`](../lotrs-rs/README.md) § *Benchmarks* —
see there for the full `N × T` grid at `N ∈ {32, 100}` and the
threshold sweep, plus the standalone RS (`T = 1`) and DualMS (`N = 1`)
edges, along with the averaging methodology (12 / 10 / 8 / 5
signing-seed samples per cell).  Empirical attempt-count means are in
the 25–40 range, materially above the estimator's `μ_total ≈ 12.3`
heuristic because the implementation also restarts on the additional
`w̃₀`-stability check used by the κ=1 commitment-compression path.

`examples/bench.rs` supports two modes:

```bash
# 1. Regression-comparison mode — named presets (TEST at d=32; and
# the d=128 presets BENCH_4OF32 / BENCH_PARAMS / PRODUCTION):
cargo run --release --example bench -- --with-prod

# 2. Grid mode — arbitrary (N, T) pairs on the shared d=128 lattice:
cargo run --release --example bench -- \
    --grid "32,1:4:8:16:32;100,1:5:10:25:50;1,2:4:8:16"
```

Full-grid wall-clock for the canonical reproduction is ~15 minutes on
a Ryzen AI 9 HX 370 (release build, single-threaded).

## Mapping to the paper

| Paper figure | Code method |
|-------------|-------------|
| Fig. 4, Setup | `LoTRS.setup()` |
| Fig. 4, KGen | `LoTRS.keygen()` |
| Fig. 4, KAgg | `LoTRS.kagg()` |
| Fig. 4, Sign_1 | `LoTRS.sign1()` |
| Fig. 4, SAgg | `LoTRS.sagg()` |
| Fig. 5, Sign_2 | `LoTRS.sign2()` — includes w̃₀ stability check |
| Fig. 5, Sign_bin | `LoTRS._sign_bin()` — w̃₀^(1) excluded from pi |
| Fig. 6, Vf | `LoTRS.verify()` — reconstructs ŵ₀^(1) from verification eq. |
| Fig. 1, Rej | `sample.rej()` |
| Fig. 1, RejOp | `sample.rej_op()` |
| Table 2 | `params.LoTRSParams` properties |
| Section 3.2, encoding | `codec.LoTRSCodec` |
| — (KAT infra) | `vectors.generate()` / `vectors.verify_vectors()` |
