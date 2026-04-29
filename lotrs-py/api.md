# api.md — Internal API Reference

Python reference implementation of the LoTRS lattice-based structured threshold ring signature scheme. **Only kappa = 1 is implemented** (the concrete parameter set from the paper).

Module dependency order:

```
sample.py, aux_ntt.py  →  ring.py  →  params.py  →  lotrs.py  →  codec.py  →  vectors.py
```

### Paper alignment

The public API mirrors the protocol figures of the accompanying paper:

| Paper            | Python                              |
|------------------|-------------------------------------|
| Fig. 4 `Setup`   | `LoTRS.setup`                       |
| Fig. 4 `KGen`    | `LoTRS.keygen`                      |
| Fig. 4 `KAgg`    | `LoTRS.kagg`                        |
| Fig. 4 `Sign₁`   | `LoTRS.sign1`                       |
| Fig. 4 `SAgg`    | `LoTRS.sagg`                        |
| Fig. 5 `Sign₂`   | `LoTRS.sign2`                       |
| Fig. 5 `Sign_bin`| `LoTRS._sign_bin`                   |
| Fig. 6 `Vf`      | `LoTRS.verify`                      |
| Fig. 1 `Rej` / `RejOp` | `sample.rej` / `sample.rej_op`|
| Table 2 bounds   | properties on `LoTRSParams`         |
| Table 3 params   | `PRODUCTION_PARAMS`                 |

Implementation notes for this artifact:

- The concrete artifact implements the `kappa = 1` parameter sets used by the
  paper.
- The ring backend uses exact CRT-NTT multiplication where available, with
  schoolbook multiplication retained as a test oracle.
- Large masking widths use the FACCT-style sampler described in
  [`../lotrs-facct-sampler.md`](../lotrs-facct-sampler.md).
- The structured public-key table is hashed once to a 256-bit digest and that
  digest is used in the `H_agg`, `H_com`, and Fiat-Shamir inputs.

---

## Shared conventions

### Data representations

| Object | Python type | Notes |
|--------|------------|-------|
| ring element | `list[int]` | d coefficients in `[0, q)` (unsigned canonical) |
| signed ring element | `list[int]` | d coefficients in `[-q//2, q//2]` |
| vector of ring elements | `list[list[int]]` | n ring elements |
| matrix of ring elements | `list[list[list[int]]]` | row-major: M[row][col] is a ring element |
| PK table | `list[list[list[list[int]]]]` | pk_table[col][row] is a vector of k ring elements |
| secret key | `list[list[int]]` | l+k ring elements (short) |
| public key | `list[list[int]]` | k ring elements (t = A_bar * s) |
| signature | `dict` | keys: `pi`, `z_tilde`, `r_tilde`, `e_tilde` |
| proof transcript pi | `dict` | keys: `B_bin_hi`, `w_tilde_hi`, `x`, `f1`, `z_b` |

### Two rings

The scheme operates over two polynomial rings:

- **R_q** (`self.Rq`): main DualMS ring, modulus `q`. Used for signing keys, public keys, commitments, aggregated responses.
- **R_qhat** (`self.Rqh`): binary proof ring, modulus `q_hat`. Used for the commitment matrix G, the binary proof transcript (B_bin, A_bin, z_b), and the g-values.

The challenge `x` and the f-values are small enough (bounded by norm checks) to be valid in both rings. They are stored as R_q elements and embedded into R_qhat via `Rqh.from_centered(Rq.centered(...))` when needed.

### Ring arithmetic

All polynomial arithmetic is negacyclic in `R_q = Z_q[X] / (X^d + 1)`. Two multiplication paths:

- **Schoolbook** — O(d²) per multiply. Always available; used by tests as a correctness reference (`Ring._mul_schoolbook`).
- **CRT-NTT** (`aux_ntt.py`) — enabled at `Ring.__init__` when `d ∈ {128, 256}`. The scheme modulus q need not be NTT-friendly; fast multiplication is run over two auxiliary 48-bit primes `p1 = 2^48 − 16383`, `p2 = 2^48 − 19967`, each supporting radix-2 negacyclic NTTs at both supported dimensions. `P = p1·p2 ≈ 2^96` is large enough that CRT reconstruction recovers the *exact* integer product coefficient before final reduction mod q.

### Seed-derived randomness

Every function that needs randomness takes an explicit seed or XOF handle. No calls to `os.urandom` or `random` during scheme operations. The pattern is:

```python
xof = make_xof(seed, tag1, tag2, ...)    # SHAKE-256 with domain separation
poly = xof_sample_gaussian(xof, cdt, lam, d)
```

Tags are bytes, strings (ASCII-encoded), or ints (4-byte little-endian).

---

## `sample.py` — XOF-based polynomial samplers

Standalone deterministic samplers. Each takes an explicit SHAKE-256 XOF handle so the caller controls domain separation and seed derivation.

### CDT construction

```python
build_cdt(sigma: float, lam: int = 128, tailcut: int = 14) -> list[int]
```
Builds a cumulative distribution table for `|D_sigma|`. `cdt[k] = floor(Pr[|X| <= k] * 2^lam)` for `k = 0, ..., ceil(tailcut * sigma)`. Last entry clamped to `2^lam` for termination guarantee. Requires `mpmath` for high-precision CDF computation.

### XOF helpers

```python
make_xof(seed: bytes, *tags) -> SHAKE256
```
Creates a SHAKE-256 instance with domain-separated seed. Tags may be `bytes`, `str` (ASCII-encoded), or `int` (4-byte LE). Used everywhere for deterministic randomness derivation.

### Prepared Gaussian samplers

The backend is **declared on the parameter set** (`mask_sampler ∈ {"cdt", "facct"}`) and the `LoTRS` constructor calls the matching builder. There is no runtime `sigma`-threshold dispatch anywhere in the scheme.

```python
build_cdt_sampler(sigma, lam=128, tailcut=14) -> GaussianSampler  # kind == "cdt"
build_facct_sampler(sigma, lam=128, tailcut=14) -> GaussianSampler # kind == "facct"

needs_large_gaussian_sampler(sigma, tailcut=14, max_entries=1_000_000) -> bool
```

`needs_large_gaussian_sampler` is a **consistency helper** used only by `LoTRSParams.check()` — it returns `True` when `ceil(14*sigma) + 1 > 1_000_000`, i.e. when a CDT would need more than a million entries and therefore the parameter set must declare `mask_sampler="facct"`.  `check()` rejects a parameter set that declares `"cdt"` while its σ₀ would cross this threshold; `LoTRS.__init__` does not call the helper.

```python
@dataclass(frozen=True)
class LoTRSParams:
    ...
    mask_sampler: str = "cdt"        # "cdt" | "facct"

TEST_PARAMS.mask_sampler       == "cdt"
BENCH_4OF32.mask_sampler       == "facct"
BENCH_PARAMS.mask_sampler      == "facct"
PRODUCTION_PARAMS.mask_sampler == "facct"
```

`LoTRS.__init__` reads this field directly:

```python
if par.mask_sampler == "cdt":
    self.gauss_0 = build_cdt_sampler(par.sigma_0, par.lam)
    self.gauss_0p = build_cdt_sampler(par.sigma_0_prime, par.lam)
elif par.mask_sampler == "facct":
    self.gauss_0 = build_facct_sampler(par.sigma_0, par.lam)
    self.gauss_0p = build_facct_sampler(par.sigma_0_prime, par.lam)
else:
    raise ValueError(f"unsupported mask_sampler {par.mask_sampler!r}")
```

Sampling calls are split matching the two builders:

```python
xof_sample_gaussian(xof, cdt, lam, d)         # consumes a CDT list
xof_sample_gaussian_facct(xof, prepared, d)   # consumes FACCT params
```

Call sites dispatch on `sampler.kind` explicitly:

```python
if self.gauss_0.kind == "cdt":
    y = xof_sample_gaussian(xof_y, self.gauss_0.data, par.lam, d)
else:                             # "facct"
    y = xof_sample_gaussian_facct(xof_y, self.gauss_0.data, d)
```

### Polynomial samplers

```python
xof_sample_uniform(xof, q: int, d: int) -> list[int]
```
Uniform polynomial in `[0, q)^d` via rejection sampling. Uses `ceil(log2(q))`-bit chunks.

```python
xof_sample_gaussian(xof, cdt: list[int], lam: int, d: int) -> list[int]
```
Discrete Gaussian `D_sigma` via CDT lookup. For each coefficient: reads `lam` bits of magnitude randomness, binary-searches `cdt` for the magnitude, then reads a separate sign byte (LSB only). The CDT is scaled to `2^lam` so the full `lam`-bit value gives an unbiased comparison. Returns signed integers.

```python
xof_sample_gaussian_facct(xof, prepared: dict, d: int) -> list[int]
```

Large-sigma FACCT-style sampler. The target law is the truncated
discrete Gaussian on `[-ceil(14*sigma), +ceil(14*sigma)]`. Runtime uses
only integer arithmetic:

- signed-uniform proposal on the truncated support
- exponent decomposition `exp(-u) = 2^{-k} * exp(-r)`
- exact Bernoulli `Ber(2^{-k})` from XOF bits
- fixed-point polynomial approximation of `exp(-r)` on `[0, ln 2)`

See the dedicated sampler specification in
[`../lotrs-facct-sampler.md`](../lotrs-facct-sampler.md).

```python
xof_sample_short(xof, eta: int, d: int) -> list[int]
```
Uniform polynomial with coefficients in `[-eta, eta]`. Rejection sampling with `2*eta+1` range.

```python
xof_sample_challenge(xof, w: int, d: int) -> list[int]
```
Challenge element `x` in `C = { x in R : ||x||_inf = 1, ||x||_1 = w }`. Picks `w` distinct positions by rejection sampling over the smallest power of 2 that contains `d` (unbiased for non-power-of-2 `d`; zero reject rate when `d` is a power of 2), then assigns each `+/-1` from a packed sign word. Returns signed integers. Callers normalize to canonical `[0, q)` form via `Rq.from_centered()`.

```python
xof_sample_ternary(xof, d: int) -> list[int]
```
Uniform ternary polynomial (coefficients in `{-1, 0, 1}`). Reads 2 bits per coefficient with rejection of the value 3.

```python
xof_sample_bounded(xof, bound: int, d: int) -> list[int]
```
Alias for `xof_sample_short(xof, bound, d)`.

### Rejection sampling (Fig. 1)

```python
rej(xof, z_flat: list[int], v_flat: list[int], phi: float, K: float) -> bool
```
Standard rejection sampling `Rej(z, v, phi, K)` from Fig. 1 (left column). Uses `mu(phi) = exp(12/phi + 1/(2*phi^2))`. Returns `True` on accept, `False` on reject.

```python
rej_op(xof, z_flat: list[int], v_flat: list[int], phi: float, K: float) -> bool
```
Optimised rejection sampling `RejOp(z, c, phi, K)` from Fig. 1 (right column). Uses `mu(phi) = exp(1/(2*phi^2))` (no 12/phi term). Returns `True` on accept, `False` on reject.

### Flat helpers

```python
_flat(vec_of_polys: list[list[int]]) -> list[int]
```
Flatten a vector of ring elements into a single coefficient list for rejection sampling.

---

## `aux_ntt.py` — CRT negacyclic multiplication via two 48-bit primes

Because `q ≡ 5 mod 8` blocks a native negacyclic NTT in `R_q`, we do fast multiplication over two auxiliary primes and reconstruct the exact integer product coefficient via CRT. The scheme's moduli are left untouched.

```python
AUX_P1 = 281474976694273           # 2^48 − 16383
AUX_P2 = 281474976690689           # 2^48 − 19967
AUX_P  = AUX_P1 * AUX_P2           # ≈ 2^96
SUPPORTED_D = (128, 256)
AUX_MERSENNE_C1 = 16383            # for C / asm Mersenne-style reduction
AUX_MERSENNE_C2 = 19967            # (unused in the Python reference)
```

Both auxiliary primes are `≡ 1 mod 512`, so each supports a radix-2 negacyclic NTT at `d ∈ {128, 256}`.

### Primitive root

```python
find_primitive_2d_root(p: int, d: int) -> int
```
Smallest `x ∈ [2, p)` with `x^d ≡ −1 (mod p)` (multiplicative order exactly `2d`). Requires `2d | p − 1`; raises `ValueError` otherwise.

### Twiddle table (bit-reversed, FIPS 204 convention)

```python
make_zeta_table(p, d, psi) -> list[int]      # zetas[k] = psi^bitrev(k, log2 d) mod p
```

### Radix-2 negacyclic transforms

```python
ntt(f, zetas, p)   -> list[int]      # natural in, bit-reversed out
intt(f, zetas, p)  -> list[int]      # bit-reversed in, natural out (includes 1/d scale)
```

### `CRTBackend(q: int, d: int)` class

The public entry point used by `Ring`. Holds one `_AuxContext` per auxiliary prime at this `d`, and the precomputed inverse `p1^−1 mod p2` used for Garner CRT.

```python
backend.mul(a, b) -> list[int]
```
Full negacyclic multiplication:
1. Centre-lift each coefficient of `a`, `b` from `[0, q)` to `(−q/2, q/2]` and reduce mod `p1` and `p2`.
2. Forward NTT both operands under each auxiliary prime.
3. Pointwise multiply.
4. Inverse NTT.
5. Garner-CRT to the centred integer in `(−P/2, P/2]`.
6. Reduce mod `q`.

Correctness rests on `P = p1·p2 > 2·max|c_k|`. For the production moduli this holds with ≥14 bits of headroom even at `d=256` with inner-product accumulation.

### Mersenne-style reduction note

Because each auxiliary prime has the form `2^48 − c` with `c ∈ {16383, 19967}`, a C / asm port can skip general division:

    x mod p  ≡  (hi · c + lo)  mod p       for  x = hi · 2^48 + lo

where `hi, lo ∈ [0, 2^48)`. The multiply constant is 14–15 bits, so the whole reduction collapses to two 64-bit fused multiply-adds on modern hardware. In pure Python the savings are negligible and we use direct `% p`.

---

## `ring.py` — Polynomial ring `R_q = Z_q[X] / (X^d + 1)`

### `Ring(q: int, d: int)` class

**Attributes:** `q`, `d`, `half_q = q // 2`, `_ntt_ctx` (a `CRTBackend` when the NTT fast path applies, else `None`).

#### Element creation

| Method | Returns | Description |
|--------|---------|-------------|
| `zero()` | `list[int]` | Additive identity |
| `one()` | `list[int]` | Multiplicative identity (constant 1) |
| `const(c)` | `list[int]` | Constant polynomial `c mod q` |

#### Coefficient conversions

| Method | Description |
|--------|-------------|
| `reduce(a)` | All coefficients into `[0, q)` |
| `centered(a)` | Unsigned `[0,q)` to signed `[-q//2, q//2]` |
| `from_centered(a)` | Signed to unsigned `[0, q)` |

#### Element arithmetic

| Method | Description |
|--------|-------------|
| `add(a, b)` | `a + b mod q` |
| `sub(a, b)` | `a - b mod q` |
| `neg(a)` | `-a mod q` |
| `scale(c, a)` | Integer scalar `c` times polynomial `a` |
| `mul(a, b)` | Negacyclic convolution `a * b` in `R_q`. Uses NTT path when `_ntt_ctx` is set, else schoolbook. |
| `_mul_schoolbook(a, b)` | Reference O(d²) schoolbook multiplication (used by tests and as the fallback path). |

#### Norms (on unsigned-canonical input)

| Method | Description |
|--------|-------------|
| `inf_norm(a)` | `max |centered(a_i)|` |
| `l2_norm_sq(a)` | Sum of squares of centered coefficients |
| `l1_norm(a)` | Sum of absolute centered coefficients |

#### Vector operations

Vectors are `list[list[int]]` — a list of ring elements.

| Method | Description |
|--------|-------------|
| `vec_zero(n)` | Zero vector of length n |
| `vec_add(u, v)` | Component-wise add |
| `vec_sub(u, v)` | Component-wise subtract |
| `vec_neg(v)` | Component-wise negate |
| `vec_scale(c_poly, v)` | Ring element `c` times each component |
| `vec_scale_int(c, v)` | Integer `c` times each component |
| `inner(u, v)` | Inner product `<u, v> = sum u_i * v_i` in `R_q` |
| `vec_inf_norm(v)` | Max `inf_norm` across components |
| `vec_l2_norm_sq(v)` | Sum of `l2_norm_sq` across components |
| `vec_concat(*vecs)` | Concatenate multiple vectors |

#### Matrix operations

Matrices are `list[list[list[int]]]` — row-major.

```python
ring.mat_vec(M, v) -> list[list[int]]    # M * v
```

#### Centered decomposition

```python
ring.centered_decompose(a, K: int) -> tuple[list[int], list[int]]
```
For each coefficient `c` (in centered form), writes `c = high * 2^K + low` with `|low| <= 2^{K-1}`. Returns `(high_poly, low_poly)` both in unsigned canonical form. Used for Bai-Galbraith compression of `B_bin`, `A_bin`, and `w_tilde`.

---

## `params.py` — Parameter sets

### `LoTRSParams` dataclass (frozen)

#### Base parameters

| Field | Type | Paper notation | Description |
|-------|------|---------------|-------------|
| `name` | str | — | Identifier |
| `d` | int | d | Ring dimension (power of 2) |
| `q` | int | q | Main modulus (prime, 5 mod 8) |
| `q_hat` | int | q-hat | Binary proof modulus (prime, 5 mod 8) |
| `kappa` | int | kappa | `log_beta(N)` — must be 1 |
| `beta` | int | beta | Base for column index encoding |
| `T` | int | T | Threshold |
| `k` | int | k | Rows of A |
| `l` | int | l | Columns of A |
| `l_prime` | int | l' | Columns of B |
| `n_hat` | int | n-hat | Rows of G |
| `k_hat` | int | k-hat | Extra columns in G |
| `w` | int | w | Challenge weight |
| `eta` | int | eta | Secret key bound |
| `phi` | float | phi_0 | Rejection sampling slack for z_u |
| `phi_a` | float | phi_a | Rejection sampling slack for f_1 |
| `phi_b` | float | phi_b | Rejection sampling slack for z_b |
| `K_A` | int | K_A | Bit-dropping for A_bin |
| `K_B` | int | K_B | Bit-dropping for B_bin |
| `K_w` | int | K_w | Bit-dropping for w_tilde (K_{w,0}) |
| `lam` | int | lambda | Security parameter (default 128) |
| `max_attempts` | int | — | Max signing restarts (default 2000) |
| `eta_prime` | int | eta' | Dual secret bound (default: use eta) |
| `tail_t` | float | t | Gaussian tail factor for aggregated l2 bounds (default 1.2) |
| `eps_tot` | float | ε_tot | Tail budget for binary-proof bound checks (default 0.01, matching `estimator/lotrs_finder.py`); each of the four checks sized for violation probability ≤ ε_tot / 4 |

#### Derived properties

| Property | Formula (Table 2 / `estimator/lotrs_finder.py`) | Description |
|----------|-------------------------------------------------|-------------|
| `N` | `beta^kappa` | Ring size |
| `eta_p` | `eta_prime` if >=0 else `eta` | Dual secret bound |
| `B_a` | `sqrt(kappa * w)` | Base mask scale |
| `B_b` | `sqrt(d * (n_hat + k_hat)) * w` | Binary proof z_b scale |
| `B_0` | `sqrt(d(l+k)) * eta * w^{kappa+1}` | z_u shift bound (kappa=1) |
| `B_0_prime` | `sqrt(d(l'+k)) * eta_p * w^{kappa+1}` | r_u shift bound |
| `sigma_0` | `phi * B_0` | Gaussian width for y_{u,0}, z_u |
| `sigma_0_prime` | `phi * B_0_prime` | Gaussian width for r_{u,0} |
| `sigma_a` | `phi_a * B_a` | Gaussian width for f_1 (binary proof, R_qhat) |
| `sigma_b` | `phi_b * B_b` | Gaussian width for z_b (binary proof, R_qhat) |
| `sigma_s` | `(2d/sqrt(2pi)) * q^{k/(l+k) + 2/(d(l+k))}` | DualMS regularity threshold |
| `sigma_s_prime` | Same with l' | Dual regularity threshold |
| `mu_phi` | `exp(12/phi + 1/(2*phi^2))` | Expected restarts per signer |
| `B_f1` | `1 + tau_{f_1} * sigma_a`, `M_{f_1} = kappa(beta-1)d` | Per-coefficient bound on transmitted f_{j,i} (paper Appendix) |
| `B_f0` | `1 + tau_{f_0} * sqrt(beta-1) * sigma_a`, `M_{f_0} = kappa*d` | Per-coefficient bound on reconstructed f_{j,0} |
| `B_g0` | `w * B_f0 + d * B_f0^2` | g_0 = f_{j,0}(x-f_{j,0}) bound |
| `B_g1` | `w * B_f1 + d * B_f1^2` | g_{j,i≥1} = f_{j,i}(x-f_{j,i}) bound |
| `B_eta_w` | `2^{K_w-1} * kappa * w^{kappa-1}` | Residual bound (analysis only) |
| `G_cols` | `n_hat + k_hat + 2*kappa*beta` | Columns of binary proof matrix G |

```python
par.check()            # basic consistency: primes mod 8, d power of 2,
                       # sigma ranges against the right modulus (q vs q_hat)
par.check_security()   # adds: range-proof condition, challenge-difference
                       # invertibility (Lemma 1).  Toy params may pass
                       # check() but not check_security().
```

#### Concrete parameter sets

All three d=128 sets share the same lattice: `k=12, l=5, l'=6, n̂=10, k̂=8`, `q = 274877906837` (largest prime ≤ 2^38 with `q ≡ 5 mod 8`), `q_hat = 8589934237` (largest prime ≤ 2^33), `K_A=20, K_B=5, K_w=5`, `phi_a=50, phi_b=4`.  Only `T`, `phi`, and `beta` vary.  All are tracked against `estimator/lotrs_estimate.py`.

- `TEST_PARAMS` — small test set (d=32, N=4, T=2, q=4194389, q_hat=7000061). Fast, not secure. Uses distinct q/q_hat to exercise both rings. NTT path is disabled at d=32; schoolbook is used. `mask_sampler = "cdt"`.
- `BENCH_4OF32` — 4-of-32 benchmark variant (d=128, β=32, T=4, `phi = 47 = 11.75·T`). `mask_sampler = "facct"` (σ₀ ≈ 2.1 × 10⁶).  Probes the small-T regime at N=32.
- `BENCH_PARAMS` — 16-of-32 benchmark set (d=128, β=32, T=16, `phi = 188 = 11.75·T`).  Signature ~23 KiB, ~12.3 expected attempts.
- `PRODUCTION_PARAMS` — 50-of-100 set (d=128, β=100, T=50, `phi = 587.5 = 11.75·T`).  Signature ~35 KiB, ~12.3 expected attempts. Matches `estimator/LoTRS-Estimate-Output-N100T50.txt`.

---

## `lotrs.py` — LoTRS scheme

### `LoTRS(par: LoTRSParams)` class

Constructor builds the Gaussian samplers declared on the parameter set — always CDT for `sigma_a`, `sigma_b`, and CDT or FACCT for `sigma_0`, `sigma_0_prime` according to `par.mask_sampler`.  At TEST parameters (mask_sampler="cdt") the CDT build is a one-time cost of a few seconds; at `BENCH_*` / `PRODUCTION` (mask_sampler="facct") the mask samplers are near-instantaneous.  Raises `NotImplementedError` if `par.kappa != 1`.

**Attributes:** `par`, `Rq` (Ring at q), `Rqh` (Ring at q_hat), `gauss_0`, `gauss_0p`, `gauss_a`, `gauss_b` — the four pre-resolved samplers.

### Public API — Fig. 4

```python
scheme.setup(seed: bytes) -> bytes
```
`Setup(1^lambda)`. Returns `pp = seed` (32 bytes). Public matrices A and G are expanded from this seed on demand.

```python
scheme.keygen(pp: bytes, seed: bytes) -> tuple[list, list]
```
`KGen(pp)`. Returns `(sk, pk)` where `sk` is `l+k` ring elements and `pk` is `k` ring elements satisfying `pk = [A | I] * sk`.

```python
scheme.kagg(pk_table) -> list
```
`KAgg(PK)`. Returns `N` aggregated column keys. The implementation first computes a 256-bit SHAKE256 digest of the canonical PK-table serialization with domain tag `pk`, then uses `alpha_u = H_agg(pk_hash, u)` so the full PK table is hashed once rather than once per row.

```python
scheme.sign1(pp, sk_u, row_u, ell, mu, pk_table, rho, attempt)
    -> tuple[dict, list]
```
`Sign_1` for one signer. Returns `(state, commitments)` where `commitments[j]` is a vector of `k` ring elements (`w_{u,j}`). The `state` dict carries all values needed by `sign2`.

```python
scheme.sagg(sigmas: list[dict]) -> dict
```
`SAgg`. Aggregates `T` partial signatures into a final signature. Checks consistency of all pi fields across signers (x, B_bin_hi, w_tilde_hi, f1, z_b). Splits each `z_u` into `(z', z'')` and each `r_u` into `(r', r'')`, then sums component-wise.

### Public API — Fig. 5

```python
scheme.sign2(state: dict, all_coms: list, pk_table) -> dict | None
```
`Sign_2` for one signer. Takes the signer's state from `sign1` and all signers' round-1 commitments. Returns `{pi, z_u, r_u}` on success, `None` on rejection (triggers restart). Internally:
1. Aggregates w_tilde and decomposes (high/low bits)
2. Calls `_sign_bin` for the binary selection proof
3. Checks w̃₀ stability: rejects if `||w̃₀^(0)||_inf > 2^{K_w-1}` (kappa=1)
4. Computes z_u response with rejection sampling
5. Computes auxiliary r_u response

### Public API — Fig. 6

```python
scheme.verify(pp, mu, sigma: dict, pk_table) -> bool
```
`Vf`. Full verification of an aggregated signature. Checks:
1. Norm bounds on `f1` (`<= B_f1`) and `z_b` (`<= 6 phi_b B_b`), centered infinity norms
2. L2 norm bounds on `z_tilde`, `r_tilde`, `e_tilde`
3. Reconstructs `f_{j,0}` and checks `||f_0||_inf <= B_f0`
4. Quadratic terms `g0`, `g1` within bounds (computed in R_qhat)
5. Bai-Galbraith low-bit check on reconstructed `A_hat_bin`
6. Reconstructs `w_hat_0` from verification equation, decomposes to get `w_hat_0^(1)`
7. Fiat-Shamir hash check: `x == H(mu, A_hat^(1), B_bin^(1), w_hat_0^(1), pk_hash)`

The w̃₀^(1) is NOT in the signature — the verifier reconstructs it.

### Convenience

```python
scheme.sign(pp, sks, ell, mu, pk_table, signing_seed) -> dict
```
Runs the full two-round ceremony with automatic restart loop. For attempt `i`, derives `rho = SHAKE256(signing_seed || "rho" || i)`. Raises `RuntimeError` after `max_attempts`.

### Internal methods

#### Matrix expansion

| Method | Description |
|--------|-------------|
| `_expand_A(pp)` | `A in R_q^{k x l}` from `SHAKE256(pp \|\| "A" \|\| i \|\| j)` |
| `_expand_G(pp)` | `G in R_{q_hat}^{n_hat x G_cols}` from `SHAKE256(pp \|\| "G" \|\| i \|\| j)` |
| `_augment_I(M, k)` | `[M \| I_k]` — append identity block |

#### Hash functions

| Method | Description |
|--------|-------------|
| `_pk_hash(pk_table)` | 256-bit SHAKE256 digest of the canonical PK-table serialization, domain-separated with tag `pk`. |
| `_hash_agg(pk_hash, u)` | `H_agg(H(PK), u)` — per-row aggregation coefficient in C. Returns canonical `[0, q)` form. |
| `_hash_com(pk_hash, mu)` | `H_com(H(PK), mu)` — commitment matrix B in R_q^{k x l'}. |
| `_hash_fs(mu, A_hi, B_hi, w_hi, pk_hash)` | Fiat-Shamir hashing with the 32-byte PK digest. Uses **8-byte signed LE** per coefficient for decomposed values and appends the digest. Returns canonical form. |

`H_agg`, `H_com`, and Fiat-Shamir all consume the 32-byte `pk_hash` rather than the full PK table. The signing and verification paths compute the digest once and reuse it. This avoids repeated serialization and hashing of the full structured public-key table, which is multi-MiB for benchmark and production parameter sets.

#### Selector polynomial

| Method | Description |
|--------|-------------|
| `_expand_a_coeffs(rho, attempt)` | Gaussian mask coefficients `a_{j,i}` from rho. Returns signed integer lists (ring-independent). |
| `_compute_p_coeffs(a_coeffs, ell)` | Selector polynomial coefficients `p_{i,j}` in R_q (Remark 2) |

#### Binary selection proof (Fig. 5, right column)

```python
scheme._sign_bin(pp, ell, mu, w_tilde_hi, pk_table, rho, attempt) -> dict | None
```
`Sign_bin`. Produces the binary proof transcript `pi` or `None` on restart. Steps:
1. Expand G (R_qhat), decompose ell into base-beta digits, form one-hot b
2. Compute a, c, d vectors in R_qhat; sample r_b (ternary) and r_a (Gaussian)
3. Commit: B_bin = G(r_b, b, c)^T, A_bin = G(r_a, a, d)^T — both in R_qhat
4. Bai-Galbraith decomposition of B_bin and A_bin
5. Fiat-Shamir challenge x (uses all w_tilde_hi including j=0 internally)
6. Response z_b = r_a + x*r_b with RejOp (R_qhat)
7. Masked opening f_{j,i} = x*delta + a_{j,i} as signed lists; Rej on f1
8. Quadratic g_{j,i} = f * (x - f) computed in Z, reduced into R_qhat; norm checks
9. A_hat_bin low-bit stability check (R_qhat)
10. **Return pi excluding w_tilde_hi[0]** — verifier reconstructs it

#### Verification helpers

| Method | Description |
|--------|-------------|
| `_reconstruct_f(x, f1_flat)` | Recover full f_{j,i} from f_1 and challenge (Vf line 8). Returns R_q elements. |
| `_compute_g(x, f_full)` | g_{j,i} = f * (x - f) in R_q. Used only internally. |
| `_compute_g_qhat(x, f_full)` | g_{j,i} = f * (x - f) computed in Z, reduced into R_qhat. Used by verify for binary proof. |
| `_interleave_g(g0, g1)` | Flatten g0/g1 into G-multiplication order |
| `_flatten_f(f_full)` | Flatten f_full[j][i] into a single list |
| `_base_digits(n, base, num_digits)` | Decompose n in base, LSB first |

### Module-level helpers

```python
_poly_mul_signed(a, b, d) -> list[int]
```
Negacyclic convolution of two signed coefficient lists over Z (not mod q). Used for selector polynomial computation and quadratic terms where the result is later reduced into the appropriate ring.

---

## `codec.py` — Compact serialization

### Bitstream helpers

```python
BitWriter()
  .write_bits(value, n)         # n bits, LSB-first
  .write_unary(value)           # value ones + stop zero
  .pad_to_byte()
  .to_bytes() -> bytes

BitReader(data: bytes)
  .read_bits(n) -> int
  .read_unary(max_val) -> int   # DoS guard: rejects if > max_val
  .check_padding()              # rejects nonzero padding bits
  .check_exhausted()            # rejects trailing bytes
```

### Packing primitives

```python
_pack_fixed(coeffs, dx, offset=0) -> bytes
_unpack_fixed(data, d, dx, offset=0) -> (list[int], bytes_consumed)
```
Fixed-width packing at `dx` bits per coefficient. Signed values use `offset` shift (caller sets offset = bound). Rejects overflow on encode, nonzero padding on decode.

```python
_pack_rice(coeffs, rice_k, bound) -> bytes
_unpack_rice(data, d, rice_k, bound) -> (list[int], bytes_consumed)
```
Golomb-Rice encoding for signed integer coefficients. Format per coefficient: `low` (rice_k bits) + unary(`high`) + sign (if nonzero). Byte-aligned per polynomial. Rejects `|coeff| > bound` and caps unary runs.

```python
optimal_rice_k(sigma) -> int
```
Optimal Rice parameter: `floor(log2(1.1774 * sigma))`.

```python
_pack_challenge(poly, w, d) -> bytes
_unpack_challenge(data, w, d) -> (list[int], bytes_consumed)
```
Challenge encoding: `w` sorted positions (ceil(log2(d)) bits each) + `w` sign bits. Rejects duplicate, non-ascending, or out-of-range positions.

### `LoTRSCodec(par: LoTRSParams)` class

Constructor derives all encoding parameters (bit widths, Rice k, bounds) from the scheme parameters.

**Encoding attributes:**

| Attribute | Value | Used for |
|-----------|-------|----------|
| `dx_pk` | ceil(log2(q)) | Public key coefficients |
| `dx_bbin` | ceil(log2(q_hat)) - K_B | B_bin high bits |
| `dx_whi` | ceil(log2(q)) - K_w | w_tilde high bits (unused for kappa=1) |
| `rice_f1` | optimal_rice_k(sigma_a) | f1 (masked binary opening) |
| `rice_zb` | optimal_rice_k(sigma_b) | z_b (binary proof response) |
| `rice_zt` | optimal_rice_k(sigma_z) | z_tilde (aggregated response) |
| `rice_rt` | optimal_rice_k(sigma_r) | r_tilde (aggregated dual response) |
| `rice_et` | optimal_rice_k(sigma_e) | e_tilde (aggregated error) |

**Public API:**

```python
codec.pp_encode(pp) -> bytes              # 32 bytes
codec.pp_decode(data) -> bytes
codec.sk_encode(seed) -> bytes            # 32 bytes
codec.sk_decode(data) -> bytes
codec.pk_encode(pk) -> bytes              # k polys, fixed-width
codec.pk_decode(data) -> list[list[int]]
codec.sig_encode(sigma) -> bytes          # full aggregated signature
codec.sig_decode(data) -> dict            # raises ValueError on malformation
codec.sizes() -> dict[str, int]           # estimated byte sizes per component
codec.total_size_estimate() -> int
codec.print_sizes() -> int                # prints breakdown, returns total
```

**Signature wire format** (all components byte-aligned):

| Field | Count (kappa=1) | Encoding | Bits per coeff |
|-------|-----------------|----------|---------------|
| B_bin^(1) | n_hat polys | fixed | log(q_hat) - K_B |
| w_tilde^(1) | **0 polys** (j >= 1 only) | fixed | log(q) - K_w |
| x | 1 challenge | challenge | w * (log(d) + 1) |
| f1 | kappa*(beta-1) polys | Rice | ~log(4.13 * sigma_a) |
| z_b | (n_hat+k_hat) polys | Rice | ~log(4.13 * sigma_b) |
| z_tilde | l polys | Rice | ~log(4.13 * sigma_z) |
| r_tilde | l' polys | Rice | ~log(4.13 * sigma_r) |
| e_tilde | k polys | Rice | ~log(4.13 * sigma_e) |

For kappa = 1, w_tilde^(1) is omitted entirely (0 bytes). The verifier reconstructs w_hat_0^(1) from the verification equation.

**Deserialization validation (all decoders):**
- Exact length match (no trailing bytes)
- All coefficients within declared range
- Nonzero padding bits rejected
- Rice unary runs capped at `(bound >> rice_k) + 1`
- Challenge positions strictly ascending and < d

---

## `vectors.py` — Deterministic test-vector generation

```python
generate(par=None) -> dict                # JSON-serialisable test vector
verify_vectors(vec: dict) -> (bool, list[str])   # re-derive and compare
```

**CLI usage:**

```bash
python vectors.py --out vectors.json          # generate
python vectors.py --verify vectors.json       # verify
```

**Fixed seeds:**
- pp: `00 01 02 ... 1f`
- signer (col, row): `bytes([col ^ 0x40, row ^ 0x80]) + b'\x00' * 30`
- signing: `aa aa ... aa`

**JSON schema (v1):** params, seeds, keygen records (per-signer pk hex), signature (hex bytes, length, challenge, norms), verification result, encoding parameters, size estimates.

`verify_vectors` re-derives every key and the signature from seeds, then checks byte-for-byte equality of all encoded objects.

---

## Dependency graph

```
lotrs.py
  ├── ring.py    (Ring)
  ├── sample.py  (all samplers, rej, rej_op, make_xof)
  └── params.py  (LoTRSParams)
codec.py
  ├── ring.py    (Ring)
  └── params.py  (LoTRSParams)
vectors.py
  ├── lotrs.py   (LoTRS)
  ├── codec.py   (LoTRSCodec)
  └── params.py  (TEST_PARAMS)
test_codec.py
  └── codec.py, lotrs.py, vectors.py, params.py
```

---

## Signature data layout

### Aggregated signature (`dict`)

| Key | Type | Dimensions | Description |
|-----|------|------------|-------------|
| `pi` | dict | — | Binary proof transcript |
| `z_tilde` | list[poly] | l elements | Aggregated signing response (R_q) |
| `r_tilde` | list[poly] | l' elements | Aggregated dual response (R_q) |
| `e_tilde` | list[poly] | k elements | Aggregated error term (R_q) |

### Proof transcript `pi` (`dict`)

| Key | Type | Dimensions | Description |
|-----|------|------------|-------------|
| `B_bin_hi` | list[poly] | n_hat elements | High bits of binary commitment (R_qhat) |
| `w_tilde_hi` | list[list[poly]] | (kappa-1) x k elements | High bits of w_tilde for j >= 1. **Empty for kappa=1.** |
| `x` | poly | 1 element | Fiat-Shamir challenge in C (R_q, canonical form) |
| `f1` | list[poly] | kappa*(beta-1) elements | Masked binary opening, i >= 1 entries (R_q) |
| `z_b` | list[poly] | n_hat+k_hat elements | Binary proof response (R_qhat) |
