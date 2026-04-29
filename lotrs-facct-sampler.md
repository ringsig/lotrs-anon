# LoTRS Large-Sigma FACCT-Style Sampler

This document specifies the large-sigma Gaussian sampler used by the
LoTRS development implementations for the masking distributions:

- `D_{sigma_0}`
- `D_{sigma_0_prime}`

It is separate from the paper draft because the paper currently does not
fix an implementation-level sampler for very large `sigma`. The goal
here is narrower:

- correct target distribution up to an explicitly stated truncation
- integer-only runtime arithmetic
- deterministic XOF-driven behavior suitable for Python/Rust interop
- practical implementation for `sigma` in the `10^7` range

This sampler is **not** used for the small binary-proof widths
`sigma_a`, `sigma_b`; those remain on the exact CDT path.

## 1. Target Distribution

For standard deviation `sigma > 0` and tail factor `tau = 14`, define

- `B = ceil(tau * sigma)`

The target law is the truncated centered discrete Gaussian on `Z`:

- `Pr[X = x]` is proportional to `exp(-x^2 / (2 sigma^2))`
- support is `x in {-B, ..., +B}`

Equivalently,

- `Pr[X = x] = rho_sigma(x) / sum_{y=-B}^{B} rho_sigma(y)`
- `rho_sigma(x) = exp(-x^2 / (2 sigma^2))`

This truncation matches the practical truncation already used by the
CDT sampler in the repo.

## 2. Randomness Source

All randomness is drawn from a caller-supplied SHAKE-256 XOF.

- integers are decoded little-endian
- rejection consumes fresh XOF output
- no ambient randomness source is used

This keeps the sampler deterministic from `(seed, tags, parameter set)`.

## 3. Proposal Distribution

Each coefficient is sampled independently.

1. Set `B = ceil(14 * sigma)`.
2. Sample `x` uniformly from `[-B, B]`.
3. Accept `x` with probability `exp(-x^2 / (2 sigma^2))`.
4. On rejection, repeat from step 2.

Because the proposal is uniform on the support, the accepted output has
exactly the truncated discrete-Gaussian law if the Bernoulli acceptance
step is exact.

## 4. FACCT-Style Bernoulli-Exp Decomposition

The runtime code avoids floating-point arithmetic. Acceptance is
performed with a FACCT-style decomposition of

- `u = x^2 / (2 sigma^2)`

Write

- `u = k * ln(2) + r`
- `k = floor(u / ln(2))`
- `0 <= r < ln(2)`

Then

- `exp(-u) = 2^{-k} * exp(-r)`

The Bernoulli acceptance is implemented as two independent tests:

1. `Ber(2^{-k})`
2. `Ber(exp(-r))`

The sample is accepted iff both tests accept.

## 5. Fixed-Point Format

The current Python reference uses:

- exponent fixed-point precision: `U_BITS = 48`
- Bernoulli threshold precision: `PROB_BITS = 64`
- polynomial degree: `20`

Constants:

- `LN2_Q = round(ln(2) * 2^U_BITS)`
- `PROB_SCALE = 2^PROB_BITS`

For a candidate `x`, the exponent is quantized as

- `u_q = floor((x^2 * 2^U_BITS) / (2 sigma^2))`

Then

- `k = floor(u_q / LN2_Q)`
- `r_q = u_q - k * LN2_Q`

Here `LN2_Q = round(ln(2) * 2^U_BITS)`, so `r_q` carries the induced
quantization error from this rounding. For the current parameters this
error is much smaller than the polynomial approximation error and is
treated as part of the fixed-point approximation budget.

## 6. Bernoulli(2^{-k})

To sample `Ber(2^{-k})`, read `k` fresh XOF bits and accept iff all are
zero.

- if `k >= 64`, process in 64-bit chunks
- the final partial chunk uses the low `k mod 64` bits

This is exact.

## 7. Bernoulli(exp(-r))

For `r in [0, ln(2))`, the implementation uses a degree-20 Taylor
polynomial for `exp(-r)`:

- `P(r) = sum_{i=0}^{20} (-1)^i * r^i / i!`

Coefficients are stored in `Q(PROB_BITS)` fixed point and evaluated with
Horner's rule. The current Python reference freezes these coefficients
as exact integers in `Q64` so ports do not depend on floating-point
rounding during setup. The resulting threshold is clamped into

- `[0, 2^PROB_BITS - 1]`

The Bernoulli draw is:

1. read a fresh 64-bit little-endian integer `c`
2. accept iff `c < floor(P(r) * 2^PROB_BITS)`

This is the only approximation step in the sampler. To keep Python and
Rust aligned, implementations should share the exact integer polynomial
coefficients, not regenerate them independently from floating-point code.

## 8. Relationship to FACCT

This is a **FACCT-style** sampler, not a line-by-line implementation of
the FACCT paper's binary base sampler / expander.

What is borrowed from FACCT is the implementation technique:

- decompose `exp(-u)` into a power-of-two part and a small remainder
- evaluate the remainder with fixed-point polynomial arithmetic
- avoid floating-point arithmetic at runtime

The proposal distribution here is simpler: uniform on the truncated
support.

## 9. Dispatch Rule in LoTRS

The repo uses two Gaussian backends:

- small `sigma`: exact CDT sampler
- large `sigma`: this FACCT-style sampler

The current threshold is based on CDT size:

- use FACCT-style when `ceil(14 * sigma) + 1 > 1_000_000`

So in practice:

- `sigma_a`, `sigma_b` use CDT
- `sigma_0`, `sigma_0_prime` use FACCT-style sampling for BENCH and
  PRODUCTION

## 10. Validation Requirements

Any implementation of this specification should satisfy:

1. deterministic output from the XOF stream
2. symmetry around zero
3. empirical variance close to `sigma^2`
4. empirical agreement with the exact CDT path at moderate `sigma`
   where both are feasible

The Python unit tests check all four.
