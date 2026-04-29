LoTRS — primitive benchmarks (release build)
Timings averaged over multiple signing seeds to smooth out
the high-variance rejection-sampling attempt count.

All rows use the shared d=128 lattice:
  d=128, κ=1, k=12, l=5, l'=6, n̂=10, k̂=8, w=31, η=1
  q = 274877906837 (largest prime ≤ 2^38 with q ≡ 5 mod 8)
  q_hat = 8589934237 (largest prime ≤ 2^33 with q_hat ≡ 5 mod 8)
  phi_a = 50, phi_b = 4, phi = 11.75·T
  K_A = 20, K_B = 5, K_w = 5, eps_tot = 0.01
  mask_sampler = facct, tail_t = 1.2

### Primitive timings

Sign / Verify are arithmetic means over the listed
number of signing seeds; KeyGen / KAgg are single-run
(deterministic given pp / PK).

| parameter set | d | N | T | samples | KeyGen | KAgg | Sign | Verify |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `N=32,T=1` | 128 | 32 | 1 | 12 | 2.51 ms | 8.73 ms | 318.6 ms | 31.39 ms |
| `N=32,T=4` | 128 | 32 | 4 | 12 | 1.63 ms | 34.83 ms | 1375.0 ms | 59.84 ms |
| `N=32,T=8` | 128 | 32 | 8 | 10 | 1.62 ms | 71.14 ms | 4286.3 ms | 99.64 ms |
| `N=32,T=16` | 128 | 32 | 16 | 8 | 1.61 ms | 144.1 ms | 6008.2 ms | 177.4 ms |
| `N=32,T=32` | 128 | 32 | 32 | 5 | 1.60 ms | 287.8 ms | 6052.2 ms | 334.5 ms |
| `N=100,T=1` | 128 | 100 | 1 | 12 | 1.63 ms | 26.89 ms | 1385.5 ms | 90.48 ms |
| `N=100,T=5` | 128 | 100 | 5 | 12 | 1.61 ms | 134.4 ms | 4600.2 ms | 208.2 ms |
| `N=100,T=10` | 128 | 100 | 10 | 10 | 1.61 ms | 276.0 ms | 12.40 s | 364.3 ms |
| `N=100,T=25` | 128 | 100 | 25 | 8 | 1.63 ms | 701.6 ms | 29.24 s | 827.5 ms |
| `N=100,T=50` | 128 | 100 | 50 | 5 | 1.63 ms | 1409.7 ms | 50.39 s | 1609.1 ms |
| `N=1,T=2` | 128 | 1 | 2 | 12 | 1.71 ms | 0.627 ms | 160.5 ms | 5.77 ms |
| `N=1,T=4` | 128 | 1 | 4 | 12 | 1.62 ms | 1.11 ms | 364.5 ms | 6.37 ms |
| `N=1,T=8` | 128 | 1 | 8 | 10 | 1.61 ms | 2.22 ms | 889.6 ms | 7.65 ms |
| `N=1,T=16` | 128 | 1 | 16 | 8 | 1.62 ms | 4.44 ms | 1395.1 ms | 10.15 ms |

### Sign / Verify breakdown — DualMS multi-sig vs RS binary proof

Sign_DualMS = sign1 (T-fold commitments) + sign2 minus the
binary proof + sagg.  Sign_RS = `sign_bin`, the binary ring
proof, summed over the T independent signers and over every
rejection-sampling attempt.  Both are wall-clock totals of
the sequential reference implementation, so
Sign_DualMS + Sign_RS ≈ Sign minus the one-shot context
expansion (expand_A/G/B, NTT prep, KAgg α_u).  `attempts` is
the mean attempt count per accepted signature.

Note: `sign_bin` runs once per signer per attempt in this
implementation (`sagg` asserts the T proofs match).  An
optimised multi-signer protocol that broadcasts a single
`pi` would have standalone RS cost ≈ Sign_RS / T.

Verify_DualMS = z̃/r̃/ẽ bounds + `A·z̃ + B·r̃ + ẽ` reconstruction
+ KAgg + ring-keyed `pk_sum` + w̃₀ recovery + closing FS hash.
Verify_RS = f1/f0/g0/g1 bounds + A_hat_bin reconstruction &
low-bit check.  pk_sum is attributed to DualMS because it's
part of the LHS=RHS multi-sig closing check.

| parameter set | attempts | Sign_DualMS | Sign_RS | Sign | Verify_DualMS | Verify_RS | Verify |
|---|---:|---:|---:|---:|---:|---:|---:|
| `N=32,T=1` | 17.67 | 213.5 ms | 93.47 ms | 318.6 ms | 19.35 ms | 10.93 ms | 31.39 ms |
| `N=32,T=4` | 25.42 | 1195.4 ms | 164.9 ms | 1375.0 ms | 45.10 ms | 10.88 ms | 59.84 ms |
| `N=32,T=8` | 42.60 | 3956.1 ms | 313.2 ms | 4286.3 ms | 81.16 ms | 10.90 ms | 99.64 ms |
| `N=32,T=16` | 30.88 | 5671.3 ms | 314.2 ms | 6008.2 ms | 151.9 ms | 10.82 ms | 177.4 ms |
| `N=32,T=32` | 15.20 | 5658.2 ms | 359.4 ms | 6052.2 ms | 294.2 ms | 10.90 ms | 334.5 ms |
| `N=100,T=1` | 33.08 | 982.4 ms | 374.7 ms | 1385.5 ms | 57.77 ms | 29.46 ms | 90.48 ms |
| `N=100,T=5` | 27.25 | 4140.9 ms | 422.1 ms | 4600.2 ms | 164.7 ms | 28.55 ms | 208.2 ms |
| `N=100,T=10` | 38.00 | 11.63 s | 722.2 ms | 12.40 s | 306.6 ms | 28.70 ms | 364.3 ms |
| `N=100,T=25` | 36.88 | 27.99 s | 1166.7 ms | 29.24 s | 725.2 ms | 29.35 ms | 827.5 ms |
| `N=100,T=50` | 31.80 | 48.41 s | 1828.9 ms | 50.39 s | 1436.0 ms | 29.08 ms | 1609.1 ms |
| `N=1,T=2` | 18.92 | 114.5 ms | 41.95 ms | 160.5 ms | 2.77 ms | 2.81 ms | 5.77 ms |
| `N=1,T=4` | 24.75 | 295.8 ms | 64.39 ms | 364.5 ms | 3.33 ms | 2.81 ms | 6.37 ms |
| `N=1,T=8` | 33.40 | 786.2 ms | 99.15 ms | 889.6 ms | 4.46 ms | 2.83 ms | 7.65 ms |
| `N=1,T=16` | 26.88 | 1273.2 ms | 117.1 ms | 1395.1 ms | 6.72 ms | 2.84 ms | 10.15 ms |

### Key / signature sizes

| parameter set | sk | pk (single signer) | ring PK table (N·T·pk) | signature |
|---|---:|---:|---:|---:|
| `N=32,T=1` | 32 B | 7.12 KiB | 228.00 KiB | 21.46 KiB |
| `N=32,T=4` | 32 B | 7.12 KiB | 912.00 KiB | 22.53 KiB |
| `N=32,T=8` | 32 B | 7.12 KiB | 1.78 MiB | 22.99 KiB |
| `N=32,T=16` | 32 B | 7.12 KiB | 3.56 MiB | 23.62 KiB |
| `N=32,T=32` | 32 B | 7.12 KiB | 7.12 MiB | 24.07 KiB |
| `N=100,T=1` | 32 B | 7.12 KiB | 712.50 KiB | 32.57 KiB |
| `N=100,T=5` | 32 B | 7.12 KiB | 3.48 MiB | 33.74 KiB |
| `N=100,T=10` | 32 B | 7.12 KiB | 6.96 MiB | 34.36 KiB |
| `N=100,T=25` | 32 B | 7.12 KiB | 17.40 MiB | 35.09 KiB |
| `N=100,T=50` | 32 B | 7.12 KiB | 34.79 MiB | 35.53 KiB |
| `N=1,T=2` | 32 B | 7.12 KiB | 14.25 KiB | 16.84 KiB |
| `N=1,T=4` | 32 B | 7.12 KiB | 28.50 KiB | 17.47 KiB |
| `N=1,T=8` | 32 B | 7.12 KiB | 57.00 KiB | 17.91 KiB |
| `N=1,T=16` | 32 B | 7.12 KiB | 114.00 KiB | 18.54 KiB |

`sk` is the 32-byte seed the signer stores — the full secret-key
material `s ∈ R_q^{l+k}` is deterministically expanded from it.
