# LoTRS benchmark summary (release build, single core)

## LoTRS combined protocol (threshold ring sig)

| N | T | Sign (s) | Verify (ms) | KAgg (ms) | sig (KiB) | single pk (KiB) | ring PK (MiB) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 32 | 4 | 1.38 | 60 | 35 | 22.53 | 7.12 | 0.89 |
| 32 | 8 | 4.29 | 100 | 71 | 22.99 | 7.12 | 1.78 |
| 32 | 16 | 6.01 | 177 | 144 | 23.62 | 7.12 | 3.56 |
| 32 | 32 | 6.05 | 334 | 288 | 24.07 | 7.12 | 7.12 |
| 100 | 5 | 4.60 | 208 | 134 | 33.74 | 7.12 | 3.48 |
| 100 | 10 | 12.40 | 364 | 276 | 34.36 | 7.12 | 6.96 |
| 100 | 25 | 29.24 | 828 | 702 | 35.09 | 7.12 | 17.40 |
| 100 | 50 | 50.39 | 1609 | 1410 | 35.53 | 7.12 | 34.79 |

## RS-alone (T=1)

Plain ring signature: one signer, ring of N keys.  Numbers are the LoTRS protocol at T=1, *not* re-tuned (φ=11.75·T is small at T=1, so attempt counts are inflated relative to a properly-tuned standalone RS).

| N | Sign (s) | Verify (ms) | sig (KiB) | attempts |
|---:|---:|---:|---:|---:|
| 32 | 0.32 | 31 | 21.46 | 17.7 |
| 100 | 1.39 | 90 | 32.57 | 33.1 |

## DualMS-alone (N=1)

Plain multi-signature, no ring hiding.  `Sign_DualMS` / `Verify_DualMS` are read from the breakdown columns to exclude the phantom β=1 binary proof overhead that still runs in this implementation.

| T | Sign$_\mathrm{DualMS}$ (ms) | Verify$_\mathrm{DualMS}$ (ms) | sig (KiB) | attempts |
|---:|---:|---:|---:|---:|
| 2 | 114 | 2.8 | 16.84 | 18.9 |
| 4 | 296 | 3.3 | 17.47 | 24.8 |
| 8 | 786 | 4.5 | 17.91 | 33.4 |
| 16 | 1273 | 6.7 | 18.54 | 26.9 |

