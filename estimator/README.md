# LoTRS Parameter Estimator

This directory contains the scripts used to reproduce the concrete
parameter estimates reported for the LoTRS paper, including the
`N = 100`, `T = 50` parameter set.

## Requirements

- SageMath with Python support available as `sage`.
- No network access or external datasets are required at run time.

## Reproducing the Paper Parameter Estimate

From the repository root, run:

```bash
cd estimator
make
```

The scripts import `sage.all`, so running them with plain `python3` will
fail unless that Python environment is the Sage Python environment.

The Makefile invokes:

```bash
sage -c "exec(open('lotrs_estimate.py').read()); main()"
```

This command is used because it runs from this directory, makes the
local helper modules importable, and explicitly invokes the script entry
point.

The script prints:

- MLWE estimates for the binary-proof and DualMS components,
- ASIS/MSIS-style estimates for the binary-proof and DualMS components,
- the concrete signature size,
- the single public-key and full ring public-key sizes,
- the expected number of rejection-sampling repetitions,
- basic parameter-condition checks.

For the paper parameter point, the expected headline values are:

```text
Signature size: about 35.14 KB
Single public key size: about 7.13 KB
Ring PK size: about 35,625 KB
Number of repetitions for rejection sampling: about 12.34
```

The checked-in reference output
`LoTRS-Estimate-Output-N100T50.txt` is intended only as a comparison
log. Before including output logs in an anonymous artifact, regenerate
or scrub them so they do not contain shell prompts, usernames,
hostnames, local paths, timestamps, or other environment metadata.

To write a fresh comparison log, run:

```bash
make reference
```

To remove generated caches and logs, run:

```bash
make clean
```

## Files

- `lotrs_estimate.py`: main script for the concrete LoTRS parameter
  point.
- `lotrs_finder.py`: LoTRS-specific parameter, size, and repetition
  helper formulas.
- `lotrs_param_checks.py`: consistency checks for the chosen moduli,
  challenge differences, regularity bounds, and range-proof condition.
- `LoTRS-Estimate-Output-N100T50.txt`: sample output for the concrete
  paper parameter point.

## Provenance of External Components

This directory includes local copies or adaptations of public estimator
code so the artifact can be evaluated without fetching dependencies:

- `estimator/`: bundled lattice/LWE estimator code used through
  `from estimator import *`.
- `kd_estimates/`: security-estimation helper scripts derived from
  public post-quantum security-estimation code.
- `ASIS_sec_estimate/`: asymmetric-SIS estimation scripts adapted from
  Dilithium/security-estimates scripts, as noted in
  `ASIS_sec_estimate/README.md`.

The LoTRS-specific entry points are `lotrs_estimate.py`,
`lotrs_finder.py`, and `lotrs_param_checks.py`. When preparing a public
artifact, preserve the provenance notes and any applicable license
notices for the bundled external estimator components.
