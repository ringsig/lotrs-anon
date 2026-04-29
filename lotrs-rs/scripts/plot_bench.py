#!/usr/bin/env python3
"""
plot_bench.py -- Parse `examples/bench` output and emit CSV + paper plots.

Inputs
------
A Markdown report produced by `cargo run --release --example bench`
(see the "Benchmarks" section of lotrs-rs/README.md). The report is
parsed from a file path argument or from stdin.

Outputs
-------
Files are written next to the input file when a path is given,
or into the current working directory when reading from stdin:

  data.csv                 one row per (N, T) cell
  summary.md               Markdown table sorted by (N, T) for the paper
  sign-vs-T.pdf            log-log plot: Sign / Verify time vs threshold T
                           (threshold cells only; excludes edge rows)
  sigsize-vs-T.pdf         linear plot: signature size (KiB) vs T
  breakdown-vs-T.pdf       Sign decomposition: Sign_DualMS vs Sign_RS vs T
  rs-alone-vs-N.pdf        T=1 sweep: RS-alone Sign / Verify / sig vs N
  dualms-alone-vs-T.pdf    N=1 sweep: DualMS-alone Sign_DualMS / Verify_DualMS vs T

Usage
-----
  # generate a consolidated grid (threshold + RS-alone + DualMS-alone):
  cargo run --release --example bench -- \\
      --grid "32,1:4:8:16:32;100,1:5:10:25:50;1,2:4:8:16" \\
      > bench-out/grid-d128.md 2>&1

  # parse + plot:
  python3 lotrs-rs/scripts/plot_bench.py bench-out/grid-d128.md

  # or via stdin:
  cat bench-out/grid-d128.md | python3 lotrs-rs/scripts/plot_bench.py

Requirements
------------
  pip install matplotlib

The script accepts the bench binary's default Markdown output verbatim;
it ignores any prelude / postscript text and pulls only the three tables
("Primitive timings" / breakdown / "Key / signature sizes") via regex.
If new columns are added to those tables in `examples/bench.rs`, the
regexes here must be updated to match.
"""

import csv
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Primitive timings: | name | d | N | T | samples | KeyGen | KAgg | Sign | Verify |
HEAD = re.compile(r"^\| `([^`]+)` \| (\d+) \| (\d+) \| (\d+) \| (\d+) "
                  r"\| ([^|]+?) \| ([^|]+?) \| ([^|]+?) \| ([^|]+?) \|$")

# Breakdown: | name | attempts | Sign_DualMS | Sign_RS | Sign | Verify_DualMS | Verify_RS | Verify |
BREAK = re.compile(r"^\| `([^`]+)` \| ([\d.]+) "
                   r"\| ([^|]+?) \| ([^|]+?) \| ([^|]+?) "
                   r"\| ([^|]+?) \| ([^|]+?) \| ([^|]+?) \|$")

# Sizes: | name | sk | pk | ring | sig |
SIZE = re.compile(r"^\| `([^`]+)` \| ([^|]+?) \| ([^|]+?) "
                  r"\| ([^|]+?) \| ([^|]+?) \|$")


def parse_ms(s):
    s = s.strip()
    if s.endswith(" s"):
        return float(s[:-2]) * 1000.0
    if s.endswith(" ms"):
        return float(s[:-3])
    raise ValueError(s)


def parse_bytes(s):
    s = s.strip()
    units = [("MiB", 1024 ** 2), ("KiB", 1024), (" B", 1)]
    for u, k in units:
        if s.endswith(u):
            return int(round(float(s[:-len(u)]) * k))
    raise ValueError(s)


def parse_report(text):
    rows = {}
    for line in text.splitlines():
        m = HEAD.match(line)
        if m:
            name, d, N, T, samples, kg, ka, sg, vf = m.groups()
            rows.setdefault(name, {"name": name})
            rows[name].update({
                "d": int(d), "N": int(N), "T": int(T),
                "samples": int(samples),
                "keygen_ms": parse_ms(kg),
                "kagg_ms": parse_ms(ka),
                "sign_ms": parse_ms(sg),
                "verify_ms": parse_ms(vf),
            })
            continue
        m = BREAK.match(line)
        if m:
            name, attempts, sd, sr, sg, vd, vr, vf = m.groups()
            rows.setdefault(name, {"name": name})
            rows[name].update({
                "attempts": float(attempts),
                "sign_dualms_ms": parse_ms(sd),
                "sign_rs_ms": parse_ms(sr),
                "sign_total_ms": parse_ms(sg),
                "verify_dualms_ms": parse_ms(vd),
                "verify_rs_ms": parse_ms(vr),
                "verify_total_ms": parse_ms(vf),
            })
            continue
        m = SIZE.match(line)
        if m:
            name, sk, pk, ring, sig = m.groups()
            rows.setdefault(name, {"name": name})
            rows[name].update({
                "sk_bytes": parse_bytes(sk),
                "pk_bytes": parse_bytes(pk),
                "ring_bytes": parse_bytes(ring),
                "sig_bytes": parse_bytes(sig),
            })
    return list(rows.values())


def write_csv(rows, path):
    cols = ["name", "d", "N", "T", "samples", "attempts",
            "keygen_ms", "kagg_ms",
            "sign_ms", "sign_dualms_ms", "sign_rs_ms",
            "verify_ms", "verify_dualms_ms", "verify_rs_ms",
            "sk_bytes", "pk_bytes", "ring_bytes", "sig_bytes"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r.get("N", 0), r.get("T", 0))):
            w.writerow({k: r.get(k, "") for k in cols})


# ---------------------------------------------------------------------------
# Row classifiers — the consolidated grid mixes three regimes.
#  - threshold:   T >= 2 and N >= 2  (the LoTRS combined protocol)
#  - rs_alone:    T == 1 and N >= 2  (one signer, ring of N)
#  - dualms_alone: N == 1 and T >= 2 (multi-sig, no ring hiding)
# ---------------------------------------------------------------------------

def is_threshold(r):
    return r.get("N", 0) >= 2 and r.get("T", 0) >= 2


def is_rs_alone(r):
    return r.get("T", 0) == 1 and r.get("N", 0) >= 2


def is_dualms_alone(r):
    return r.get("N", 0) == 1 and r.get("T", 0) >= 2


def split_by_N(rows):
    out = {}
    for r in rows:
        out.setdefault(r["N"], []).append(r)
    for N in out:
        out[N].sort(key=lambda r: r["T"])
    return out


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_combined(rows, path):
    """Sign + Verify vs T for the threshold-protocol rows."""
    grouped = split_by_N([r for r in rows if is_threshold(r)])
    if not grouped:
        return
    fig, ax = plt.subplots(figsize=(5.5, 3.8))
    for N, group in sorted(grouped.items()):
        T = [r["T"] for r in group]
        sg = [r["sign_ms"] / 1000.0 for r in group]
        vf = [r["verify_ms"] / 1000.0 for r in group]
        ax.plot(T, sg, marker="o", linewidth=1.4,
                label=f"Sign,   N={N}")
        ax.plot(T, vf, marker="s", linewidth=1.0, linestyle="--",
                label=f"Verify, N={N}")
    ax.set_xlabel("threshold $T$")
    ax.set_ylabel("time per call (s)")
    ax.set_title("LoTRS primitive timings (Rust, single core)")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize="small", framealpha=0.9, loc="upper left")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_breakdown(rows, path):
    """Sign decomposition: Sign_DualMS vs Sign_RS vs T, faceted by N.

    Each N gets one color; three line styles distinguish total / DualMS / RS.
    """
    grouped = split_by_N([r for r in rows if is_threshold(r)
                          and "sign_dualms_ms" in r])
    if not grouped:
        return
    fig, ax = plt.subplots(figsize=(5.8, 4.0))
    cmap = plt.get_cmap("tab10")
    for idx, (N, group) in enumerate(sorted(grouped.items())):
        color = cmap(idx)
        T = [r["T"] for r in group]
        tot = [r["sign_ms"] / 1000.0 for r in group]
        dm = [r["sign_dualms_ms"] / 1000.0 for r in group]
        rs = [r["sign_rs_ms"] / 1000.0 for r in group]
        ax.plot(T, tot, marker="o", color=color, linewidth=1.6,
                label=f"Sign total,   N={N}")
        ax.plot(T, dm, marker="^", color=color, linewidth=1.1,
                linestyle="--", label=f"Sign DualMS, N={N}")
        ax.plot(T, rs, marker="v", color=color, linewidth=1.1,
                linestyle=":", label=f"Sign RS,     N={N}")
    ax.set_xlabel("threshold $T$")
    ax.set_ylabel("time per call (s)")
    ax.set_title("Sign decomposition: DualMS multi-sig vs RS binary proof")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize="x-small", framealpha=0.9, loc="upper left", ncol=2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_sigsize(rows, path):
    """Signature size vs T for threshold rows."""
    grouped = split_by_N([r for r in rows if is_threshold(r)])
    if not grouped:
        return
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    for N, group in sorted(grouped.items()):
        T = [r["T"] for r in group]
        sz = [r["sig_bytes"] / 1024.0 for r in group]
        ax.plot(T, sz, marker="o", linewidth=1.4, label=f"N = {N}")
    ax.set_xlabel("threshold $T$")
    ax.set_ylabel("signature size (KiB)")
    ax.set_title("LoTRS signature size vs threshold")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize="small", framealpha=0.9)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_rs_alone(rows, path):
    """T=1 sweep: RS-alone Sign / Verify / sig-size vs ring size N."""
    rs_rows = sorted([r for r in rows if is_rs_alone(r)],
                     key=lambda r: r["N"])
    if not rs_rows:
        return
    fig, (ax_t, ax_s) = plt.subplots(1, 2, figsize=(8.4, 3.4))
    N = [r["N"] for r in rs_rows]
    sg = [r["sign_ms"] / 1000.0 for r in rs_rows]
    vf = [r["verify_ms"] for r in rs_rows]
    sz = [r["sig_bytes"] / 1024.0 for r in rs_rows]
    ax_t.plot(N, sg, marker="o", linewidth=1.5, label="Sign (s)")
    ax_t.plot(N, vf, marker="s", linewidth=1.2, linestyle="--",
              label="Verify (ms)")
    ax_t.set_xlabel("ring size $N$")
    ax_t.set_ylabel("time")
    ax_t.set_title("RS-alone (T=1): time vs N")
    ax_t.set_yscale("log")
    ax_t.set_xscale("log")
    ax_t.grid(True, which="both", alpha=0.3)
    ax_t.legend(fontsize="small", framealpha=0.9)
    ax_s.plot(N, sz, marker="o", linewidth=1.5)
    ax_s.set_xlabel("ring size $N$")
    ax_s.set_ylabel("signature size (KiB)")
    ax_s.set_title("RS-alone (T=1): sig size vs N")
    ax_s.set_xscale("log")
    ax_s.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_dualms_alone(rows, path):
    """N=1 sweep: DualMS-alone Sign_DualMS / Verify_DualMS vs T.

    Reads the DualMS-only sub-times to exclude phantom β=1 binary-proof
    overhead from the headline numbers.
    """
    dm_rows = sorted([r for r in rows if is_dualms_alone(r)
                      and "sign_dualms_ms" in r],
                     key=lambda r: r["T"])
    if not dm_rows:
        return
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    T = [r["T"] for r in dm_rows]
    sg = [r["sign_dualms_ms"] / 1000.0 for r in dm_rows]
    vf = [r["verify_dualms_ms"] for r in dm_rows]
    ax.plot(T, sg, marker="o", linewidth=1.5, label="Sign$_\\mathrm{DualMS}$ (s)")
    ax.plot(T, vf, marker="s", linewidth=1.2, linestyle="--",
            label="Verify$_\\mathrm{DualMS}$ (ms)")
    ax.set_xlabel("number of signers $T$")
    ax.set_ylabel("time")
    ax.set_title("DualMS-alone (N=1): time vs T")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize="small", framealpha=0.9, loc="upper left")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def emit_summary_table(rows, path):
    """Threshold-protocol cells in one table; standalone cells split out."""
    rows_sorted = sorted(rows, key=lambda r: (r.get("N", 0), r.get("T", 0)))
    threshold = [r for r in rows_sorted if is_threshold(r)]
    rs_alone = [r for r in rows_sorted if is_rs_alone(r)]
    dualms_alone = [r for r in rows_sorted if is_dualms_alone(r)]
    with open(path, "w") as f:
        f.write("# LoTRS benchmark summary (release build, single core)\n\n")

        if threshold:
            f.write("## LoTRS combined protocol (threshold ring sig)\n\n")
            f.write("| N | T | Sign (s) | Verify (ms) | KAgg (ms) "
                    "| sig (KiB) | single pk (KiB) | ring PK (MiB) |\n")
            f.write("|---:|---:|---:|---:|---:|---:|---:|---:|\n")
            for r in threshold:
                ring_mib = r["ring_bytes"] / (1024 ** 2)
                f.write(f"| {r['N']} | {r['T']} | "
                        f"{r['sign_ms']/1000:.2f} | "
                        f"{r['verify_ms']:.0f} | "
                        f"{r['kagg_ms']:.0f} | "
                        f"{r['sig_bytes']/1024:.2f} | "
                        f"{r['pk_bytes']/1024:.2f} | "
                        f"{ring_mib:.2f} |\n")
            f.write("\n")

        if rs_alone:
            f.write("## RS-alone (T=1)\n\n")
            f.write("Plain ring signature: one signer, ring of N keys.  "
                    "Numbers are the LoTRS protocol at T=1, *not* re-tuned "
                    "(φ=11.75·T is small at T=1, so attempt counts are "
                    "inflated relative to a properly-tuned standalone RS).\n\n")
            f.write("| N | Sign (s) | Verify (ms) | sig (KiB) | attempts |\n")
            f.write("|---:|---:|---:|---:|---:|\n")
            for r in rs_alone:
                f.write(f"| {r['N']} | "
                        f"{r['sign_ms']/1000:.2f} | "
                        f"{r['verify_ms']:.0f} | "
                        f"{r['sig_bytes']/1024:.2f} | "
                        f"{r.get('attempts', float('nan')):.1f} |\n")
            f.write("\n")

        if dualms_alone:
            f.write("## DualMS-alone (N=1)\n\n")
            f.write("Plain multi-signature, no ring hiding.  "
                    "`Sign_DualMS` / `Verify_DualMS` are read from the "
                    "breakdown columns to exclude the phantom β=1 binary "
                    "proof overhead that still runs in this implementation.\n\n")
            f.write("| T | Sign$_\\mathrm{DualMS}$ (ms) "
                    "| Verify$_\\mathrm{DualMS}$ (ms) "
                    "| sig (KiB) | attempts |\n")
            f.write("|---:|---:|---:|---:|---:|\n")
            for r in dualms_alone:
                f.write(f"| {r['T']} | "
                        f"{r['sign_dualms_ms']:.0f} | "
                        f"{r['verify_dualms_ms']:.1f} | "
                        f"{r['sig_bytes']/1024:.2f} | "
                        f"{r.get('attempts', float('nan')):.1f} |\n")
            f.write("\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) > 1:
        text = Path(sys.argv[1]).read_text()
    else:
        text = sys.stdin.read()
    rows = parse_report(text)
    if not rows:
        print("no data parsed", file=sys.stderr)
        sys.exit(2)
    out = Path(sys.argv[1]).parent if len(sys.argv) > 1 else Path(".")
    out.mkdir(exist_ok=True)
    write_csv(rows, out / "data.csv")
    plot_combined(rows, out / "sign-vs-T.pdf")
    plot_sigsize(rows, out / "sigsize-vs-T.pdf")
    plot_breakdown(rows, out / "breakdown-vs-T.pdf")
    plot_rs_alone(rows, out / "rs-alone-vs-N.pdf")
    plot_dualms_alone(rows, out / "dualms-alone-vs-T.pdf")
    emit_summary_table(rows, out / "summary.md")
    print(f"wrote {out}/data.csv, summary.md, "
          f"sign-vs-T.pdf, sigsize-vs-T.pdf, "
          f"breakdown-vs-T.pdf, rs-alone-vs-N.pdf, dualms-alone-vs-T.pdf")
    print(f"parsed {len(rows)} cells: {[r['name'] for r in rows]}")


if __name__ == "__main__":
    main()
