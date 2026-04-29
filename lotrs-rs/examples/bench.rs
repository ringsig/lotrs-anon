//! Simple benchmarks for the LoTRS high-level primitives.
//!
//! Measures wall-clock time for every public primitive alongside the
//! concrete public-key / secret-key / signature byte sizes.  Output is
//! a Markdown table so the numbers can be pasted into artifact notes.
//!
//! Usage:
//!
//!   cargo run --release --example bench                      # TEST + BENCH_4OF32 + BENCH_PARAMS
//!   cargo run --release --example bench -- --with-prod       # + PRODUCTION_PARAMS (slow)
//!   cargo run --release --example bench -- --skip-test
//!   cargo run --release --example bench -- --grid 25,5:13:25;50,10:25:50;100,25:50
//!
//! `--grid` accepts a `;`-separated list of `N,T1:T2:…` groups.  Each
//! group shares a ring size `N` and instantiates one signer set per
//! listed threshold `T`.  All grid entries use the d=128 lattice
//! declared by the paper (`k=12, l=5, l'=6, n̂=10, k̂=8, phi_a=50,
//! phi=11.75·T`, `q = largest prime ≤ 2^38`, `q_hat = largest prime
//! ≤ 2^33`, both ≡ 5 mod 8; `mask_sampler = facct`).

use std::env;
use std::time::{Duration, Instant};

use lotrs::lotrs::LoTRS;
use lotrs::{
    LoTRSCodec, LoTRSParams, MaskSamplerKind, BENCH_4OF32, BENCH_PARAMS, PRODUCTION_PARAMS,
    TEST_PARAMS,
};

// -------------------------------------------------------------------------

struct Row {
    name: String,
    par: LoTRSParams,
    keygen: Duration,
    kagg: Duration,
    sign: Duration,
    /// DualMS portion of `sign` — sign1 (T-fold commitments) +
    /// sign2 minus the binary proof + sagg.  Summed across all
    /// rejection-sampling attempts in a single call, then averaged
    /// across seeds.
    sign_dualms: Duration,
    /// Ring-signature portion of `sign` — `sign_bin` only.  Summed
    /// across all attempts then averaged across seeds.
    sign_rs: Duration,
    /// Mean number of rejection-sampling attempts per accepted
    /// signature.  Mostly a sanity check (μ_total ≈ 12 at production).
    sign_attempts: f64,
    verify: Duration,
    /// DualMS portion of `verify` — z̃/r̃/ẽ bounds, A·z̃ + B·r̃ + ẽ
    /// reconstruction, KAgg, ring-keyed pk_sum, w̃₀ recovery, FS hash.
    verify_dualms: Duration,
    /// Ring-signature portion of `verify` — f1/f0/g0/g1 bounds,
    /// A_hat_bin reconstruction & low-bit check.
    verify_rs: Duration,
    pk_bytes: usize,
    sk_bytes: usize,
    sig_bytes: usize,
    sign_samples: usize,
}

/// Build a d=128 grid parameter set for an arbitrary `(N, T)`.  Shares
/// the lattice, bit-drops, and `phi_a / phi_b / eps_tot` with the
/// paper-aligned `PRODUCTION_PARAMS`; only `beta = N`, `T`, and
/// `phi = 11.75·T` change.  The `name` field borrows from the
/// `lotrs-bench-16of32` whitelist entry so `resolve_cdt` picks the
/// right (identical across 4of32 / 16of32 / 50of100) `sigma_a` CDT.
fn make_grid_params(n: usize, t: usize) -> LoTRSParams {
    LoTRSParams {
        name: "lotrs-bench-16of32",
        d: 128,
        q: 274_877_906_837,
        q_hat: 8_589_934_237,
        kappa: 1,
        beta: n, // κ = 1 ⇒ N = β
        T: t,
        k: 12,
        l: 5,
        l_prime: 6,
        n_hat: 10,
        k_hat: 8,
        w: 31,
        eta: 1,
        phi: 11.75 * t as f64,
        phi_a: 50.0,
        phi_b: 4.0,
        K_A: 20,
        K_B: 5,
        K_w: 5,
        lam: 128,
        max_attempts: 200,
        eta_prime: -1,
        tail_t: 1.2,
        mask_sampler: MaskSamplerKind::Facct,
        eps_tot: 0.01,
    }
}

fn bench_one(name: String, par: LoTRSParams) -> Row {
    let scheme = LoTRS::new(par);
    let codec = LoTRSCodec::new(par);
    let pp = scheme.setup(&[0u8; 32]);

    // ---- KeyGen — average over a few runs ----------------------------
    let keygen_samples = if par.d == 32 { 64 } else { 4 };
    let t0 = Instant::now();
    for i in 0..keygen_samples {
        let mut seed = [0u8; 32];
        seed[0] = i as u8;
        let _ = scheme.keygen(&pp, &seed);
    }
    let keygen = t0.elapsed() / keygen_samples as u32;

    // Build the full PK/SK table deterministically.
    let mut pk_table: Vec<Vec<Vec<Vec<u64>>>> = Vec::with_capacity(par.N());
    let mut sk_table: Vec<Vec<Vec<Vec<u64>>>> = Vec::with_capacity(par.N());
    let mut pk_bytes_tbl: Vec<Vec<Vec<u8>>> = Vec::with_capacity(par.N());
    for col in 0..par.N() {
        let mut col_pk = Vec::with_capacity(par.T);
        let mut col_sk = Vec::with_capacity(par.T);
        let mut col_pkb = Vec::with_capacity(par.T);
        for row in 0..par.T {
            let mut seed = [0u8; 32];
            seed[0] = col as u8;
            seed[1] = row as u8;
            let (sk, pk) = scheme.keygen(&pp, &seed);
            let pkb = codec.pk_encode(&pk).expect("pk_encode");
            col_pk.push(pk);
            col_sk.push(sk);
            col_pkb.push(pkb);
        }
        pk_table.push(col_pk);
        sk_table.push(col_sk);
        pk_bytes_tbl.push(col_pkb);
    }

    // ---- KAgg --------------------------------------------------------
    let t0 = Instant::now();
    let _ = scheme.kagg(&pk_table);
    let kagg = t0.elapsed();

    // ---- Sign & Verify — average over multiple signing seeds --------
    // One attempt of sign() takes essentially the same wall-clock
    // regardless of seed, but the number of rejection-sampling attempts
    // until first acceptance is geometrically distributed with mean
    // μ_total ≈ 12 (per the estimator).  A single sample can vary
    // by 3-5× around the mean, so we average several runs.
    // Sample counts sized to keep a full d=128 grid (through T=50 at
    // N=100) under ~25 min total wall-clock while smoothing rejection-
    // sampling attempt variance below ~10% for each cell.  Lower-T rows
    // have cheaper individual calls and get proportionally more reps.
    let sign_samples = match par.d {
        32 => 32, // TEST — cheap
        _ if par.T <= 5 => 12,
        _ if par.T <= 10 => 10,
        _ if par.T <= 25 => 8,
        _ => 5, // T >= 50
    };
    let ell = 0;
    let sks: Vec<_> = (0..par.T).map(|u| sk_table[ell][u].clone()).collect();
    let mu = b"bench";
    let mut last_sig: Vec<u8> = Vec::new();
    let mut sign_dualms_sum = Duration::ZERO;
    let mut sign_rs_sum = Duration::ZERO;
    let mut total_attempts: u64 = 0;
    let t0 = Instant::now();
    for i in 0..sign_samples {
        let mut seed = [0xaa; 32];
        seed[1] = i as u8;
        let (sig_bytes, st) = scheme
            .sign_with_timings(&pp, &sks, ell, mu, &pk_table, &seed)
            .expect("sign");
        last_sig = sig_bytes;
        // DualMS = sign1 + sign2_rest + sagg; RS = sign_bin.
        sign_dualms_sum += st.sign1 + st.sign2_rest + st.sagg;
        sign_rs_sum += st.sign_bin;
        total_attempts += st.attempts as u64;
    }
    let sign = t0.elapsed() / sign_samples as u32;
    let sign_dualms = sign_dualms_sum / sign_samples as u32;
    let sign_rs = sign_rs_sum / sign_samples as u32;
    let sign_attempts = total_attempts as f64 / sign_samples as f64;

    // ---- Verify — average too; verify() cost is deterministic given
    // pp/μ/PK but the same sample hedge keeps numbers stable under load.
    let mut verify_dualms_sum = Duration::ZERO;
    let mut verify_rs_sum = Duration::ZERO;
    let t0 = Instant::now();
    for _ in 0..sign_samples {
        let (ok, vt) = scheme.verify_with_timings(&pp, mu, &last_sig, &pk_bytes_tbl);
        assert!(ok, "bench signature failed to verify");
        verify_dualms_sum += vt.verify_dualms;
        verify_rs_sum += vt.verify_bin;
    }
    let verify = t0.elapsed() / sign_samples as u32;
    let verify_dualms = verify_dualms_sum / sign_samples as u32;
    let verify_rs = verify_rs_sum / sign_samples as u32;
    let sig_bytes_vec = last_sig;

    Row {
        name,
        par,
        keygen,
        kagg,
        sign,
        sign_dualms,
        sign_rs,
        sign_attempts,
        verify,
        verify_dualms,
        verify_rs,
        pk_bytes: pk_bytes_tbl[0][0].len(),
        sk_bytes: 32, // codec::sk_encode output
        sig_bytes: sig_bytes_vec.len(),
        sign_samples,
    }
}

// -------------------------------------------------------------------------

fn fmt_ms(d: Duration) -> String {
    let ms = d.as_secs_f64() * 1000.0;
    if ms < 1.0 {
        format!("{:.3} ms", ms)
    } else if ms < 100.0 {
        format!("{:.2} ms", ms)
    } else if ms < 10_000.0 {
        format!("{:.1} ms", ms)
    } else {
        format!("{:.2} s", ms / 1000.0)
    }
}

fn fmt_bytes(n: usize) -> String {
    if n < 1024 {
        format!("{} B", n)
    } else if n < 1024 * 1024 {
        format!("{:.2} KiB", n as f64 / 1024.0)
    } else {
        format!("{:.2} MiB", n as f64 / (1024.0 * 1024.0))
    }
}

fn print_report(rows: &[Row]) {
    println!();
    println!("### Primitive timings");
    println!();
    println!("Sign / Verify are arithmetic means over the listed");
    println!("number of signing seeds; KeyGen / KAgg are single-run");
    println!("(deterministic given pp / PK).");
    println!();
    println!("| parameter set | d | N | T | samples | KeyGen | KAgg | Sign | Verify |");
    println!("|---|---:|---:|---:|---:|---:|---:|---:|---:|");
    for r in rows {
        println!(
            "| `{}` | {} | {} | {} | {} | {} | {} | {} | {} |",
            r.name,
            r.par.d,
            r.par.N(),
            r.par.T,
            r.sign_samples,
            fmt_ms(r.keygen),
            fmt_ms(r.kagg),
            fmt_ms(r.sign),
            fmt_ms(r.verify),
        );
    }

    println!();
    println!("### Sign / Verify breakdown — DualMS multi-sig vs RS binary proof");
    println!();
    println!("Sign_DualMS = sign1 (T-fold commitments) + sign2 minus the");
    println!("binary proof + sagg.  Sign_RS = `sign_bin`, the binary ring");
    println!("proof, summed over the T independent signers and over every");
    println!("rejection-sampling attempt.  Both are wall-clock totals of");
    println!("the sequential reference implementation, so");
    println!("Sign_DualMS + Sign_RS ≈ Sign minus the one-shot context");
    println!("expansion (expand_A/G/B, NTT prep, KAgg α_u).  `attempts` is");
    println!("the mean attempt count per accepted signature.");
    println!();
    println!("Note: `sign_bin` runs once per signer per attempt in this");
    println!("implementation (`sagg` asserts the T proofs match).  An");
    println!("optimised multi-signer protocol that broadcasts a single");
    println!("`pi` would have standalone RS cost ≈ Sign_RS / T.");
    println!();
    println!("Verify_DualMS = z̃/r̃/ẽ bounds + `A·z̃ + B·r̃ + ẽ` reconstruction");
    println!("+ KAgg + ring-keyed `pk_sum` + w̃₀ recovery + closing FS hash.");
    println!("Verify_RS = f1/f0/g0/g1 bounds + A_hat_bin reconstruction &");
    println!("low-bit check.  pk_sum is attributed to DualMS because it's");
    println!("part of the LHS=RHS multi-sig closing check.");
    println!();
    println!(
        "| parameter set | attempts | Sign_DualMS | Sign_RS | Sign | Verify_DualMS | Verify_RS | Verify |"
    );
    println!("|---|---:|---:|---:|---:|---:|---:|---:|");
    for r in rows {
        println!(
            "| `{}` | {:.2} | {} | {} | {} | {} | {} | {} |",
            r.name,
            r.sign_attempts,
            fmt_ms(r.sign_dualms),
            fmt_ms(r.sign_rs),
            fmt_ms(r.sign),
            fmt_ms(r.verify_dualms),
            fmt_ms(r.verify_rs),
            fmt_ms(r.verify),
        );
    }

    println!();
    println!("### Key / signature sizes");
    println!();
    println!("| parameter set | sk | pk (single signer) | ring PK table (N·T·pk) | signature |");
    println!("|---|---:|---:|---:|---:|");
    for r in rows {
        let ring_bytes = r.pk_bytes * r.par.N() * r.par.T;
        println!(
            "| `{}` | {} | {} | {} | {} |",
            r.name,
            fmt_bytes(r.sk_bytes),
            fmt_bytes(r.pk_bytes),
            fmt_bytes(ring_bytes),
            fmt_bytes(r.sig_bytes),
        );
    }
    println!();
    println!("`sk` is the 32-byte seed the signer stores — the full secret-key");
    println!("material `s ∈ R_q^{{l+k}}` is deterministically expanded from it.");
}

/// Parse `--grid N1,T1:T2:T3;N2,Ta:Tb` into a list of `(N, T)` pairs.
fn parse_grid(spec: &str) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for group in spec.split(';').filter(|s| !s.is_empty()) {
        let (n_s, ts_s) = group
            .split_once(',')
            .unwrap_or_else(|| panic!("grid group {group:?} is not N,T1:T2:…"));
        let n: usize = n_s
            .trim()
            .parse()
            .unwrap_or_else(|_| panic!("bad N in {group:?}"));
        for t_s in ts_s.split(':') {
            let t: usize = t_s
                .trim()
                .parse()
                .unwrap_or_else(|_| panic!("bad T in {group:?}"));
            out.push((n, t));
        }
    }
    out
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let with_prod = args.iter().any(|a| a == "--with-prod");
    let skip_test = args.iter().any(|a| a == "--skip-test");
    let grid_spec: Option<&str> = args
        .iter()
        .position(|a| a == "--grid")
        .and_then(|i| args.get(i + 1).map(|s| s.as_str()));

    let mut targets: Vec<(String, LoTRSParams)> = Vec::new();

    if let Some(spec) = grid_spec {
        let grid = parse_grid(spec);
        for (n, t) in grid {
            let name = format!("N={n},T={t}");
            targets.push((name, make_grid_params(n, t)));
        }
    } else {
        if !skip_test {
            targets.push(("TEST".into(), TEST_PARAMS));
        }
        targets.push(("BENCH_4OF32".into(), BENCH_4OF32));
        targets.push(("BENCH_PARAMS".into(), BENCH_PARAMS));
        if with_prod {
            targets.push(("PRODUCTION".into(), PRODUCTION_PARAMS));
        }
    }

    println!("LoTRS — primitive benchmarks (release build)");
    println!("Timings averaged over multiple signing seeds to smooth out");
    println!("the high-variance rejection-sampling attempt count.");
    println!();
    if grid_spec.is_some() {
        let par0 = targets[0].1;
        println!("All rows use the shared d=128 lattice:");
        println!(
            "  d={}, κ=1, k={}, l={}, l'={}, n̂={}, k̂={}, w={}, η={}",
            par0.d, par0.k, par0.l, par0.l_prime, par0.n_hat, par0.k_hat, par0.w, par0.eta
        );
        println!("  q = {} (largest prime ≤ 2^38 with q ≡ 5 mod 8)", par0.q);
        println!(
            "  q_hat = {} (largest prime ≤ 2^33 with q_hat ≡ 5 mod 8)",
            par0.q_hat
        );
        println!(
            "  phi_a = {}, phi_b = {}, phi = 11.75·T",
            par0.phi_a, par0.phi_b
        );
        println!(
            "  K_A = {}, K_B = {}, K_w = {}, eps_tot = {}",
            par0.K_A, par0.K_B, par0.K_w, par0.eps_tot
        );
        println!("  mask_sampler = facct, tail_t = {}", par0.tail_t);
    }

    let mut rows = Vec::new();
    for (name, par) in targets {
        eprintln!("running {name} ({}-of-{})…", par.T, par.N());
        let r = bench_one(name.clone(), par);
        eprintln!(
            "  sign {} (DualMS {} + RS {}), verify {} (DualMS {} + RS {}), sig {}",
            fmt_ms(r.sign),
            fmt_ms(r.sign_dualms),
            fmt_ms(r.sign_rs),
            fmt_ms(r.verify),
            fmt_ms(r.verify_dualms),
            fmt_ms(r.verify_rs),
            fmt_bytes(r.sig_bytes),
        );
        rows.push(r);
    }
    print_report(&rows);
}
