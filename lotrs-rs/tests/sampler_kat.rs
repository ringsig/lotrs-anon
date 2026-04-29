//! Cross-language Gaussian-sampler KAT.
//!
//! Loads `tests/sampler_kat.json` — emitted by
//! `scripts/gen_sampler_kat.py` against the Python reference — and
//! verifies that the Rust backends produce byte-identical samples for
//! each entry.  The KAT covers:
//!
//! * forced-FACCT at moderate sigma (two seeds)
//! * FACCT at large sigma (the BENCH `sigma_0` regime)
//! * CDT at small sigma
//!
//! The test is skipped gracefully if the KAT file is missing so a
//! fresh clone doesn't fail `cargo test` before the user has generated
//! the file.

use std::fs;
use std::path::PathBuf;

use lotrs::sample::{prepare_facct, xof_sample_gaussian, xof_sample_gaussian_facct, Tag, Xof};

#[derive(serde::Deserialize)]
struct Kat {
    kat: Vec<Entry>,
}

#[derive(serde::Deserialize)]
struct Entry {
    name: String,
    backend: String,
    sigma: f64,
    d: usize,
    seed: String,
    tags: Vec<String>,
    samples: Vec<i64>,
    /// Optional: names a shipped `cdt::CDT_*` constant.  Present on
    /// CDT KAT entries whose table we ship; absent on informational
    /// CDT entries (e.g. on-the-fly small-sigma tables Rust can't
    /// rebuild).
    #[serde(default)]
    cdt_name: Option<String>,
}

fn kat_path() -> Option<PathBuf> {
    let here = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let p = here.join("tests").join("sampler_kat.json");
    if p.exists() {
        Some(p)
    } else {
        None
    }
}

fn build_xof(seed_hex: &str, tag_hex_list: &[String]) -> Xof {
    let seed = hex::decode(seed_hex).expect("seed hex");
    // Convert each tag from hex into a Tag::Bytes slice.  We store the
    // decoded bytes in a Vec and take slice references in a second pass
    // because Tag borrows its contents.
    let tag_bytes: Vec<Vec<u8>> = tag_hex_list
        .iter()
        .map(|t| hex::decode(t).expect("tag hex"))
        .collect();
    let tags: Vec<Tag<'_>> = tag_bytes.iter().map(|b| Tag::Bytes(b)).collect();
    Xof::new(&seed, &tags)
}

#[test]
fn cross_language_sampler_kat_matches_python() {
    let path = match kat_path() {
        Some(p) => p,
        None => {
            eprintln!(
                "tests/sampler_kat.json not found — skipping. \
                       Regenerate via `python scripts/gen_sampler_kat.py > tests/sampler_kat.json`"
            );
            return;
        }
    };
    let raw = fs::read_to_string(path).expect("read KAT");
    let kat: Kat = serde_json::from_str(&raw).expect("parse KAT");

    for entry in kat.kat {
        let mut xof = build_xof(&entry.seed, &entry.tags);
        let rust_samples: Vec<i64> = match entry.backend.as_str() {
            "facct" => {
                let p = prepare_facct(entry.sigma).expect("KAT sigma must be finite positive");
                xof_sample_gaussian_facct(&mut xof, &p, entry.d)
            }
            "cdt" => {
                // Resolve the shipped CDT by `cdt_name` so every CDT
                // fixture is directly cross-language pinned.  Entries
                // without a `cdt_name` are informational (no shipped
                // table) and are skipped with a notice.
                let table_name = match &entry.cdt_name {
                    Some(n) => n.as_str(),
                    None => {
                        eprintln!("skipping CDT KAT {:?} (no cdt_name)", entry.name);
                        continue;
                    }
                };
                let cdt = match lookup_cdt_by_name(table_name) {
                    Some(t) => t,
                    None => panic!(
                        "KAT {:?} references unknown cdt_name {:?}",
                        entry.name, table_name
                    ),
                };
                xof_sample_gaussian(&mut xof, cdt, 128, entry.d)
            }
            other => panic!("unknown backend {other:?} in KAT"),
        };
        assert_eq!(
            rust_samples, entry.samples,
            "KAT {:?} diverged between Rust and Python",
            entry.name
        );
    }
}

/// Resolve a shipped `cdt::CDT_*` constant by its unqualified name.
/// Every entry in `tests/sampler_kat.json` with a `cdt_name` must
/// match one of the arms below.
fn lookup_cdt_by_name(name: &str) -> Option<&'static [u128]> {
    use lotrs::cdt;
    match name {
        "CDT_SIGMA_0_TEST" => Some(cdt::CDT_SIGMA_0_TEST),
        "CDT_SIGMA_0_PRIME_TEST" => Some(cdt::CDT_SIGMA_0_PRIME_TEST),
        "CDT_SIGMA_A_TEST" => Some(cdt::CDT_SIGMA_A_TEST),
        "CDT_SIGMA_B_TEST" => Some(cdt::CDT_SIGMA_B_TEST),
        "CDT_SIGMA_A_BENCH_4OF32" => Some(cdt::CDT_SIGMA_A_BENCH_4OF32),
        "CDT_SIGMA_A_BENCH" => Some(cdt::CDT_SIGMA_A_BENCH),
        "CDT_SIGMA_A_PRODUCTION" => Some(cdt::CDT_SIGMA_A_PRODUCTION),
        "CDT_SIGMA_B_BENCH_PRODUCTION" => Some(cdt::CDT_SIGMA_B_BENCH_PRODUCTION),
        _ => None,
    }
}
