//! Byte-exact interoperability test against `../lotrs-py/vectors.json`.
//!
//! Runs the Rust implementation with the same deterministic seeds used
//! by the Python reference and asserts that every serialized artefact
//! matches byte-for-byte.  This test is the real guarantor of
//! test-vector compatibility.
//!
//! The integration test is skipped gracefully if the vectors file is
//! missing (so that a fresh checkout doesn't fail `cargo test` before
//! the user has generated the file).

use std::fs;
use std::path::PathBuf;

use lotrs::lotrs::LoTRS;
use lotrs::{LoTRSCodec, TEST_PARAMS};

#[derive(serde::Deserialize)]
struct Vectors {
    schema_version: u32,
    params: String,
    d: usize,
    #[serde(rename = "N")]
    n: usize,
    #[serde(rename = "T")]
    t: usize,
    pp_seed: String,
    pp_bytes: String,
    signing_seed: String,
    ell: usize,
    message: String,
    keygen: Vec<KeygenEntry>,
    signature: SigEntry,
}

#[derive(serde::Deserialize)]
struct KeygenEntry {
    col: usize,
    row: usize,
    seed: String,
    sk_bytes: String,
    pk_bytes: String,
}

#[derive(serde::Deserialize)]
struct SigEntry {
    bytes: String,
    byte_length: usize,
}

fn vectors_path() -> Option<PathBuf> {
    // `CARGO_MANIFEST_DIR` is set to lotrs-rs/ when running tests.
    let here = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let p = here.parent()?.join("lotrs-py").join("vectors.json");
    if p.exists() {
        Some(p)
    } else {
        None
    }
}

fn load_vectors() -> Option<Vectors> {
    let p = vectors_path()?;
    let raw = fs::read_to_string(p).ok()?;
    serde_json::from_str(&raw).ok()
}

#[test]
fn pp_and_keypairs_are_byte_identical() {
    let v = match load_vectors() {
        Some(v) => v,
        None => {
            eprintln!(
                "vectors.json not found — skipping interop test. \
                       Generate it with `python vectors.py --out vectors.json` \
                       in lotrs-py/."
            );
            return;
        }
    };
    assert_eq!(v.schema_version, 1, "vectors schema version");
    assert_eq!(v.params, TEST_PARAMS.name);
    assert_eq!(v.d, TEST_PARAMS.d);
    assert_eq!(v.n, TEST_PARAMS.N());
    assert_eq!(v.t, TEST_PARAMS.T);

    let scheme = LoTRS::new(TEST_PARAMS);
    let codec = LoTRSCodec::new(TEST_PARAMS);

    // ---- pp  round-trip ------------------------------------------------
    let pp_seed = hex::decode(&v.pp_seed).expect("pp_seed hex");
    let pp = scheme.setup(&pp_seed);
    let pp_enc = codec.pp_encode(&pp).expect("pp_encode");
    let expected_pp_bytes = hex::decode(&v.pp_bytes).expect("pp_bytes hex");
    assert_eq!(pp_enc, expected_pp_bytes, "pp bytes mismatch");

    // ---- per-keypair sk / pk  byte-exact -------------------------------
    for kg in &v.keygen {
        let seed = hex::decode(&kg.seed).expect("seed hex");
        let expected_sk = hex::decode(&kg.sk_bytes).expect("sk hex");
        let expected_pk = hex::decode(&kg.pk_bytes).expect("pk hex");

        // sk bytes
        let sk_enc = codec.sk_encode(&seed).expect("sk_encode");
        assert_eq!(
            sk_enc, expected_sk,
            "sk byte mismatch at col={}, row={}",
            kg.col, kg.row
        );

        // full keygen (pp, seed) → (sk_object, pk_object)
        let (_, pk) = scheme.keygen(&pp, &seed);
        let pk_enc = codec.pk_encode(&pk).expect("pk_encode");
        assert_eq!(
            pk_enc, expected_pk,
            "pk byte mismatch at col={}, row={}",
            kg.col, kg.row
        );
    }
}

#[test]
fn signature_decodes_and_reencodes_to_same_bytes() {
    let v = match load_vectors() {
        Some(v) => v,
        None => {
            eprintln!("skipping (no vectors.json)");
            return;
        }
    };
    let codec = LoTRSCodec::new(TEST_PARAMS);
    let sig_bytes = hex::decode(&v.signature.bytes).expect("sig hex");
    assert_eq!(sig_bytes.len(), v.signature.byte_length);

    let sig = codec.sig_decode(&sig_bytes).expect("sig_decode");
    let reencoded = codec.sig_encode(&sig).expect("sig_encode");
    assert_eq!(
        reencoded, sig_bytes,
        "canonical form broken: decoded-then-encoded bytes differ"
    );
}

#[test]
fn signature_decoder_rejects_truncated() {
    let v = match load_vectors() {
        Some(v) => v,
        None => {
            eprintln!("skipping (no vectors.json)");
            return;
        }
    };
    let codec = LoTRSCodec::new(TEST_PARAMS);
    let sig_bytes = hex::decode(&v.signature.bytes).expect("sig hex");
    // drop the last 4 bytes → should fail with a codec error, never panic
    let truncated = &sig_bytes[..sig_bytes.len() - 4];
    assert!(codec.sig_decode(truncated).is_err());
}

#[test]
fn signature_decoder_rejects_trailing_bytes() {
    let v = match load_vectors() {
        Some(v) => v,
        None => {
            eprintln!("skipping (no vectors.json)");
            return;
        }
    };
    let codec = LoTRSCodec::new(TEST_PARAMS);
    let mut sig_bytes = hex::decode(&v.signature.bytes).expect("sig hex");
    sig_bytes.push(0x00);
    assert!(matches!(
        codec.sig_decode(&sig_bytes),
        Err(lotrs::CodecError::TrailingBytes)
    ));
}

#[test]
fn python_signature_verifies_in_rust() {
    let v = match load_vectors() {
        Some(v) => v,
        None => {
            eprintln!("skipping (no vectors.json)");
            return;
        }
    };

    let scheme = LoTRS::new(TEST_PARAMS);
    let pp_seed = hex::decode(&v.pp_seed).expect("pp_seed hex");
    let pp = scheme.setup(&pp_seed);
    let mu = v.message.as_bytes();
    let sig_bytes = hex::decode(&v.signature.bytes).expect("sig hex");

    // rebuild pk_table[col][row] = pk_bytes
    let mut pk_table_bytes: Vec<Vec<Vec<u8>>> = vec![Vec::new(); v.n];
    for kg in &v.keygen {
        let pk_bytes = hex::decode(&kg.pk_bytes).expect("pk hex");
        let col = &mut pk_table_bytes[kg.col];
        while col.len() <= kg.row {
            col.push(Vec::new());
        }
        col[kg.row] = pk_bytes;
    }

    assert!(
        scheme.verify(&pp, mu, &sig_bytes, &pk_table_bytes),
        "verify failed on an honest Python-generated signature"
    );
}

#[test]
fn verify_rejects_tampered_message() {
    let v = match load_vectors() {
        Some(v) => v,
        None => {
            eprintln!("skipping (no vectors.json)");
            return;
        }
    };

    let scheme = LoTRS::new(TEST_PARAMS);
    let pp_seed = hex::decode(&v.pp_seed).expect("pp_seed hex");
    let pp = scheme.setup(&pp_seed);
    let sig_bytes = hex::decode(&v.signature.bytes).expect("sig hex");

    let mut pk_table_bytes: Vec<Vec<Vec<u8>>> = vec![Vec::new(); v.n];
    for kg in &v.keygen {
        let pk_bytes = hex::decode(&kg.pk_bytes).expect("pk hex");
        let col = &mut pk_table_bytes[kg.col];
        while col.len() <= kg.row {
            col.push(Vec::new());
        }
        col[kg.row] = pk_bytes;
    }

    // wrong message
    let wrong_mu = b"not the original message";
    assert!(!scheme.verify(&pp, wrong_mu, &sig_bytes, &pk_table_bytes));
}

#[test]
fn verify_rejects_flipped_signature_byte() {
    let v = match load_vectors() {
        Some(v) => v,
        None => {
            eprintln!("skipping (no vectors.json)");
            return;
        }
    };

    let scheme = LoTRS::new(TEST_PARAMS);
    let pp_seed = hex::decode(&v.pp_seed).expect("pp_seed hex");
    let pp = scheme.setup(&pp_seed);
    let mu = v.message.as_bytes();
    let mut sig_bytes = hex::decode(&v.signature.bytes).expect("sig hex");

    // flip one bit in the middle
    sig_bytes[100] ^= 0x01;

    let mut pk_table_bytes: Vec<Vec<Vec<u8>>> = vec![Vec::new(); v.n];
    for kg in &v.keygen {
        let pk_bytes = hex::decode(&kg.pk_bytes).expect("pk hex");
        let col = &mut pk_table_bytes[kg.col];
        while col.len() <= kg.row {
            col.push(Vec::new());
        }
        col[kg.row] = pk_bytes;
    }

    assert!(!scheme.verify(&pp, mu, &sig_bytes, &pk_table_bytes));
}

#[test]
fn rust_sign_produces_same_bytes_as_python() {
    let v = match load_vectors() {
        Some(v) => v,
        None => {
            eprintln!("skipping (no vectors.json)");
            return;
        }
    };
    let scheme = LoTRS::new(TEST_PARAMS);
    let pp_seed = hex::decode(&v.pp_seed).expect("pp_seed hex");
    let pp = scheme.setup(&pp_seed);
    let signing_seed = hex::decode(&v.signing_seed).expect("signing_seed hex");
    let mu = v.message.as_bytes();

    // rebuild pk_table[col][row] = pk (decoded) and sks = sks[col][row] = decoded sk
    let codec = lotrs::LoTRSCodec::new(TEST_PARAMS);
    let par = TEST_PARAMS;
    let mut pk_table: Vec<Vec<Vec<Vec<u64>>>> = vec![vec![]; v.n];
    let mut sk_table: Vec<Vec<Vec<Vec<u64>>>> = vec![vec![]; v.n];
    for kg in &v.keygen {
        let seed = hex::decode(&kg.seed).expect("seed hex");
        let pk = codec
            .pk_decode(&hex::decode(&kg.pk_bytes).expect("pk hex"))
            .expect("pk_decode");
        // sk is the keygen output (object form), rebuilt deterministically
        let (sk, _pk_check) = scheme.keygen(&pp, &seed);
        let col = &mut pk_table[kg.col];
        let col_sk = &mut sk_table[kg.col];
        while col.len() <= kg.row {
            col.push(vec![]);
            col_sk.push(vec![]);
        }
        col[kg.row] = pk;
        col_sk[kg.row] = sk;
    }

    // sks for the signing column ell
    let sks: Vec<Vec<Vec<u64>>> = (0..par.T).map(|u| sk_table[v.ell][u].clone()).collect();

    let sig_bytes = scheme
        .sign(&pp, &sks, v.ell, mu, &pk_table, &signing_seed)
        .expect("Rust sign should succeed");

    let expected = hex::decode(&v.signature.bytes).expect("sig hex");
    assert_eq!(
        sig_bytes.len(),
        expected.len(),
        "signature byte lengths differ: rust={} vs python={}",
        sig_bytes.len(),
        expected.len()
    );
    assert_eq!(sig_bytes, expected, "signature bytes differ");
}

#[test]
fn rust_sign_output_verifies() {
    let v = match load_vectors() {
        Some(v) => v,
        None => {
            eprintln!("skipping (no vectors.json)");
            return;
        }
    };
    let scheme = LoTRS::new(TEST_PARAMS);
    let pp_seed = hex::decode(&v.pp_seed).expect("pp_seed hex");
    let pp = scheme.setup(&pp_seed);
    let signing_seed = hex::decode(&v.signing_seed).expect("signing_seed hex");
    let mu = v.message.as_bytes();

    // decoded pk_table and sks from vectors
    let codec = lotrs::LoTRSCodec::new(TEST_PARAMS);
    let par = TEST_PARAMS;
    let mut pk_table: Vec<Vec<Vec<Vec<u64>>>> = vec![vec![]; v.n];
    let mut sk_table: Vec<Vec<Vec<Vec<u64>>>> = vec![vec![]; v.n];
    let mut pk_table_bytes: Vec<Vec<Vec<u8>>> = vec![vec![]; v.n];
    for kg in &v.keygen {
        let seed = hex::decode(&kg.seed).expect("seed hex");
        let pk_b = hex::decode(&kg.pk_bytes).expect("pk hex");
        let pk = codec.pk_decode(&pk_b).expect("pk_decode");
        let (sk, _) = scheme.keygen(&pp, &seed);
        let col = &mut pk_table[kg.col];
        let col_sk = &mut sk_table[kg.col];
        let col_b = &mut pk_table_bytes[kg.col];
        while col.len() <= kg.row {
            col.push(vec![]);
            col_sk.push(vec![]);
            col_b.push(vec![]);
        }
        col[kg.row] = pk;
        col_sk[kg.row] = sk;
        col_b[kg.row] = pk_b;
    }
    let sks: Vec<Vec<Vec<u64>>> = (0..par.T).map(|u| sk_table[v.ell][u].clone()).collect();

    let sig_bytes = scheme
        .sign(&pp, &sks, v.ell, mu, &pk_table, &signing_seed)
        .expect("Rust sign should succeed");

    assert!(
        scheme.verify(&pp, mu, &sig_bytes, &pk_table_bytes),
        "Rust-generated signature must verify"
    );
}

/// BENCH-level signing at `d = 128` with 16 signers takes ~1 min per sign
/// on a modern laptop in release mode.  Marked `#[ignore]` so the
/// default `cargo test` stays fast; run explicitly with:
///
///     cargo test --release --test interop -- --ignored bench_and_production_signing_round_trip
#[test]
#[ignore]
fn bench_and_production_signing_round_trip() {
    // Smoke test: BENCH / PRODUCTION parameter sets should now sign +
    // self-verify via the FACCT mask path.  We don't have a Python
    // signature to compare against (Python can't build those CDTs),
    // but we can confirm the Rust pipeline is consistent with itself.
    //
    // Intentionally use the smaller BENCH set; PRODUCTION's larger
    // sigmas make this slow — covered by mask_sampler_dispatch_matches_sigma
    // for the static property, and by the Python / Rust FACCT KAT for
    // correctness of the sampler itself.
    use lotrs::BENCH_PARAMS;
    let par = BENCH_PARAMS;
    let scheme = LoTRS::new(par);
    let pp = scheme.setup(&[0u8; 32]);

    // deterministic keygen for every signer in the ring
    let mut pk_table: Vec<Vec<Vec<Vec<u64>>>> = Vec::with_capacity(par.N());
    let mut sk_table: Vec<Vec<Vec<Vec<u64>>>> = Vec::with_capacity(par.N());
    let mut pk_table_bytes: Vec<Vec<Vec<u8>>> = Vec::with_capacity(par.N());
    let codec = lotrs::LoTRSCodec::new(par);
    for col in 0..par.N() {
        let mut col_pks = Vec::with_capacity(par.T);
        let mut col_sks = Vec::with_capacity(par.T);
        let mut col_pkb = Vec::with_capacity(par.T);
        for row in 0..par.T {
            let mut seed = [0u8; 32];
            seed[0] = col as u8;
            seed[1] = row as u8;
            let (sk, pk) = scheme.keygen(&pp, &seed);
            let pk_b = codec.pk_encode(&pk).expect("pk_encode");
            col_pks.push(pk);
            col_sks.push(sk);
            col_pkb.push(pk_b);
        }
        pk_table.push(col_pks);
        sk_table.push(col_sks);
        pk_table_bytes.push(col_pkb);
    }

    let ell = 3usize;
    let sks: Vec<Vec<Vec<u64>>> = (0..par.T).map(|u| sk_table[ell][u].clone()).collect();
    let mu = b"bench round trip";

    let sig_bytes = scheme
        .sign(&pp, &sks, ell, mu, &pk_table, &[0xaa; 32])
        .expect("BENCH sign should succeed");
    assert!(
        scheme.verify(&pp, mu, &sig_bytes, &pk_table_bytes),
        "BENCH self-round-trip failed"
    );
}

#[test]
fn verify_rejects_wrong_pk_table_shape() {
    let scheme = LoTRS::new(TEST_PARAMS);
    let pp = [0u8; 32];
    // wrong pk_table shape (too few columns) → false, no panic
    assert!(!scheme.verify(&pp, b"x", b"y", &vec![]));
    // one column with wrong number of rows → false
    assert!(!scheme.verify(&pp, b"x", b"y", &vec![vec![]; TEST_PARAMS.N()]));
}

#[test]
fn pk_decoder_rejects_malformed() {
    let codec = LoTRSCodec::new(TEST_PARAMS);
    // too short
    assert!(codec.pk_decode(&[0u8; 4]).is_err());
    // right length but filled with 0xFF → out-of-range
    let par = TEST_PARAMS;
    let poly_bytes = ((codec.dx_pk as usize * par.d) + 7) / 8;
    let total = par.k * poly_bytes;
    let bogus = vec![0xFFu8; total];
    assert!(codec.pk_decode(&bogus).is_err());
}
