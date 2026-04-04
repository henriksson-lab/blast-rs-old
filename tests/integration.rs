//! Integration tests for blast-rs.
//!
//! Each test builds a small in-memory BLAST database, runs a search, and
//! validates the results.

use tempfile::TempDir;

use blast_rs::db::index::SeqType;
use blast_rs::{
    BlastDb, BlastDbBuilder, SequenceEntry, SearchParams,
    blastp, blastn, blastx, tblastn, tblastx, psiblast, PsiblastParams,
    parse_fasta, reverse_complement, six_frame_translate,
    dc_megablast_search,
};

// ── Helpers ──────────────────────────────────────────────────────────────────

/// Build a protein database in a temporary directory, returning (TempDir, BlastDb).
/// The TempDir must be kept alive for the lifetime of the BlastDb.
fn build_protein_db(entries: Vec<SequenceEntry>) -> (TempDir, BlastDb) {
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("testdb");
    let mut builder = BlastDbBuilder::new(SeqType::Protein, "test protein db");
    for e in entries { builder.add(e); }
    builder.write(&base).unwrap();
    let db = BlastDb::open(&base).unwrap();
    (tmp, db)
}

/// Build a nucleotide database in a temporary directory.
fn build_nucleotide_db(entries: Vec<SequenceEntry>) -> (TempDir, BlastDb) {
    let tmp = TempDir::new().unwrap();
    let base = tmp.path().join("testdb");
    let mut builder = BlastDbBuilder::new(SeqType::Nucleotide, "test nt db");
    for e in entries { builder.add(e); }
    builder.write(&base).unwrap();
    let db = BlastDb::open(&base).unwrap();
    (tmp, db)
}

fn protein_entry(acc: &str, title: &str, seq: &str) -> SequenceEntry {
    SequenceEntry {
        title: title.to_string(),
        accession: acc.to_string(),
        sequence: seq.as_bytes().to_vec(),
        taxid: None,
    }
}

fn nt_entry(acc: &str, title: &str, seq: &str) -> SequenceEntry {
    SequenceEntry {
        title: title.to_string(),
        accession: acc.to_string(),
        sequence: seq.as_bytes().to_vec(),
        taxid: None,
    }
}

// ── BLASTP tests ─────────────────────────────────────────────────────────────

#[test]
fn blastp_exact_match() {
    let seq = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "exact match protein", seq),
        protein_entry("P002", "unrelated protein", "WWWWWWWWWWWWWWWWWWWWWWWWWWWWW"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);

    let results = blastp(&db, seq.as_bytes(), &params);

    assert!(!results.is_empty(), "should find at least one hit");
    let best = &results[0];
    assert_eq!(best.subject_accession, "P001");
    let hsp = &best.hsps[0];
    assert!((hsp.percent_identity() - 100.0).abs() < 0.01, "exact match should be 100% identity");
    assert_eq!(hsp.num_gaps, 0);
    assert_eq!(hsp.alignment_length, seq.len());
}

#[test]
fn blastp_no_hit_for_unrelated() {
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "all alanine", "AAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
    ]);
    let params = SearchParams::blastp()
        .evalue(1e-10)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);

    let results = blastp(&db, b"WWWWWWWWWWWWWWWWWWWWWWWWWWWW", &params);
    assert!(results.is_empty(), "unrelated sequences should not produce hits at strict evalue");
}

#[test]
fn blastp_finds_similar_sequence() {
    // Two sequences differing by a few residues
    let query   = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let subject = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let mutated = "MKFLILLFNILCLFPVLAAENHGVSMNAS"; // D→E, one mismatch
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "near identical", subject),
        protein_entry("P002", "one mismatch", mutated),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);

    let results = blastp(&db, query.as_bytes(), &params);
    assert!(results.len() >= 2, "should find both similar sequences");
    // The exact match should have best (lowest) evalue
    let exact = results.iter().find(|r| r.subject_accession == "P001").unwrap();
    let approx = results.iter().find(|r| r.subject_accession == "P002").unwrap();
    assert!(exact.best_evalue() <= approx.best_evalue());
}

#[test]
fn blastp_max_target_seqs_limits_results() {
    let query = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let entries: Vec<_> = (0..10)
        .map(|i| protein_entry(&format!("P{:03}", i), "copy", query))
        .collect();
    let (_tmp, db) = build_protein_db(entries);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .max_target_seqs(3)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);

    let results = blastp(&db, query.as_bytes(), &params);
    assert!(results.len() <= 3, "max_target_seqs=3 should limit results to at most 3, got {}", results.len());
}

#[test]
fn blastp_evalue_filter() {
    let query = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "exact", query),
    ]);
    // Very strict threshold
    let strict = SearchParams::blastp()
        .evalue(1e-100)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);
    // Loose threshold
    let loose = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);

    let strict_results = blastp(&db, query.as_bytes(), &strict);
    let loose_results = blastp(&db, query.as_bytes(), &loose);
    assert!(loose_results.len() >= strict_results.len(),
        "loose evalue should find at least as many hits as strict");
}

#[test]
fn blastp_hsp_coordinates() {
    let query = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "target", query),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);

    let results = blastp(&db, query.as_bytes(), &params);
    assert!(!results.is_empty());
    let hsp = &results[0].hsps[0];
    assert_eq!(hsp.query_start, 0);
    assert_eq!(hsp.query_end, query.len());
    assert_eq!(hsp.subject_start, 0);
    assert_eq!(hsp.subject_end, query.len());
    assert_eq!(hsp.query_aln.len(), hsp.alignment_length);
    assert_eq!(hsp.subject_aln.len(), hsp.alignment_length);
    assert_eq!(hsp.midline.len(), hsp.alignment_length);
}

#[test]
fn blastp_multiple_hsps() {
    // A subject with two separate matching regions separated by unrelated sequence
    let region1 = "MKFLILLFNILCLFPVLAAD";
    let spacer = "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW";
    let region2 = "NHGVSMNASQRDHFKLAEV";
    let subject = format!("{}{}{}", region1, spacer, region2);
    let query = format!("{}{}", region1, region2);

    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "two regions", &subject),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);

    let results = blastp(&db, query.as_bytes(), &params);
    // Should find the hit; may produce 1 or 2 HSPs depending on algorithm
    assert!(!results.is_empty(), "should find the target with matching regions");
}

// ── BLASTN tests ─────────────────────────────────────────────────────────────

#[test]
fn blastn_exact_match() {
    // Use a non-repetitive, biologically plausible sequence (>100bp to exceed word_size=11)
    let seq = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "exact nt", seq),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, seq.as_bytes(), &params);
    assert!(!results.is_empty(), "should find exact nucleotide match");
    let hsp = &results[0].hsps[0];
    assert!((hsp.percent_identity() - 100.0).abs() < 0.01);
}

#[test]
fn blastn_reverse_complement_hit() {
    let seq = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let rc = String::from_utf8(reverse_complement(seq.as_bytes())).unwrap();
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "forward strand", seq),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .strand("both");

    // Query with reverse complement should still find the hit
    let results = blastn(&db, rc.as_bytes(), &params);
    assert!(!results.is_empty(), "should find hit on reverse complement strand");
}

#[test]
fn blastn_no_hit_unrelated() {
    // Two unrelated complex sequences with no shared 11-mers
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "seq A", "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC"),
    ]);
    let params = SearchParams::blastn()
        .evalue(1e-10)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, b"TTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGC", &params);
    assert!(results.is_empty(), "unrelated nucleotide sequences should not match");
}

#[test]
fn blastn_mismatch_scoring() {
    let seq     = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let mutated = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACAATCGTAAGGCCTTAGCAGTCAATGC"; // G→A at pos 76
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "original", seq),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, mutated.as_bytes(), &params);
    assert!(!results.is_empty());
    let hsp = &results[0].hsps[0];
    assert!(hsp.percent_identity() > 90.0, "one mismatch in ~100bp should be >90% identity");
    assert!(hsp.percent_identity() < 100.0, "should not be 100% with mismatch");
}

// ── BLASTX tests ─────────────────────────────────────────────────────────────

#[test]
fn blastx_finds_translated_hit() {
    // Create a nucleotide query that encodes a known protein
    // ATG=M, AAA=K, TTT=F, CTG=L, ATT=I, CTG=L, CTG=L, TTT=F
    let nt_query = "ATGAAATTTCTGATTCTGCTGTTT";
    let protein = "MKFLILLF";

    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "target protein", protein),
    ]);
    let params = SearchParams::blastp() // blastx uses protein-style params
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);

    let results = blastx(&db, nt_query.as_bytes(), &params);
    assert!(!results.is_empty(), "blastx should find the translated protein");
    let hsp = &results[0].hsps[0];
    assert!(hsp.percent_identity() > 80.0);
    assert!(hsp.query_frame != 0, "blastx HSP should have a non-zero query frame");
}

// ── TBLASTN tests ────────────────────────────────────────────────────────────

#[test]
fn tblastn_finds_protein_in_nt_db() {
    // Protein query against a nucleotide database (translated subject)
    let protein_query = "MKFLILLF";
    // The nucleotide that encodes this protein (frame +1)
    let nt_subject = "ATGAAATTTCTGATTCTGCTGTTT";

    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "coding region", nt_subject),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);

    let results = tblastn(&db, protein_query.as_bytes(), &params);
    assert!(!results.is_empty(), "tblastn should find protein in translated nucleotide db");
    let hsp = &results[0].hsps[0];
    assert!(hsp.subject_frame != 0, "tblastn HSP should have non-zero subject frame");
}

// ── TBLASTX tests ────────────────────────────────────────────────────────────

#[test]
fn tblastx_translated_vs_translated() {
    // Both query and subject are nucleotide, both get translated
    let nt_seq = "ATGAAATTTCTGATTCTGCTGTTTAACATTCTGTGCCTGTTC";
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "coding nt", nt_seq),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);

    let results = tblastx(&db, nt_seq.as_bytes(), &params);
    assert!(!results.is_empty(), "tblastx should find self-hit in translated mode");
    let hsp = &results[0].hsps[0];
    assert!(hsp.query_frame != 0, "tblastx should set query frame");
    assert!(hsp.subject_frame != 0, "tblastx should set subject frame");
}

// ── PSI-BLAST test ───────────────────────────────────────────────────────────

#[test]
fn psiblast_iterates_and_returns_pssm() {
    let query = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let entries = vec![
        protein_entry("P001", "exact", query),
        protein_entry("P002", "similar", "MKFLILLFNILCLFPVLAAENHGVSMNAS"),
        protein_entry("P003", "unrelated", "WWWWWWWWWWWWWWWWWWWWWWWWWWWWW"),
    ];
    let (_tmp, db) = build_protein_db(entries);
    let search_params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);
    let psi_params = PsiblastParams::new(search_params).num_iterations(2);

    let (results, pssm) = psiblast(&db, query.as_bytes(), &psi_params);
    assert!(!results.is_empty(), "psiblast should return results");
    // PSSM should have dimensions matching the query
    assert_eq!(pssm.scores.len(), query.len(),
        "PSSM should have one row per query position");
}

// ── Database round-trip tests ────────────────────────────────────────────────

#[test]
fn protein_db_roundtrip() {
    let seq = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "test seq", seq),
    ]);
    assert_eq!(db.num_sequences(), 1);
    assert_eq!(db.seq_type(), SeqType::Protein);
    let decoded = db.get_sequence_protein(0).unwrap();
    assert_eq!(String::from_utf8(decoded).unwrap(), seq);
}

#[test]
fn nucleotide_db_roundtrip() {
    let seq = "ATGGCTAGCGATCGATCGATCGATCG";
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "test nt", seq),
    ]);
    assert_eq!(db.num_sequences(), 1);
    assert_eq!(db.seq_type(), SeqType::Nucleotide);
    let decoded = db.get_sequence_nucleotide(0).unwrap();
    assert_eq!(String::from_utf8(decoded).unwrap(), seq);
}

#[test]
fn db_multiple_sequences() {
    let entries = vec![
        protein_entry("P001", "first", "MKFLILLFNILCLFPVLAAD"),
        protein_entry("P002", "second", "NHGVSMNASQRDHFKLAEV"),
        protein_entry("P003", "third", "ACDEFGHIKLMNPQRSTVWY"),
    ];
    let (_tmp, db) = build_protein_db(entries);
    assert_eq!(db.num_sequences(), 3);

    // Verify each sequence can be retrieved
    for oid in 0..3 {
        let seq = db.get_sequence_protein(oid).unwrap();
        assert!(!seq.is_empty());
        let deflines = db.get_headers(oid).unwrap();
        assert!(!deflines.is_empty());
    }
}

#[test]
fn db_header_roundtrip() {
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("ACC123", "my protein title", "ACDEFGHIKLMNPQRSTVWY"),
    ]);
    let deflines = db.get_headers(0).unwrap();
    assert_eq!(deflines.len(), 1);
    assert_eq!(deflines[0].accession, "ACC123");
    assert!(deflines[0].title.contains("my protein title"));
}

// ── Reverse complement test ─────────────────────────────────────────────────

#[test]
fn reverse_complement_roundtrip() {
    let seq = b"ATGGCTAGCGATCG";
    let rc = reverse_complement(seq);
    let rc2 = reverse_complement(&rc);
    assert_eq!(&rc2, seq, "double reverse complement should return original");
}

// ── Six-frame translation test ──────────────────────────────────────────────

#[test]
fn six_frame_translate_produces_six_frames() {
    let seq = b"ATGGCTAGCGATCGATCGATCGATCG";
    let frames = six_frame_translate(seq);
    assert_eq!(frames.len(), 6);
    // Frames +1, +2, +3, -1, -2, -3
    let frame_nums: Vec<i32> = frames.iter().map(|f| f.frame).collect();
    assert_eq!(frame_nums, vec![1, 2, 3, -1, -2, -3]);
    // Frame +1 starts with M (ATG = Met)
    assert_eq!(frames[0].protein[0], b'M');
}

// ── FASTA parsing edge cases ─────────────────────────────────────────────────

#[test]
fn parse_fasta_with_blank_lines() {
    let input = b">seq1\n\nACGT\n\nTGCA\n\n>seq2\nAAAA\n";
    let seqs = parse_fasta(input);
    assert_eq!(seqs.len(), 2);
    assert_eq!(&seqs[0].1, b"ACGTTGCA".as_slice());
    assert_eq!(&seqs[1].1, b"AAAA".as_slice());
}

#[test]
fn parse_fasta_no_trailing_newline() {
    let input = b">seq1\nACGT";
    let seqs = parse_fasta(input);
    assert_eq!(seqs.len(), 1);
    assert_eq!(&seqs[0].1, b"ACGT".as_slice());
}

// ── Multi-threaded search ────────────────────────────────────────────────────

#[test]
fn blastp_multithreaded_same_results() {
    let query = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let entries: Vec<_> = (0..20)
        .map(|i| protein_entry(&format!("P{:03}", i), "copy", query))
        .collect();
    let (_tmp, db) = build_protein_db(entries);

    let single = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);
    let multi = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(4)
        .filter_low_complexity(false)
        .comp_adjust(false);

    let r1 = blastp(&db, query.as_bytes(), &single);
    let r2 = blastp(&db, query.as_bytes(), &multi);
    assert_eq!(r1.len(), r2.len(),
        "single-threaded and multi-threaded should return same number of results");
}

// ── Scoring matrix tests ────────────────────────────────────────────────────

#[test]
fn blastp_different_matrices() {
    let query = "MKFLILLFNILCLFPVLAADNHGVSMNAS";
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "target", query),
    ]);

    for matrix in &[
        blast_rs::MatrixType::Blosum62,
        blast_rs::MatrixType::Blosum45,
        blast_rs::MatrixType::Blosum80,
    ] {
        let params = SearchParams::blastp()
            .evalue(10.0)
            .matrix(*matrix)
            .num_threads(1)
            .filter_low_complexity(false)
            .comp_adjust(false);

        let results = blastp(&db, query.as_bytes(), &params);
        assert!(!results.is_empty(),
            "should find exact match with {:?} matrix", matrix);
    }
}

// ── Discontiguous megablast tests ───────────────────────────────────────────

#[test]
fn dc_megablast_exact_match() {
    let seq = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "exact nt", seq),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = dc_megablast_search(&db, seq.as_bytes(), &params, 0, 18);
    assert!(!results.is_empty(), "dc_megablast should find exact match");
    let hsp = &results[0].hsps[0];
    assert!((hsp.percent_identity() - 100.0).abs() < 0.01);
    // Verify alignment strings contain proper ASCII nucleotides
    assert!(hsp.query_aln.iter().all(|&b| b == b'-' || b"ACGTN".contains(&b)),
        "alignment should contain ASCII nucleotide letters, got: {:?}",
        String::from_utf8_lossy(&hsp.query_aln));
}

#[test]
fn dc_megablast_mismatch() {
    let seq     = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let mutated = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACAATCGTAAGGCCTTAGCAGTCAATGC";
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "original", seq),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = dc_megablast_search(&db, mutated.as_bytes(), &params, 0, 18);
    assert!(!results.is_empty(), "dc_megablast should find near-match");
    let hsp = &results[0].hsps[0];
    assert!(hsp.percent_identity() > 90.0);
    assert!(hsp.percent_identity() < 100.0);
}

// ── Edge case tests ─────────────────────────────────────────────────────────

#[test]
fn blastp_short_query_below_word_size() {
    // Query shorter than default word_size=3 should not panic
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "target", "MKFLILLFNILCLFPVLAADNHGVSMNAS"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);

    // 2-residue query, below word_size=3
    let results = blastp(&db, b"MK", &params);
    // May or may not find hits, but should not panic
    let _ = results;
}

#[test]
fn blastp_single_residue_query() {
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "target", "MKFLILLFNILCLFPVLAADNHGVSMNAS"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);

    let results = blastp(&db, b"M", &params);
    let _ = results; // should not panic
}

#[test]
fn blastn_short_query_below_word_size() {
    // Query shorter than word_size=11
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "target", "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC"),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, b"ATGCGT", &params);
    let _ = results; // should not panic
}

#[test]
fn blastp_empty_query() {
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "target", "MKFLILLFNILCLFPVLAADNHGVSMNAS"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);

    let results = blastp(&db, b"", &params);
    assert!(results.is_empty(), "empty query should produce no results");
}

#[test]
fn blastn_empty_query() {
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "target", "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC"),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, b"", &params);
    assert!(results.is_empty(), "empty query should produce no results");
}

#[test]
fn blastp_query_with_ambiguous_residues() {
    // X is an ambiguous amino acid
    let query = "MKFLXLLFNILCLFPVLAADNHGVSMNAS";
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "target", "MKFLILLFNILCLFPVLAADNHGVSMNAS"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);

    let results = blastp(&db, query.as_bytes(), &params);
    // Should still find a hit despite one X
    assert!(!results.is_empty(), "query with X should still find hits");
}

#[test]
fn blastn_query_with_ambiguous_bases() {
    // N is an ambiguous nucleotide
    let query   = "ATGCGTACCTGAAAGCTTCAGNACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let subject = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "target", subject),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    // Should not panic; may or may not find hit depending on N handling
    let results = blastn(&db, query.as_bytes(), &params);
    let _ = results;
}

#[test]
fn blastp_short_subject_in_db() {
    // Subject shorter than word_size
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "tiny", "MK"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);

    let results = blastp(&db, b"MKFLILLFNILCLFPVLAADNHGVSMNAS", &params);
    let _ = results; // should not panic
}

#[test]
fn blastn_alignment_strings_are_ascii() {
    let seq = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "target", seq),
    ]);
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false);

    let results = blastn(&db, seq.as_bytes(), &params);
    assert!(!results.is_empty());
    let hsp = &results[0].hsps[0];
    // Verify alignment strings contain proper ASCII letters, not 2-bit codes
    for &b in &hsp.query_aln {
        assert!(b == b'-' || b.is_ascii_alphabetic(),
            "query_aln byte {} is not ASCII letter or gap", b);
    }
    for &b in &hsp.subject_aln {
        assert!(b == b'-' || b.is_ascii_alphabetic(),
            "subject_aln byte {} is not ASCII letter or gap", b);
    }
}

#[test]
fn blastx_empty_nt_query() {
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "target", "MKFLILLFNILCLFPVLAAD"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);

    let results = blastx(&db, b"", &params);
    assert!(results.is_empty());
}

#[test]
fn tblastn_empty_protein_query() {
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "coding", "ATGAAATTTCTGATTCTGCTGTTT"),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);

    let results = tblastn(&db, b"", &params);
    assert!(results.is_empty());
}

#[test]
fn blastn_plus_strand_only() {
    let seq = "ATGCGTACCTGAAAGCTTCAGTACGGTAATCCTGAACGTTAGCCAATGCTTGAAGTCAACGTATCGCAAGCTTAACGATCGTAAGGCCTTAGCAGTCAATGC";
    let rc = String::from_utf8(reverse_complement(seq.as_bytes())).unwrap();
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "forward only", seq),
    ]);

    // With strand=plus, querying the reverse complement should find no hit
    let params = SearchParams::blastn()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .strand("plus");

    let results = blastn(&db, rc.as_bytes(), &params);
    assert!(results.is_empty(), "plus-strand-only should not find reverse complement hit");
}

#[test]
fn blastp_all_twenty_amino_acids() {
    // Ensure all standard amino acids are handled correctly
    let query = "ACDEFGHIKLMNPQRSTVWY";
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "all aa", query),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);

    let results = blastp(&db, query.as_bytes(), &params);
    assert!(!results.is_empty());
    let hsp = &results[0].hsps[0];
    assert!((hsp.percent_identity() - 100.0).abs() < 0.01);
    assert_eq!(hsp.alignment_length, 20);
}

#[test]
fn blastp_stop_codon_in_sequence() {
    // * represents stop codon / selenocysteine in protein
    let query   = "MKFLILLFNILCLFPVLAAD";
    let subject = "MKFLILLF*ILCLFPVLAAD"; // stop codon in middle
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "with stop", subject),
    ]);
    let params = SearchParams::blastp()
        .evalue(10.0)
        .num_threads(1)
        .filter_low_complexity(false)
        .comp_adjust(false);

    // Should not panic
    let results = blastp(&db, query.as_bytes(), &params);
    let _ = results;
}

#[test]
fn db_single_base_nucleotide() {
    // Minimal nucleotide sequence
    let (_tmp, db) = build_nucleotide_db(vec![
        nt_entry("N001", "single", "A"),
    ]);
    assert_eq!(db.num_sequences(), 1);
    let decoded = db.get_sequence_nucleotide(0).unwrap();
    assert_eq!(String::from_utf8(decoded).unwrap(), "A");
}

#[test]
fn db_single_residue_protein() {
    let (_tmp, db) = build_protein_db(vec![
        protein_entry("P001", "single", "M"),
    ]);
    assert_eq!(db.num_sequences(), 1);
    let decoded = db.get_sequence_protein(0).unwrap();
    assert_eq!(String::from_utf8(decoded).unwrap(), "M");
}
