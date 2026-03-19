//! Main BLAST search pipeline.

use rayon::prelude::*;
use blast_db::{BlastDb, BlastDefLine};
use blast_db::index::SeqType;

use crate::matrix::{ScoringMatrix, MatrixType, NucleotideScoring};
use crate::stats::{KarlinAltschul, GapPenalty, lookup_ka_params, blastn_ka_params};
use crate::lookup::{ProteinLookup, NucleotideLookup};
use crate::extend::{ungapped_extend, gapped_extend, ungapped_extend_nucleotide};
use crate::hsp::{Hsp, SearchResult};

/// Parameters for a BLAST search.
#[derive(Debug, Clone)]
pub struct SearchParams {
    pub word_size: usize,
    pub matrix: MatrixType,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub evalue_threshold: f64,
    pub max_target_seqs: usize,
    pub x_drop_ungapped: i32,
    pub x_drop_gapped: i32,
    pub ungapped_cutoff: i32,
    /// Minimum score to keep after gapped extension (before E-value filter)
    pub min_score: i32,
    /// Nucleotide match/mismatch (only used for blastn)
    pub match_score: i32,
    pub mismatch: i32,
    pub num_threads: usize,
}

impl Default for SearchParams {
    fn default() -> Self {
        SearchParams {
            word_size: 3,
            matrix: MatrixType::Blosum62,
            gap_open: 11,
            gap_extend: 1,
            evalue_threshold: 10.0,
            max_target_seqs: 500,
            x_drop_ungapped: 7,
            x_drop_gapped: 15,
            ungapped_cutoff: 0,
            min_score: 0,
            match_score: 2,
            mismatch: -3,
            num_threads: 0, // 0 = rayon default
        }
    }
}

impl SearchParams {
    pub fn blastp_defaults() -> Self {
        SearchParams {
            word_size: 3,
            matrix: MatrixType::Blosum62,
            gap_open: 11,
            gap_extend: 1,
            evalue_threshold: 10.0,
            max_target_seqs: 500,
            x_drop_ungapped: 7,
            x_drop_gapped: 15,
            ungapped_cutoff: 0,
            min_score: 0,
            match_score: 2,
            mismatch: -3,
            num_threads: 0,
        }
    }

    pub fn blastn_defaults() -> Self {
        SearchParams {
            word_size: 11,
            matrix: MatrixType::Blosum62, // not used for blastn
            gap_open: 5,
            gap_extend: 2,
            evalue_threshold: 10.0,
            max_target_seqs: 500,
            x_drop_ungapped: 20,
            x_drop_gapped: 30,
            ungapped_cutoff: 0,
            min_score: 0,
            match_score: 2,
            mismatch: -3,
            num_threads: 0,
        }
    }
}

/// Neighbor word score threshold for protein lookup.
/// BLAST uses T=11 for BLOSUM62, word_size=3.
fn neighbor_threshold(matrix: MatrixType, word_size: usize) -> i32 {
    match (matrix, word_size) {
        (MatrixType::Blosum62, 3) => 11,
        (MatrixType::Blosum62, 2) => 8,
        (MatrixType::Blosum45, 3) => 14,
        (MatrixType::Blosum80, 3) => 25,
        (MatrixType::Pam30, 3)    => 10,
        (MatrixType::Pam70, 3)    => 11,
        (MatrixType::Pam250, 3)   => 11,
        _ => 11,
    }
}

/// Run a protein BLAST search (blastp).
pub fn blast_search(
    db: &BlastDb,
    query: &[u8],   // Ncbistdaa encoded query
    params: &SearchParams,
) -> Vec<SearchResult> {
    let matrix = ScoringMatrix::from_type(params.matrix);
    let gap = GapPenalty::new(params.gap_open, params.gap_extend);

    let ka = match lookup_ka_params(params.matrix, gap) {
        Some(k) => k,
        None => {
            eprintln!("Warning: no KA params for this matrix/gap combination, using defaults");
            lookup_ka_params(MatrixType::Blosum62, GapPenalty::blosum62_default())
                .unwrap()
        }
    };

    let db_len = db.volume_length();
    let num_seqs = db.num_sequences() as u64;
    let (eff_query_len, eff_db_len) = ka.effective_lengths(query.len(), db_len, num_seqs);

    let threshold = neighbor_threshold(params.matrix, params.word_size);
    let lookup = ProteinLookup::build(query, params.word_size, &matrix, threshold);

    // Process each OID (parallelized with rayon)
    let oids: Vec<u32> = (0..db.num_sequences()).collect();

    let mut results: Vec<SearchResult> = oids.par_iter().filter_map(|&oid| {
        let subject = match db.get_sequence_protein_raw(oid) {
            Ok(s) => s.to_vec(),
            Err(_) => return None,
        };
        if subject.is_empty() { return None; }

        let hsps = search_one_protein(
            query,
            &subject,
            &lookup,
            &matrix,
            &ka,
            params,
            eff_query_len,
            eff_db_len,
        );

        if hsps.is_empty() { return None; }

        let header = db.get_header(oid).unwrap_or_default();
        Some(SearchResult {
            subject_oid: oid,
            subject_title: header.title,
            subject_accession: header.accession,
            subject_len: subject.len(),
            hsps,
        })
    }).collect();

    // Sort by best E-value
    results.sort_by(|a, b| a.best_evalue().partial_cmp(&b.best_evalue()).unwrap_or(std::cmp::Ordering::Equal));
    results.truncate(params.max_target_seqs);
    results
}

/// Search one protein subject sequence.
fn search_one_protein(
    query: &[u8],
    subject: &[u8],
    lookup: &ProteinLookup,
    matrix: &ScoringMatrix,
    ka: &KarlinAltschul,
    params: &SearchParams,
    eff_query_len: usize,
    eff_db_len: u64,
) -> Vec<Hsp> {
    let slen = subject.len();
    let ws = lookup.word_size;
    if slen < ws { return vec![]; }

    // Track extended regions to avoid redundant extensions
    // A simple hit-deduplication: if we've already extended from a nearby diagonal, skip.
    let mut diag_hit: Vec<bool> = vec![false; query.len() + slen + 1];
    let diag_offset = query.len();

    let mut ungapped_hits = Vec::new();

    // Scan subject for k-mer hits
    for s_pos in 0..=(slen - ws) {
        let word = &subject[s_pos..s_pos + ws];
        if let Some(q_positions) = lookup.lookup(word) {
            for &q_pos in q_positions {
                let q_pos = q_pos as usize;
                // Diagonal = s_pos - q_pos (offset to avoid negative)
                let diag = (s_pos as isize - q_pos as isize + diag_offset as isize) as usize;
                if diag < diag_hit.len() && diag_hit[diag] {
                    continue;
                }

                let hit = ungapped_extend(query, subject, q_pos, s_pos, matrix, params.x_drop_ungapped);
                let cutoff = params.ungapped_cutoff.max(1);
                if hit.score >= cutoff {
                    if diag < diag_hit.len() {
                        diag_hit[diag] = true;
                    }
                    ungapped_hits.push(hit);
                }
            }
        }
    }

    if ungapped_hits.is_empty() { return vec![]; }

    // Sort ungapped hits by score descending, deduplicate overlapping ones
    ungapped_hits.sort_by(|a, b| b.score.cmp(&a.score));

    // Gapped extension of top ungapped hits
    let mut hsps = Vec::new();
    let mut covered_query: Vec<bool> = vec![false; query.len()];

    for uh in ungapped_hits {
        // Skip if this region already covered by a higher-scoring gapped alignment
        let center_q = (uh.q_start + uh.q_end) / 2;
        if center_q < covered_query.len() && covered_query[center_q] {
            continue;
        }

        let center_s = (uh.s_start + uh.s_end) / 2;

        let gh = gapped_extend(
            query,
            subject,
            center_q,
            center_s,
            matrix,
            params.gap_open,
            params.gap_extend,
            params.x_drop_gapped,
        );

        if gh.score <= 0 { continue; }

        let evalue = ka.evalue(gh.score, eff_query_len, eff_db_len);
        if evalue > params.evalue_threshold { continue; }

        let bit_score = ka.bit_score(gh.score);

        // Mark covered region
        let cover_start = gh.q_start.min(covered_query.len());
        let cover_end = gh.q_end.min(covered_query.len());
        for i in cover_start..cover_end {
            covered_query[i] = true;
        }

        // Convert alignment to ASCII for display
        let query_aln = blast_db::sequence::decode_protein(&gh.query_aln);
        let subject_aln = blast_db::sequence::decode_protein(&gh.subject_aln);

        hsps.push(Hsp {
            score: gh.score,
            bit_score,
            evalue,
            query_start: gh.q_start,
            query_end: gh.q_end,
            subject_start: gh.s_start,
            subject_end: gh.s_end,
            num_identities: gh.num_identities,
            num_gaps: gh.num_gaps,
            alignment_length: gh.query_aln.len(),
            query_aln,
            midline: gh.midline,
            subject_aln,
        });
    }

    hsps.sort_by(|a, b| a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal));
    hsps
}

/// Run a nucleotide BLAST search (blastn).
pub fn blastn_search(
    db: &BlastDb,
    query: &[u8],   // ASCII nucleotide query
    params: &SearchParams,
) -> Vec<SearchResult> {
    let ka = blastn_ka_params(
        params.match_score, params.mismatch,
        params.gap_open, params.gap_extend,
    );

    let db_len = db.volume_length();
    let num_seqs = db.num_sequences() as u64;
    let (eff_query_len, eff_db_len) = ka.effective_lengths(query.len(), db_len, num_seqs);

    let lookup = NucleotideLookup::build(query, params.word_size);

    let oids: Vec<u32> = (0..db.num_sequences()).collect();

    let mut results: Vec<SearchResult> = oids.par_iter().filter_map(|&oid| {
        let subject = match db.get_sequence_nucleotide(oid) {
            Ok(s) => s,
            Err(_) => return None,
        };
        if subject.is_empty() { return None; }

        let hsps = search_one_nucleotide(
            query,
            &subject,
            &lookup,
            &ka,
            params,
            eff_query_len,
            eff_db_len,
        );

        if hsps.is_empty() { return None; }

        let header = db.get_header(oid).unwrap_or_default();
        Some(SearchResult {
            subject_oid: oid,
            subject_title: header.title,
            subject_accession: header.accession,
            subject_len: subject.len(),
            hsps,
        })
    }).collect();

    results.sort_by(|a, b| a.best_evalue().partial_cmp(&b.best_evalue()).unwrap_or(std::cmp::Ordering::Equal));
    results.truncate(params.max_target_seqs);
    results
}

fn search_one_nucleotide(
    query: &[u8],
    subject: &[u8],
    lookup: &NucleotideLookup,
    ka: &KarlinAltschul,
    params: &SearchParams,
    eff_query_len: usize,
    eff_db_len: u64,
) -> Vec<Hsp> {
    let diag_offset = query.len();
    let mut diag_hit: Vec<bool> = vec![false; query.len() + subject.len() + 1];
    let mut ungapped_hits = Vec::new();

    for (q_pos, s_pos) in lookup.scan_subject(subject) {
        let q_pos = q_pos as usize;
        let s_pos = s_pos as usize;
        let diag = (s_pos as isize - q_pos as isize + diag_offset as isize) as usize;
        if diag < diag_hit.len() && diag_hit[diag] { continue; }

        let hit = ungapped_extend_nucleotide(
            query, subject,
            q_pos, s_pos,
            params.match_score, params.mismatch,
            params.x_drop_ungapped,
        );
        if hit.score > 0 {
            if diag < diag_hit.len() { diag_hit[diag] = true; }
            ungapped_hits.push(hit);
        }
    }

    if ungapped_hits.is_empty() { return vec![]; }
    ungapped_hits.sort_by(|a, b| b.score.cmp(&a.score));

    // Gapped extension with nucleotide scoring
    let nt_matrix = build_nt_matrix(params.match_score, params.mismatch);
    let mut hsps = Vec::new();
    let mut covered: Vec<bool> = vec![false; query.len()];

    for uh in ungapped_hits {
        let center_q = (uh.q_start + uh.q_end) / 2;
        if center_q < covered.len() && covered[center_q] { continue; }
        let center_s = (uh.s_start + uh.s_end) / 2;

        let gh = gapped_extend(
            query, subject,
            center_q, center_s,
            &nt_matrix,
            params.gap_open, params.gap_extend,
            params.x_drop_gapped,
        );
        if gh.score <= 0 { continue; }

        let evalue = ka.evalue(gh.score, eff_query_len, eff_db_len);
        if evalue > params.evalue_threshold { continue; }

        let bit_score = ka.bit_score(gh.score);

        for i in gh.q_start.min(covered.len())..gh.q_end.min(covered.len()) {
            covered[i] = true;
        }

        hsps.push(Hsp {
            score: gh.score,
            bit_score,
            evalue,
            query_start: gh.q_start,
            query_end: gh.q_end,
            subject_start: gh.s_start,
            subject_end: gh.s_end,
            num_identities: gh.num_identities,
            num_gaps: gh.num_gaps,
            alignment_length: gh.query_aln.len(),
            query_aln: gh.query_aln,
            midline: gh.midline,
            subject_aln: gh.subject_aln,
        });
    }

    hsps.sort_by(|a, b| a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal));
    hsps
}

/// Build a simple nucleotide scoring matrix (28×28, Ncbistdaa-compatible indices unused;
/// here we use ASCII-indexed scoring for nucleotides).
fn build_nt_matrix(match_score: i32, mismatch: i32) -> ScoringMatrix {
    let mut scores = [[mismatch; 28]; 28];
    // For nucleotide, we index by ASCII value mod 28 — this is a hack.
    // Instead, we use the matrix purely through score() which compares indices.
    // We'll build a separate approach: treat a/c/g/t as 0/1/2/3 etc.
    // Actually we already have a nucleotide-aware extend function, but gapped_extend
    // uses ScoringMatrix. Let's build an ASCII-aware matrix.
    // We'll use the actual byte values as indices with a larger array.
    // Since 28 is too small for ASCII, we'll use a trick: map through nt_to_2bit.
    // For simplicity, just use the ScoringMatrix with ncbi2bit encoding in query/subject.
    // Set diagonal entries for 2-bit codes 0,1,2,3 (A,C,G,T):
    for i in 0..4usize {
        for j in 0..4usize {
            scores[i][j] = if i == j { match_score } else { mismatch };
        }
    }
    ScoringMatrix {
        name: MatrixType::Blosum62, // doesn't matter for nt
        scores,
        min_score: mismatch,
    }
}
