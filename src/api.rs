//! High-level public API for running BLAST searches.
//!
//! # Example
//! ```no_run
//! use blast_rs::{BlastDb, SearchParams, blastp, parse_fasta};
//!
//! let db = BlastDb::open("nr".as_ref()).unwrap();
//! let fasta = std::fs::read("query.faa").unwrap();
//! let sequences = parse_fasta(&fasta);
//!
//! let params = SearchParams::blastp()
//!     .evalue(1e-5)
//!     .max_target_seqs(50)
//!     .num_threads(8);
//!
//! for (title, seq) in &sequences {
//!     let results = blastp(&db, seq, &params);
//!     for r in &results {
//!         println!("{}\t{}\t{:.2e}", title, r.subject_accession, r.best_evalue());
//!     }
//! }
//! ```

use crate::db::BlastDb;
use crate::{SearchResult, SearchParams, Pssm};
use crate::search::{blast_search, blastn_search, blastx_search, tblastn_search, tblastx_search, aa_to_ncbistdaa};
use crate::pssm::psiblast_search;

/// PSI-BLAST-specific parameters.
#[derive(Debug, Clone)]
pub struct PsiblastParams {
    pub search: SearchParams,
    /// Number of PSI-BLAST iterations (default: 3).
    pub num_iterations: u32,
    /// E-value threshold for including hits in PSSM construction (default: 0.001).
    pub inclusion_evalue: f64,
}

impl PsiblastParams {
    pub fn new(search: SearchParams) -> Self {
        PsiblastParams {
            search,
            num_iterations: 3,
            inclusion_evalue: 0.001,
        }
    }

    pub fn num_iterations(mut self, v: u32) -> Self { self.num_iterations = v; self }
    pub fn inclusion_evalue(mut self, v: f64) -> Self { self.inclusion_evalue = v; self }
}

// ── Search functions ──────────────────────────────────────────────────────────

/// Protein query vs protein database (blastp).
///
/// `query` is an ASCII amino-acid sequence (as read from a FASTA file).
/// Encoding to Ncbistdaa is handled internally.
pub fn blastp(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    let query_ncbi = aa_to_ncbistdaa(query);
    blast_search(db, &query_ncbi, params)
}

/// Nucleotide query vs nucleotide database (blastn).
///
/// `query` is an ASCII nucleotide sequence.
pub fn blastn(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    blastn_search(db, query, params)
}

/// Translated nucleotide query vs protein database, 6 frames (blastx).
///
/// `query` is an ASCII nucleotide sequence.
pub fn blastx(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    blastx_search(db, query, params)
}

/// Protein query vs translated nucleotide database, 6 frames (tblastn).
///
/// `query` is an ASCII amino-acid sequence.
pub fn tblastn(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    tblastn_search(db, query, params)
}

/// Translated nucleotide query vs translated nucleotide database (tblastx).
///
/// `query` is an ASCII nucleotide sequence.
pub fn tblastx(db: &BlastDb, query: &[u8], params: &SearchParams) -> Vec<SearchResult> {
    tblastx_search(db, query, params)
}

/// Iterative protein search with PSSM (psiblast).
///
/// `query` is an ASCII amino-acid sequence.
/// Returns final-round hits and the resulting PSSM.
pub fn psiblast(db: &BlastDb, query: &[u8], params: &PsiblastParams) -> (Vec<SearchResult>, Pssm) {
    let query_ncbi = aa_to_ncbistdaa(query);
    psiblast_search(db, &query_ncbi, &params.search, params.num_iterations, params.inclusion_evalue)
}

// ── FASTA parser ──────────────────────────────────────────────────────────────

/// Parse a multi-FASTA byte slice into `(title, sequence)` pairs.
///
/// The returned sequences are bare ASCII bytes (no `>` header), suitable for
/// passing directly to [`blastp`], [`blastn`], etc.
pub fn parse_fasta(input: &[u8]) -> Vec<(String, Vec<u8>)> {
    let mut sequences = Vec::new();
    let mut current_title = String::new();
    let mut current_seq: Vec<u8> = Vec::new();

    for line in input.split(|&b| b == b'\n') {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.is_empty() { continue; }
        if line.starts_with(b">") {
            if !current_title.is_empty() {
                sequences.push((current_title.clone(), current_seq.clone()));
                current_seq.clear();
            }
            current_title = String::from_utf8_lossy(&line[1..]).trim().to_string();
        } else {
            current_seq.extend_from_slice(line);
        }
    }
    if !current_title.is_empty() {
        sequences.push((current_title, current_seq));
    }
    sequences
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_fasta_single() {
        let input = b">seq1 a protein\nMKTLLLTLVV\nVTIYCLLL\n";
        let seqs = parse_fasta(input);
        assert_eq!(seqs.len(), 1);
        assert_eq!(seqs[0].0, "seq1 a protein");
        assert_eq!(&seqs[0].1, b"MKTLLLTLVVVTIYCLLL".as_slice());
    }

    #[test]
    fn test_parse_fasta_multi() {
        let input = b">alpha\nACGT\n>beta\nTGCA\nACGT\n";
        let seqs = parse_fasta(input);
        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].0, "alpha");
        assert_eq!(&seqs[0].1, b"ACGT".as_slice());
        assert_eq!(seqs[1].0, "beta");
        assert_eq!(&seqs[1].1, b"TGCAACGT".as_slice());
    }

    #[test]
    fn test_parse_fasta_crlf() {
        let input = b">seq1\r\nACGT\r\n";
        let seqs = parse_fasta(input);
        assert_eq!(seqs.len(), 1);
        assert_eq!(&seqs[0].1, b"ACGT".as_slice());
    }

    #[test]
    fn test_parse_fasta_empty() {
        let seqs = parse_fasta(b"");
        assert!(seqs.is_empty());
    }
}
