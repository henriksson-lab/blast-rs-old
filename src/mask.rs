//! Low-complexity region masking: DUST (nucleotide) and SEG (protein).
//!
//! Both functions return a `Vec<bool>` mask (true = masked) and provide
//! helper functions to apply the mask in-place.

use crate::translate::nt_to_2bit;

// ─── DUST ────────────────────────────────────────────────────────────────────

/// Default DUST window size.
pub const DUST_WINDOW: usize = 64;
/// Default DUST score threshold.
pub const DUST_THRESHOLD: f64 = 2.5;

/// Compute the DUST low-complexity score for a window of nucleotide sequence.
///
/// Score = Σ c_k*(c_k-1)/2  /  (W-2)
/// where c_k is the count of each distinct 3-mer in the window.
/// Returns `f64::NAN` if the window is too short.
fn dust_score_window(bits: &[u8]) -> f64 {
    let n = bits.len();
    if n < 3 { return f64::NAN; }
    let mut counts = [0u32; 64];
    for i in 0..n - 2 {
        let (b1, b2, b3) = (bits[i], bits[i+1], bits[i+2]);
        if b1 < 4 && b2 < 4 && b3 < 4 {
            counts[(b1 as usize) * 16 + (b2 as usize) * 4 + b3 as usize] += 1;
        }
    }
    let sum: f64 = counts.iter().map(|&c| (c as f64) * (c as f64 - 1.0) / 2.0).sum();
    sum / (n - 2) as f64
}

/// Compute a DUST mask for an ASCII nucleotide sequence.
///
/// Positions inside windows with score ≥ `threshold` are marked as masked.
/// Uses a sliding window of `window` bases (default: 64).
pub fn dust_mask(seq: &[u8], window: usize, threshold: f64) -> Vec<bool> {
    let n = seq.len();
    let mut masked = vec![false; n];
    if n < 3 { return masked; }

    // Pre-convert to 2-bit
    let bits: Vec<u8> = seq.iter().map(|&c| nt_to_2bit(c)).collect();

    for start in 0..n {
        let end = (start + window).min(n);
        if end - start < 3 { break; }
        let score = dust_score_window(&bits[start..end]);
        if !score.is_nan() && score >= threshold {
            for i in start..end {
                masked[i] = true;
            }
        }
    }
    masked
}

/// Apply DUST masking in-place, replacing masked positions with `b'N'`.
pub fn apply_dust(seq: &mut Vec<u8>) {
    apply_dust_opts(seq, DUST_WINDOW, DUST_THRESHOLD);
}

/// Apply DUST masking with custom parameters.
pub fn apply_dust_opts(seq: &mut Vec<u8>, window: usize, threshold: f64) {
    let mask = dust_mask(seq, window, threshold);
    for (i, &m) in mask.iter().enumerate() {
        if m { seq[i] = b'N'; }
    }
}

// ─── SEG ─────────────────────────────────────────────────────────────────────

/// Default SEG window size.
pub const SEG_WINDOW: usize = 12;
/// Default SEG low-entropy trigger threshold (bits).
pub const SEG_K1: f64 = 2.2;
/// Default SEG extension threshold (bits).
pub const SEG_K2: f64 = 2.5;

/// Compute the Shannon entropy (in bits) of the amino acid composition
/// in `window`. Uses all ASCII letters (case-insensitive, 26 buckets).
fn seg_entropy(window: &[u8]) -> f64 {
    let n = window.len();
    if n == 0 { return 0.0; }
    let mut counts = [0u32; 26];
    for &b in window {
        let idx = b.to_ascii_uppercase().wrapping_sub(b'A') as usize;
        if idx < 26 { counts[idx] += 1; }
    }
    let nf = n as f64;
    counts.iter().fold(0.0f64, |h, &c| {
        if c == 0 { h } else {
            let p = c as f64 / nf;
            h - p * p.log2()
        }
    })
}

/// Compute a SEG mask for an ASCII amino acid sequence.
///
/// Windows with Shannon entropy < `k1` are marked; the mask is then extended
/// as long as the extended window entropy stays below `k2`.
pub fn seg_mask(seq: &[u8], window: usize, k1: f64, k2: f64) -> Vec<bool> {
    let n = seq.len();
    let mut masked = vec![false; n];
    if n < window { return masked; }

    // Phase 1: identify trigger windows
    let mut triggers: Vec<bool> = vec![false; n];
    for start in 0..=(n - window) {
        if seg_entropy(&seq[start..start + window]) < k1 {
            for i in start..start + window {
                triggers[i] = true;
            }
        }
    }

    // Phase 2: extend triggered regions while entropy < k2
    let mut i = 0;
    while i < n {
        if triggers[i] {
            // Find the end of this triggered block
            let mut seg_start = i;
            let mut seg_end = i;
            while seg_end < n && triggers[seg_end] { seg_end += 1; }

            // Trim: shrink from both ends while extending would exceed k2
            // Extend left
            while seg_start > 0 {
                let new_start = seg_start - 1;
                let w = (seg_end - new_start).min(window * 2);
                if seg_entropy(&seq[new_start..new_start + w]) < k2 {
                    seg_start = new_start;
                } else {
                    break;
                }
            }
            // Extend right
            while seg_end < n {
                let w = (seg_end + 1 - seg_start).min(window * 2);
                let rstart = seg_end + 1 - w;
                if seg_entropy(&seq[rstart..seg_end + 1]) < k2 {
                    seg_end += 1;
                } else {
                    break;
                }
            }

            for j in seg_start..seg_end {
                masked[j] = true;
            }
            i = seg_end;
        } else {
            i += 1;
        }
    }
    masked
}

/// Apply SEG masking in-place, replacing masked positions with ASCII `b'X'`.
pub fn apply_seg(seq: &mut Vec<u8>) {
    apply_seg_opts(seq, SEG_WINDOW, SEG_K1, SEG_K2);
}

/// Apply SEG masking with custom parameters.
pub fn apply_seg_opts(seq: &mut Vec<u8>, window: usize, k1: f64, k2: f64) {
    let mask = seg_mask(seq, window, k1, k2);
    for (i, &m) in mask.iter().enumerate() {
        if m { seq[i] = b'X'; }
    }
}

/// Apply SEG masking to a Ncbistdaa-encoded protein sequence (replaces with code 21 = X).
pub fn apply_seg_ncbistdaa(seq: &mut Vec<u8>) {
    // Decode to ASCII, mask, then find which positions were masked
    let ascii: Vec<u8> = seq.iter().map(|&c| {
        crate::db::sequence::decode_protein(&[c]).into_iter().next().unwrap_or(b'X')
    }).collect();
    let mask = seg_mask(&ascii, SEG_WINDOW, SEG_K1, SEG_K2);
    for (i, &m) in mask.iter().enumerate() {
        if m { seq[i] = 21; } // Ncbistdaa X
    }
}

// ─── Repeat Masker (simplified WindowMasker-like) ───────────────────────────

/// Default repeat window size.
pub const REPEAT_WINDOW: usize = 256;
/// Default repeat n-mer size for counting.
pub const REPEAT_NMER: usize = 15;
/// Default repeat frequency threshold (mask n-mers occurring more than this fraction).
pub const REPEAT_THRESHOLD: f64 = 0.01;

/// Compute a repeat mask for a nucleotide sequence using n-mer frequency analysis.
///
/// This is a simplified version of NCBI's WindowMasker. It identifies
/// positions covered by overrepresented n-mers (repeats/low-complexity regions).
///
/// Algorithm:
/// 1. Count all n-mers in the sequence
/// 2. Positions covered by n-mers occurring > threshold * total are masked
pub fn repeat_mask(seq: &[u8], nmer_size: usize, threshold: f64) -> Vec<bool> {
    use std::collections::HashMap;
    let n = seq.len();
    let mut masked = vec![false; n];
    if n < nmer_size { return masked; }

    let bits: Vec<u8> = seq.iter().map(|&c| nt_to_2bit(c)).collect();

    // Count n-mers
    let mut counts: HashMap<u64, u32> = HashMap::new();
    let mut total = 0u32;

    for i in 0..=(n - nmer_size) {
        let mut code = 0u64;
        let mut valid = true;
        for j in 0..nmer_size {
            let b = bits[i + j];
            if b > 3 { valid = false; break; }
            code = (code << 2) | b as u64;
        }
        if valid {
            *counts.entry(code).or_insert(0) += 1;
            total += 1;
        }
    }

    if total == 0 { return masked; }

    // Determine frequency cutoff
    let cutoff = (threshold * total as f64).max(2.0) as u32;

    // Mark positions covered by overrepresented n-mers
    for i in 0..=(n - nmer_size) {
        let mut code = 0u64;
        let mut valid = true;
        for j in 0..nmer_size {
            let b = bits[i + j];
            if b > 3 { valid = false; break; }
            code = (code << 2) | b as u64;
        }
        if valid {
            if let Some(&count) = counts.get(&code) {
                if count >= cutoff {
                    for j in i..i + nmer_size {
                        masked[j] = true;
                    }
                }
            }
        }
    }
    masked
}

/// Apply repeat masking in-place, replacing masked positions with `b'N'`.
pub fn apply_repeat_mask(seq: &mut Vec<u8>) {
    apply_repeat_mask_opts(seq, REPEAT_NMER, REPEAT_THRESHOLD);
}

/// Apply repeat masking with custom parameters.
pub fn apply_repeat_mask_opts(seq: &mut Vec<u8>, nmer_size: usize, threshold: f64) {
    let mask = repeat_mask(seq, nmer_size, threshold);
    for (i, &m) in mask.iter().enumerate() {
        if m { seq[i] = b'N'; }
    }
}

// ─── Lowercase masking helper ───────────────────────────────────────────────

/// Create a mask from lowercase characters in a sequence.
/// Returns true for positions that were lowercase in the original.
pub fn lowercase_mask(seq: &[u8]) -> Vec<bool> {
    seq.iter().map(|&b| b.is_ascii_lowercase()).collect()
}

/// Apply lowercase masking: replace lowercase positions with X (protein) or N (nucleotide).
pub fn apply_lowercase_mask_protein(seq: &mut Vec<u8>) {
    for b in seq.iter_mut() {
        if b.is_ascii_lowercase() {
            *b = b'X';
        }
    }
}

pub fn apply_lowercase_mask_nucleotide(seq: &mut Vec<u8>) {
    for b in seq.iter_mut() {
        if b.is_ascii_lowercase() {
            *b = b'N';
        }
    }
}
