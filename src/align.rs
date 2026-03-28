//! Alignment helpers: convert ncbistdaa to ASCII for display,
//! handle nucleotide complement, etc.

/// Convert Ncbistdaa sequence back to single-letter AA codes for alignment display.
pub fn ncbistdaa_to_ascii(seq: &[u8]) -> Vec<u8> {
    crate::db::sequence::decode_protein(seq)
}

/// Complement a nucleotide base (ASCII).
pub fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'N' => b'N',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        other => other,
    }
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}
