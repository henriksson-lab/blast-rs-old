//! Lookup tables for k-mer seeding.

use crate::matrix::ScoringMatrix;

/// Protein lookup table using neighboring words.
///
/// Uses NCBI-style shift-based indexing (charsize=5 bits per residue) for fast
/// rolling hash computation via shift+or+mask instead of multiply+modulo.
/// Table size: 2^(word_size * 5) = 32768 for word_size=3.
///
/// Cache-friendly packed layout:
/// - `presence`: bitfield (~4 KB for word_size=3) — fits in L1, quick rejection
/// - `table`: compact (offset, len) pairs for direct indexing
/// - `hits`: flat packed array of all query positions
pub struct ProteinLookup {
    pub word_size: usize,
    /// Number of bits per residue (5 for Ncbistdaa alphabet of 28)
    pub charsize: u32,
    /// Mask for rolling hash: (1 << (word_size * charsize)) - 1
    pub mask: u32,
    /// Presence bitfield: bit `code` is set if any query word matches word `code`.
    presence: Vec<u64>,
    /// Flat packed array of query positions for all words, concatenated.
    hits: Vec<u32>,
    /// Direct index: table[code] = (offset_into_hits, count)
    pub table: Vec<(u32, u32)>,
    capacity: usize,
}

/// Encode a protein word using shift-based encoding (NCBI style).
/// charsize bits per residue, matching ComputeTableIndex.
#[inline]
pub fn encode_protein_word(residues: &[u8]) -> u32 {
    let mut code = 0u32;
    for &r in residues {
        code = (code << 5) | (r as u32 & 0x1F);
    }
    code
}

impl ProteinLookup {
    /// Build lookup table from query (Ncbistdaa encoded).
    /// Uses charsize=5 bit encoding for shift-based rolling hash.
    pub fn build(query: &[u8], word_size: usize, matrix: &ScoringMatrix, threshold: i32) -> Self {
        let charsize = 5u32; // ilog2(28) + 1 = 5, matching NCBI
        let capacity = 1usize << (word_size as u32 * charsize);
        let mask = (capacity - 1) as u32;
        let qlen = query.len();

        if qlen < word_size {
            return ProteinLookup {
                word_size, charsize, mask,
                presence: vec![0u64; capacity.div_ceil(64)],
                hits: Vec::new(),
                table: vec![(0, 0); capacity],
                capacity,
            };
        }

        // Phase 1: collect hits using temporary Vec per code
        let mut tmp = vec![Vec::new(); capacity];
        for q_pos in 0..=(qlen - word_size) {
            let query_word = &query[q_pos..q_pos + word_size];
            // Enumerate neighbors but encode with shift-based scheme
            enumerate_neighbors_shift(query_word, word_size, charsize, matrix, threshold, &mut |code| {
                tmp[code as usize].push(q_pos as u32);
            });
        }

        // Phase 2: pack into flat arrays
        let mut presence = vec![0u64; capacity.div_ceil(64)];
        let mut total = 0u32;
        let mut offsets = vec![0u32; capacity + 1];
        for code in 0..capacity {
            offsets[code] = total;
            if !tmp[code].is_empty() {
                presence[code / 64] |= 1u64 << (code % 64);
            }
            total += tmp[code].len() as u32;
        }
        offsets[capacity] = total;

        let mut hits = vec![0u32; total as usize];
        for code in 0..capacity {
            let start = offsets[code] as usize;
            for (i, &pos) in tmp[code].iter().enumerate() {
                hits[start + i] = pos;
            }
        }

        let table: Vec<(u32, u32)> = (0..capacity)
            .map(|code| (offsets[code], offsets[code + 1] - offsets[code]))
            .collect();

        ProteinLookup { word_size, charsize, mask, presence, hits, table, capacity }
    }

    /// Get the hit slice for a word code. Returns empty slice if no hits.
    #[inline]
    pub fn get_hits(&self, code: u32) -> &[u32] {
        let code = code as usize;
        if code < self.capacity {
            let word = code / 64;
            let bit = code % 64;
            if unsafe { *self.presence.get_unchecked(word) } & (1u64 << bit) == 0 {
                return &[];
            }
            let (start, len) = unsafe { *self.table.get_unchecked(code) };
            unsafe { self.hits.get_unchecked(start as usize..(start + len) as usize) }
        } else {
            &[]
        }
    }

    /// Look up query positions for a subject word (Ncbistdaa encoded).
    #[inline]
    pub fn lookup(&self, subject_word: &[u8]) -> Option<&[u32]> {
        let code = encode_protein_word(subject_word) as usize;
        if code < self.capacity {
            let (start, len) = self.table[code];
            if len > 0 {
                Some(&self.hits[start as usize..(start + len) as usize])
            } else {
                None
            }
        } else {
            None
        }
    }
}

/// Enumerate neighbors using shift-based (charsize=5) encoding.
fn enumerate_neighbors_shift(
    query_word: &[u8],
    word_size: usize,
    charsize: u32,
    matrix: &ScoringMatrix,
    threshold: i32,
    callback: &mut impl FnMut(u32),
) {
    let mut max_suffix = vec![0i32; word_size + 1];
    for i in (0..word_size).rev() {
        let q = query_word[i];
        let best = (0u8..28).map(|r| matrix.score(q, r)).max().unwrap_or(0);
        max_suffix[i] = max_suffix[i + 1] + best;
    }
    let mut score_rows = vec![[0i32; 28]; word_size];
    for (i, row) in score_rows.iter_mut().enumerate() {
        let q = query_word[i];
        for r in 0u8..28 {
            row[r as usize] = matrix.score(q, r);
        }
    }
    enumerate_rec_shift(&score_rows, &max_suffix, word_size, charsize, threshold, 0, 0, 0, callback);
}

/// Recursive neighbor enumeration with shift-based encoding.
#[allow(clippy::too_many_arguments)]
fn enumerate_rec_shift(
    score_rows: &[[i32; 28]],
    max_suffix: &[i32],
    word_size: usize,
    charsize: u32,
    threshold: i32,
    pos: usize,
    current_score: i32,
    current_code: u32,
    callback: &mut impl FnMut(u32),
) {
    if pos == word_size {
        if current_score >= threshold {
            callback(current_code);
        }
        return;
    }
    if current_score + max_suffix[pos] < threshold { return; }
    let scores = &score_rows[pos];
    for r in 0u8..28 {
        let new_score = current_score + scores[r as usize];
        if new_score + max_suffix[pos + 1] < threshold { continue; }
        // Shift-based encoding: (code << charsize) | residue
        let new_code = (current_code << charsize) | r as u32;
        enumerate_rec_shift(score_rows, max_suffix, word_size, charsize, threshold, pos + 1, new_score, new_code, callback);
    }
}

/// Nucleotide lookup table.
/// Word size typically 11. Uses a flat array indexed by 2-bit packed word.
pub struct NucleotideLookup {
    pub word_size: usize,
    /// positions[word_code] = list of query positions
    pub table: Vec<Vec<u32>>,
    pub capacity: usize,
}

impl NucleotideLookup {
    pub fn build(query: &[u8], word_size: usize) -> Self {
        assert!(word_size <= 16, "word_size must be <= 16 for 32-bit index");
        let capacity = 1usize << (2 * word_size);
        let mut table = vec![Vec::new(); capacity];

        if query.len() < word_size {
            return NucleotideLookup { word_size, table, capacity };
        }

        // Convert query to 2-bit encoded
        let encoded: Vec<u8> = query.iter().map(|&b| crate::matrix::nt_to_2bit(b)).collect();

        // Build initial word from first (word_size-1) bases
        let mut word_code: u32 = 0;
        let mask = (capacity - 1) as u32;

        for &b in encoded.iter().take(word_size - 1) {
            if b > 3 {
                // Ambiguous base — reset
                word_code = 0;
                continue;
            }
            word_code = ((word_code << 2) | b as u32) & mask;
        }

        let mut valid = true;
        // Check if first word_size-1 bases are valid
        for &b in encoded.iter().take(word_size - 1) {
            if b > 3 { valid = false; break; }
        }

        for i in word_size - 1..query.len() {
            let b = encoded[i];
            if b > 3 {
                valid = false;
                word_code = 0;
                continue;
            }
            if !valid {
                // Try to rebuild
                let start = i + 1 - word_size;
                valid = true;
                word_code = 0;
                for &eb in encoded.iter().take(i + 1).skip(start) {
                    if eb > 3 { valid = false; break; }
                    word_code = ((word_code << 2) | eb as u32) & mask;
                }
                if !valid { continue; }
            } else {
                word_code = ((word_code << 2) | b as u32) & mask;
            }
            let q_pos = (i + 1 - word_size) as u32;
            table[word_code as usize].push(q_pos);
        }

        NucleotideLookup { word_size, table, capacity }
    }

    /// Returns query positions for a given 2-bit encoded word.
    pub fn lookup(&self, word_code: u32) -> &[u32] {
        &self.table[word_code as usize]
    }

    /// Scan subject sequence for hits against this lookup table.
    /// Returns iterator of (query_pos, subject_pos) hits.
    pub fn scan_subject<'a>(&'a self, subject: &'a [u8]) -> impl Iterator<Item = (u32, u32)> + 'a {
        let mask = (self.capacity - 1) as u32;
        let ws = self.word_size;
        let encoded: Vec<u8> = subject.iter().map(|&b| crate::matrix::nt_to_2bit(b)).collect();

        let mut hits = Vec::new();
        if encoded.len() < ws { return hits.into_iter(); }

        let mut word_code: u32 = 0;
        let mut valid = true;

        for &b in encoded.iter().take(ws - 1) {
            if b > 3 { valid = false; word_code = 0; continue; }
            word_code = ((word_code << 2) | b as u32) & mask;
        }
        // Recheck first window validity
        for &b in encoded.iter().take(ws - 1) {
            if b > 3 { valid = false; break; }
        }

        for i in ws - 1..encoded.len() {
            let b = encoded[i];
            if b > 3 {
                valid = false;
                word_code = 0;
                continue;
            }
            if !valid {
                let start = i + 1 - ws;
                valid = true;
                word_code = 0;
                for &eb in encoded.iter().take(i + 1).skip(start) {
                    if eb > 3 { valid = false; break; }
                    word_code = ((word_code << 2) | eb as u32) & mask;
                }
                if !valid { continue; }
            } else {
                word_code = ((word_code << 2) | b as u32) & mask;
            }
            let s_pos = (i + 1 - ws) as u32;
            for &q_pos in &self.table[word_code as usize] {
                hits.push((q_pos, s_pos));
            }
        }
        hits.into_iter()
    }
}

/// NCBI discontiguous megablast templates.
/// Template type 0: coding (optimized for coding sequences)
/// Template type 1: optimal (maximizes sensitivity)
/// Template type 2: two simultaneous templates
pub fn get_discontiguous_template(template_type: u8, template_length: usize) -> Vec<bool> {
    match (template_type, template_length) {
        // Coding template, length 21, 11 care positions
        (0, 21) => vec![
            true, true, true, false, true, true, false, true, false, true, false,
            false, true, true, false, true, false, true, true, true, true,
        ],
        // Optimal template, length 21, 11 care positions
        (1, 21) => vec![
            true, true, false, true, false, true, true, false, false, true, false,
            true, false, false, true, true, false, true, true, true, true,
        ],
        // Coding template, length 18, 11 care positions
        (0, 18) => vec![
            true, true, true, false, true, true, false, true, false, true,
            false, true, true, false, true, true, true, true,
        ],
        // Optimal template, length 18, 11 care positions
        (1, 18) => vec![
            true, true, false, true, false, true, true, false, true, false,
            true, false, true, true, true, true, true, true,
        ],
        // Default: contiguous (all care)
        _ => vec![true; template_length],
    }
}

/// Discontiguous megablast lookup table.
/// Uses spaced seed templates for more sensitive nucleotide searching.
pub struct DiscontiguousLookup {
    pub template_length: usize,
    pub word_size: usize, // number of care positions
    /// Template mask: true = care position, false = don't care
    pub template: Vec<bool>,
    /// table[code] = list of query positions
    pub table: Vec<Vec<u32>>,
    pub capacity: usize,
}

impl DiscontiguousLookup {
    pub fn build(query: &[u8], template_type: u8, template_length: usize) -> Self {
        let template = get_discontiguous_template(template_type, template_length);
        let word_size = template.iter().filter(|&&b| b).count();
        assert!(word_size <= 16, "word_size must be <= 16 for 32-bit index");
        let capacity = 1usize << (2 * word_size);
        let mut table = vec![Vec::new(); capacity];

        let encoded: Vec<u8> = query.iter().map(|&b| crate::matrix::nt_to_2bit(b)).collect();

        if encoded.len() < template_length {
            return DiscontiguousLookup { template_length, word_size, template, table, capacity };
        }

        for pos in 0..=(encoded.len() - template_length) {
            let mut code = 0u32;
            let mut valid = true;
            for (i, &care) in template.iter().enumerate() {
                if care {
                    let b = encoded[pos + i];
                    if b > 3 { valid = false; break; }
                    code = (code << 2) | b as u32;
                }
            }
            if valid {
                table[code as usize].push(pos as u32);
            }
        }

        DiscontiguousLookup { template_length, word_size, template, table, capacity }
    }

    pub fn scan_subject<'a>(&'a self, subject: &'a [u8]) -> Vec<(u32, u32)> {
        let encoded: Vec<u8> = subject.iter().map(|&b| crate::matrix::nt_to_2bit(b)).collect();
        let mut hits = Vec::new();

        if encoded.len() < self.template_length { return hits; }

        for pos in 0..=(encoded.len() - self.template_length) {
            let mut code = 0u32;
            let mut valid = true;
            for (i, &care) in self.template.iter().enumerate() {
                if care {
                    let b = encoded[pos + i];
                    if b > 3 { valid = false; break; }
                    code = (code << 2) | b as u32;
                }
            }
            if valid {
                for &q_pos in &self.table[code as usize] {
                    hits.push((q_pos, pos as u32));
                }
            }
        }
        hits
    }
}
