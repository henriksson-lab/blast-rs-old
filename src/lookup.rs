//! Lookup tables for k-mer seeding.

use crate::matrix::ScoringMatrix;

/// Protein lookup table using neighboring words.
///
/// Cache-friendly packed layout:
/// - `presence`: bitfield (2.7 KB for word_size=3) — fits in L1, quick rejection
/// - `offsets`: compact offset table (86 KB) — fits in L2
/// - `hits`: flat packed array of all query positions
///
/// For word_size=3 the naive Vec<Vec<u32>> layout uses 514 KB (doesn't fit in L2),
/// causing cache misses on every subject word lookup. This packed layout reduces
/// working set by ~6x.
pub struct ProteinLookup {
    pub word_size: usize,
    /// Presence bitfield: bit `code` is set if any query word matches word `code`.
    presence: Vec<u64>,
    /// offsets[code] = start index into `hits` for word `code`.
    /// offsets[capacity] = total length of `hits` (sentinel). Used during build.
    #[allow(dead_code)]
    offsets: Vec<u32>,
    /// Flat packed array of query positions for all words, concatenated.
    hits: Vec<u32>,
    /// For direct indexing from search loop: table[code] = slice into hits.
    /// This is a Vec of (start, len) pairs, 8 bytes each = 171 KB for word_size=3.
    /// Still better than 514 KB for Vec<Vec<u32>>.
    pub table: Vec<(u32, u32)>, // (offset, length)
    capacity: usize,
}

/// Encode a protein word of length `word_size` from Ncbistdaa residues.
/// Uses base-28 encoding.
#[inline]
pub fn encode_protein_word(residues: &[u8]) -> u32 {
    let mut code = 0u32;
    for &r in residues {
        code = code * 28 + (r as u32 % 28);
    }
    code
}

impl ProteinLookup {
    /// Build lookup table from query (Ncbistdaa encoded).
    /// Includes neighboring words scoring >= threshold.
    pub fn build(query: &[u8], word_size: usize, matrix: &ScoringMatrix, threshold: i32) -> Self {
        let capacity = 28usize.pow(word_size as u32);
        let qlen = query.len();

        if qlen < word_size {
            return ProteinLookup {
                word_size,
                presence: vec![0u64; capacity.div_ceil(64)],
                offsets: vec![0u32; capacity + 1],
                hits: Vec::new(),
                table: vec![(0, 0); capacity],
                capacity,
            };
        }

        // Phase 1: collect hits into temporary Vec<Vec<u32>> (build phase only)
        let mut tmp = vec![Vec::new(); capacity];
        for q_pos in 0..=(qlen - word_size) {
            let query_word = &query[q_pos..q_pos + word_size];
            enumerate_neighbors(query_word, word_size, matrix, threshold, &mut |neighbor_code| {
                tmp[neighbor_code as usize].push(q_pos as u32);
            });
        }

        // Phase 2: pack into flat arrays
        let mut presence = vec![0u64; capacity.div_ceil(64)];
        let mut offsets = vec![0u32; capacity + 1];
        let mut total = 0u32;
        for code in 0..capacity {
            offsets[code] = total;
            let len = tmp[code].len() as u32;
            if len > 0 {
                presence[code / 64] |= 1u64 << (code % 64);
            }
            total += len;
        }
        offsets[capacity] = total;

        let mut hits = vec![0u32; total as usize];
        for code in 0..capacity {
            let start = offsets[code] as usize;
            for (i, &pos) in tmp[code].iter().enumerate() {
                hits[start + i] = pos;
            }
        }

        // Build table for direct indexing
        let table: Vec<(u32, u32)> = (0..capacity)
            .map(|code| {
                let start = offsets[code];
                let len = offsets[code + 1] - start;
                (start, len)
            })
            .collect();

        ProteinLookup { word_size, presence, offsets, hits, table, capacity }
    }

    /// Get the hit slice for a word code. Returns empty slice if no hits.
    #[inline]
    pub fn get_hits(&self, code: u32) -> &[u32] {
        let code = code as usize;
        if code < self.capacity {
            // Quick presence check (L1 cache — 2.7 KB bitfield)
            let word = code / 64;
            let bit = code % 64;
            if self.presence[word] & (1u64 << bit) == 0 {
                return &[];
            }
            let (start, len) = self.table[code];
            &self.hits[start as usize..(start + len) as usize]
        } else {
            &[]
        }
    }

    /// Look up query positions for a subject word (Ncbistdaa encoded).
    /// Kept for backward compatibility with tblastn/tblastx/psiblast.
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

/// Enumerate all words of length `word_size` that score >= threshold against the query word.
fn enumerate_neighbors(
    query_word: &[u8],
    word_size: usize,
    matrix: &ScoringMatrix,
    threshold: i32,
    callback: &mut impl FnMut(u32),
) {
    // Precompute suffix max scores: max_suffix[i] = max possible score from positions i..word_size
    // This avoids recomputing at every recursive call.
    let mut max_suffix = vec![0i32; word_size + 1];
    for i in (0..word_size).rev() {
        let q = query_word[i];
        let best = (0u8..28).map(|r| matrix.score(q, r)).max().unwrap_or(0);
        max_suffix[i] = max_suffix[i + 1] + best;
    }

    // Precompute score rows: for each query position, the score for each possible residue
    let mut score_rows = vec![[0i32; 28]; word_size];
    for (i, row) in score_rows.iter_mut().enumerate() {
        let q = query_word[i];
        for r in 0u8..28 {
            row[r as usize] = matrix.score(q, r);
        }
    }

    enumerate_rec(&score_rows, &max_suffix, word_size, threshold, 0, 0, 0, callback);
}

#[allow(clippy::too_many_arguments)]
fn enumerate_rec(
    score_rows: &[[i32; 28]],
    max_suffix: &[i32],
    word_size: usize,
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

    // Prune: even with best possible remaining scores, can't reach threshold
    if current_score + max_suffix[pos] < threshold {
        return;
    }

    let scores = &score_rows[pos];
    for r in 0u8..28 {
        let s = scores[r as usize];
        let new_score = current_score + s;
        // Per-residue pruning: check if this choice can still reach threshold
        if new_score + max_suffix[pos + 1] < threshold {
            continue;
        }
        let new_code = current_code * 28 + r as u32;
        enumerate_rec(
            score_rows,
            max_suffix,
            word_size,
            threshold,
            pos + 1,
            new_score,
            new_code,
            callback,
        );
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
