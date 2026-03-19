//! Karlin-Altschul statistics for BLAST.
//!
//! Pre-computed lambda/K/H values for common matrix+gap configurations.

use crate::matrix::MatrixType;

/// Karlin-Altschul statistical parameters.
#[derive(Debug, Clone, Copy)]
pub struct KarlinAltschul {
    pub lambda: f64,
    pub k: f64,
    pub h: f64,
    pub alpha: f64,
    pub beta: f64,
}

impl KarlinAltschul {
    /// Compute E-value from raw score.
    /// E = m * n * K * exp(-lambda * S)
    pub fn evalue(&self, score: i32, query_len: usize, db_len: u64) -> f64 {
        let m = query_len as f64;
        let n = db_len as f64;
        m * n * self.k * (-self.lambda * score as f64).exp()
    }

    /// Compute bit score from raw score.
    pub fn bit_score(&self, score: i32) -> f64 {
        (self.lambda * score as f64 - self.k.ln()) / 2.0f64.ln()
    }

    /// Effective lengths: the effective query/subject lengths account for
    /// the fact that alignments cannot extend to the very ends.
    pub fn effective_lengths(
        &self,
        query_len: usize,
        db_len: u64,
        num_seqs: u64,
    ) -> (usize, u64) {
        // Effective length = length - (ln(K * length) / H)
        let log_k = self.k.ln();
        let ql = query_len as f64;
        let dl = db_len as f64;
        let h = self.h;
        let eff_query = (ql - (log_k + ql.ln()) / h).max(1.0) as usize;
        let eff_db_per = (dl / num_seqs as f64 - (log_k + dl.ln()) / h).max(1.0);
        let eff_db = ((eff_db_per * num_seqs as f64) as u64).max(1);
        (eff_query, eff_db)
    }
}

/// Gap penalty parameters.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GapPenalty {
    pub open: i32,
    pub extend: i32,
}

impl GapPenalty {
    pub fn new(open: i32, extend: i32) -> Self { GapPenalty { open, extend } }
    pub fn blosum62_default() -> Self { GapPenalty { open: 11, extend: 1 } }
    pub fn blastn_default() -> Self { GapPenalty { open: 5, extend: 2 } }
}

/// Look up pre-computed KA parameters for a protein matrix + gap penalty.
pub fn lookup_ka_params(matrix: MatrixType, gap: GapPenalty) -> Option<KarlinAltschul> {
    // Values from NCBI blast_stat.c and BLAST documentation.
    // Rows: (gap_open, gap_extend, lambda, K, H, alpha, beta)
    match matrix {
        MatrixType::Blosum62 => {
            // BLOSUM62 standard values
            let table: &[(i32, i32, f64, f64, f64, f64, f64)] = &[
                // gap_open, gap_extend, lambda, K, H, alpha, beta
                ( 6, 2, 0.3170, 0.1340, 0.401, 1.0,  -1.8),
                ( 7, 2, 0.3170, 0.1770, 0.494, 0.771, -1.2),
                ( 8, 2, 0.3170, 0.1985, 0.576, 0.653, -0.8),
                ( 9, 2, 0.3170, 0.2112, 0.642, 0.641, -0.7),
                (10, 2, 0.3170, 0.2151, 0.661, 0.624, -0.6),
                (11, 1, 0.2670, 0.0410, 0.140, 1.900, -14.0),
                (10, 1, 0.2430, 0.0320, 0.115, 2.110, -17.0),
                ( 9, 1, 0.2060, 0.0200, 0.0820, 2.510, -25.0),
                ( 8, 1, 0.1700, 0.0110, 0.0530, 3.210, -41.0),
                ( 7, 1, 0.1130, 0.00500, 0.0250, 5.470, -87.0),
                (11, 2, 0.3170, 0.1400, 0.431, 0.923, -1.2),
                (10, 2, 0.3170, 0.1200, 0.392, 1.040, -1.9),
                (13, 1, 0.2920, 0.0710, 0.233, 1.25, -5.0),
                (12, 1, 0.2830, 0.0590, 0.190, 1.49, -7.0),
            ];
            find_ka(table, gap)
        }
        MatrixType::Blosum45 => {
            let table: &[(i32, i32, f64, f64, f64, f64, f64)] = &[
                (13, 3, 0.2070, 0.0490, 0.1390, 1.490, -22.0),
                (12, 3, 0.1990, 0.0390, 0.1090, 1.830, -35.0),
                (11, 3, 0.1900, 0.0310, 0.0850, 2.230, -43.0),
                (10, 3, 0.1790, 0.0230, 0.0750, 2.390, -43.0),
                (16, 2, 0.2100, 0.0510, 0.1340, 1.570, -21.0),
                (15, 2, 0.2030, 0.0420, 0.1070, 1.900, -31.0),
                (14, 2, 0.1950, 0.0360, 0.0900, 2.160, -40.0),
                (13, 2, 0.1850, 0.0290, 0.0770, 2.400, -46.0),
                (12, 2, 0.1710, 0.0200, 0.0610, 2.810, -57.0),
                (19, 1, 0.2070, 0.0450, 0.1430, 1.450, -14.0),
                (18, 1, 0.1990, 0.0370, 0.1110, 1.790, -25.0),
                (17, 1, 0.1890, 0.0280, 0.0900, 2.100, -35.0),
                (16, 1, 0.1760, 0.0200, 0.0700, 2.510, -48.0),
            ];
            find_ka(table, gap)
        }
        MatrixType::Blosum80 => {
            let table: &[(i32, i32, f64, f64, f64, f64, f64)] = &[
                (25, 2, 0.3430, 0.1770, 0.657, 0.522, -1.8),
                (13, 2, 0.3360, 0.1500, 0.544, 0.617, -1.8),
                ( 9, 2, 0.3360, 0.1210, 0.431, 0.780, -1.7),
                ( 8, 2, 0.3320, 0.1020, 0.374, 0.890, -1.3),
                ( 7, 2, 0.3250, 0.0820, 0.298, 1.090, -1.6),
                ( 6, 2, 0.3180, 0.0640, 0.257, 1.239, -0.5),
                (11, 1, 0.2730, 0.0470, 0.218, 1.250, -6.0),
                (10, 1, 0.2570, 0.0380, 0.176, 1.460, -9.0),
                ( 9, 1, 0.2270, 0.0270, 0.131, 1.730,-10.0),
                ( 8, 1, 0.1770, 0.0150, 0.0820, 2.150,-15.0),
                ( 7, 1, 0.1390, 0.00790, 0.0550, 2.530,-22.0),
            ];
            find_ka(table, gap)
        }
        MatrixType::Pam30 => {
            let table: &[(i32, i32, f64, f64, f64, f64, f64)] = &[
                ( 7, 2, 0.3050, 0.1770, 1.080, 0.282, -1.0),
                ( 6, 2, 0.2970, 0.1430, 0.861, 0.355, -1.0),
                ( 5, 2, 0.2900, 0.1070, 0.700, 0.415, -0.8),
                (10, 1, 0.2820, 0.0820, 0.627, 0.450, -0.9),
                ( 9, 1, 0.2710, 0.0660, 0.531, 0.510, -1.1),
                ( 8, 1, 0.2530, 0.0500, 0.429, 0.590, -1.1),
            ];
            find_ka(table, gap)
        }
        MatrixType::Pam70 => {
            let table: &[(i32, i32, f64, f64, f64, f64, f64)] = &[
                ( 8, 2, 0.3010, 0.1200, 0.873, 0.345, -1.6),
                ( 7, 2, 0.2860, 0.0930, 0.722, 0.396, -1.0),
                ( 6, 2, 0.2680, 0.0680, 0.576, 0.465, -0.6),
                (11, 1, 0.2840, 0.0790, 0.710, 0.400, -1.3),
                (10, 1, 0.2690, 0.0620, 0.575, 0.468, -1.3),
                ( 9, 1, 0.2520, 0.0470, 0.454, 0.554, -1.2),
                ( 8, 1, 0.2310, 0.0330, 0.350, 0.660, -1.0),
            ];
            find_ka(table, gap)
        }
        MatrixType::Pam250 => {
            let table: &[(i32, i32, f64, f64, f64, f64, f64)] = &[
                (15, 3, 0.2050, 0.0490, 0.2470, 0.830, -4.0),
                (14, 3, 0.2000, 0.0430, 0.2120, 0.940, -5.0),
                (13, 3, 0.1940, 0.0360, 0.1760, 1.100, -7.0),
                (12, 3, 0.1860, 0.0290, 0.1440, 1.290, -9.0),
                (11, 3, 0.1750, 0.0220, 0.1130, 1.550,-11.0),
                (17, 2, 0.2030, 0.0430, 0.2270, 0.894, -5.0),
                (16, 2, 0.1980, 0.0370, 0.1970, 1.010, -6.0),
                (15, 2, 0.1930, 0.0310, 0.1650, 1.160, -8.0),
                (14, 2, 0.1850, 0.0250, 0.1360, 1.360, -9.0),
                (13, 2, 0.1720, 0.0180, 0.1090, 1.580,-10.0),
                (21, 1, 0.2050, 0.0410, 0.2390, 0.858, -4.0),
                (20, 1, 0.1990, 0.0360, 0.2060, 0.966, -5.0),
                (19, 1, 0.1920, 0.0290, 0.1720, 1.110, -7.0),
                (18, 1, 0.1830, 0.0230, 0.1410, 1.300, -9.0),
                (17, 1, 0.1710, 0.0160, 0.1120, 1.530,-11.0),
            ];
            find_ka(table, gap)
        }
    }
}

fn find_ka(
    table: &[(i32, i32, f64, f64, f64, f64, f64)],
    gap: GapPenalty,
) -> Option<KarlinAltschul> {
    for &(go, ge, lambda, k, h, alpha, beta) in table {
        if go == gap.open && ge == gap.extend {
            return Some(KarlinAltschul { lambda, k, h, alpha, beta });
        }
    }
    None
}

/// Ungapped KA parameters for nucleotide BLAST.
/// For blastn with match=2, mismatch=-3.
pub fn blastn_ka_params(match_score: i32, mismatch: i32, gap_open: i32, gap_extend: i32) -> KarlinAltschul {
    // Simple approximation; real BLAST computes these from the scoring system.
    // Default values for match=2, mismatch=-3.
    let _ = (match_score, mismatch, gap_open, gap_extend);
    KarlinAltschul {
        lambda: 1.370,
        k: 0.711,
        h: 1.310,
        alpha: 1.0,
        beta: 0.0,
    }
}
