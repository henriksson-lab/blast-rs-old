#[cfg(test)]
mod tests {
    use crate::matrix::{ScoringMatrix, MatrixType};
    use crate::stats::{KarlinAltschul, GapPenalty, lookup_ka_params};

    #[test]
    fn test_blosum62_diagonal() {
        let m = ScoringMatrix::blosum62();
        // A=1 in Ncbistdaa, BLOSUM62 A-A = 4
        assert_eq!(m.score(1, 1), 4, "A-A score");
        // W=20 in Ncbistdaa, BLOSUM62 W-W = 11
        assert_eq!(m.score(20, 20), 11, "W-W score");
        // * = 25, BLOSUM62 *-* = 1
        assert_eq!(m.score(25, 25), 1, "*-* score");
    }

    #[test]
    fn test_blosum62_off_diagonal() {
        let m = ScoringMatrix::blosum62();
        // A(1) vs R(16) in BLOSUM62 = -1
        assert_eq!(m.score(1, 16), -1, "A-R score");
        // W(20) vs Y(22) = 2
        assert_eq!(m.score(20, 22), 2, "W-Y score");
    }

    #[test]
    fn test_ka_params_lookup() {
        let ka = lookup_ka_params(MatrixType::Blosum62, GapPenalty::new(11, 1));
        assert!(ka.is_some(), "Should find KA params for BLOSUM62 gap=11,1");
        let ka = ka.unwrap();
        assert!((ka.lambda - 0.267).abs() < 0.001, "lambda ~0.267, got {}", ka.lambda);
        assert!((ka.k - 0.041).abs() < 0.001, "K ~0.041, got {}", ka.k);
    }

    #[test]
    fn test_evalue_calculation() {
        let ka = KarlinAltschul { lambda: 0.267, k: 0.041, h: 0.14, alpha: 1.9, beta: -14.0 };
        // Short alignment: query=100, db=1_000_000, score=100
        let ev = ka.evalue(100, 100, 1_000_000);
        // Should be a very small number: 100 * 1e6 * 0.041 * exp(-0.267*100) ~ 1.7e-7
        assert!(ev < 0.001, "E-value for score 100 should be < 0.001, got {}", ev);
        // score=50 with these params gives ~6.5 (above threshold) - correct behavior
        let ev50 = ka.evalue(50, 100, 1_000_000);
        assert!(ev50 > 1.0, "E-value for score 50 expected > 1, got {}", ev50);
        // Bit score
        let bs = ka.bit_score(100);
        assert!(bs > 20.0, "bit score should be > 20, got {}", bs);
    }
}
