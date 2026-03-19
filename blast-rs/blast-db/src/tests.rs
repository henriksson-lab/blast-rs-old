#[cfg(test)]
mod tests {
    use crate::sequence::{decode_protein, decode_nucleotide, NCBISTDAA_TO_AA};
    use crate::header::parse_def_line_set;

    #[test]
    fn test_protein_decode_roundtrip() {
        // Ncbistdaa code 1 = A, 11 = L, 5 = E
        let raw = &[1u8, 11, 5, 7, 16]; // A L E G R
        let decoded = decode_protein(raw);
        assert_eq!(&decoded, b"ALEGR");
    }

    #[test]
    fn test_nucleotide_decode_tggttacaac() {
        // From spec example: TGGTTACAAC -> (235, 196, 18)
        // T=3, G=2, G=2, T=3, T=3, A=0, C=1, A=0, A=0, C=1
        // Full bytes: 0xEB (11101011), 0xC4 (11000100), 0x12 (00010010)
        // 0x12 = 00 01 00 10 -> base A(00), C(01), then count=10(2)
        let packed = &[0xEB, 0xC4, 0x12];
        let decoded = decode_nucleotide(packed, &[]);
        assert_eq!(&decoded, b"TGGTTACAAC", "Got: {:?}", std::str::from_utf8(&decoded).unwrap_or("?"));
    }

    #[test]
    fn test_nucleotide_decode_tacg() {
        // TACG: T=3, A=0, C=1, G=2
        // One byte: (3<<6)|(0<<4)|(1<<2)|2 = 198 = 0xC6
        // Plus sentinel byte 0 (remainder=0, meaning 0 bases in sentinel)
        let packed = &[0xC6u8, 0x00];
        let decoded = decode_nucleotide(packed, &[]);
        assert_eq!(&decoded, b"TACG", "Got: {:?}", std::str::from_utf8(&decoded).unwrap_or("?"));
    }

    #[test]
    fn test_nucleotide_single_base_a() {
        // A: 0b00_00_00_01 with remainder=1
        // bits 7,6 = A(0,0); bits 1,0 = count(01)=1
        // byte = 0b00_000001 = 0x01
        let packed = &[0x01u8];
        let decoded = decode_nucleotide(packed, &[]);
        assert_eq!(&decoded, b"A", "Got: {:?}", std::str::from_utf8(&decoded).unwrap_or("?"));
    }

    #[test]
    fn test_ncbistdaa_table() {
        // Standard checks
        assert_eq!(NCBISTDAA_TO_AA[1], b'A');
        assert_eq!(NCBISTDAA_TO_AA[11], b'L');
        assert_eq!(NCBISTDAA_TO_AA[25], b'*');
    }
}
