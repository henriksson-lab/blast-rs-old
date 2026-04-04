//! Minimal BER decoder for Blast-def-line-set ASN.1 objects.
//!
//! We only decode the fields we need: title, seqid (accession), taxid.
//!
//! ASN.1 schema:
//!   Blast-def-line-set ::= SEQUENCE OF Blast-def-line
//!   Blast-def-line ::= SEQUENCE {
//!       title VisibleString OPTIONAL,
//!       seqid SEQUENCE OF Seq-id,
//!       taxid INTEGER OPTIONAL,
//!       memberships SEQUENCE OF INTEGER OPTIONAL,
//!       links SEQUENCE OF INTEGER OPTIONAL,
//!       other-info SEQUENCE OF INTEGER OPTIONAL
//!   }

use crate::db::error::{DbError, Result};

/// A parsed defline from the BLAST database header.
#[derive(Debug, Clone, Default)]
pub struct BlastDefLine {
    pub title: String,
    /// Primary accession string (first seqid that looks useful)
    pub accession: String,
    pub taxid: u32,
}

// ---- BER primitive helpers ----

/// Returns (tag_byte, length, data_slice, rest)
fn read_tlv(data: &[u8]) -> Result<(u8, &[u8], &[u8])> {
    if data.is_empty() {
        return Err(DbError::AsnParse("unexpected end of data".into()));
    }
    let tag = data[0];
    let mut pos = 1usize;

    let len = if data[pos] & 0x80 == 0 {
        let l = data[pos] as usize;
        pos += 1;
        l
    } else {
        let num_bytes = (data[pos] & 0x7f) as usize;
        pos += 1;
        if num_bytes == 0 || num_bytes > 4 {
            return Err(DbError::AsnParse(format!("unsupported length encoding: num_bytes={}", num_bytes)));
        }
        let mut l = 0usize;
        for _ in 0..num_bytes {
            if pos >= data.len() {
                return Err(DbError::AsnParse("truncated length".into()));
            }
            l = (l << 8) | data[pos] as usize;
            pos += 1;
        }
        l
    };

    if pos + len > data.len() {
        return Err(DbError::AsnParse(format!(
            "data too short: need {} bytes, have {}",
            len,
            data.len() - pos
        )));
    }

    Ok((tag, &data[pos..pos + len], &data[pos + len..]))
}

fn decode_integer(data: &[u8]) -> i64 {
    if data.is_empty() {
        return 0;
    }
    let sign_extend = data[0] as i8 as i64;
    let mut val = sign_extend;
    for &b in &data[1..] {
        val = (val << 8) | b as i64;
    }
    val
}

/// Parse a Seq-id CHOICE and return a string representation.
/// Seq-id is a CHOICE with many alternatives; we handle the most common ones.
fn parse_seq_id(data: &[u8]) -> String {
    // Seq-id is a CHOICE — the tag tells us which alternative
    if data.is_empty() {
        return String::new();
    }
    let tag = data[0];
    let rest = match read_tlv(data) {
        Ok((_, content, _)) => content,
        Err(_) => return String::new(),
    };

    match tag {
        // [0] local — Object-id (CHOICE string/id)
        0xa0 => parse_object_id(rest),
        // [1] gibbsq — INTEGER
        0xa1 => format!("gibbsq:{}", decode_integer(rest)),
        // [2] gibbmt — INTEGER
        0xa2 => format!("gibbmt:{}", decode_integer(rest)),
        // [3] giim (not common)
        // [4] embl
        0xa4 => parse_textseq_id(rest, "embl"),
        // [5] pir
        0xa5 => parse_textseq_id(rest, "pir"),
        // [6] swissprot
        0xa6 => parse_textseq_id(rest, "sp"),
        // [7] patent
        // [8] other = refseq
        0xa8 => parse_textseq_id(rest, "ref"),
        // [9] general — Dbtag
        0xa9 => parse_dbtag(rest),
        // [10] gi — INTEGER
        0xaa => {
            // gi: context [10] IMPLICIT INTEGER
            if let Ok((_, content, _)) = read_tlv(data) {
                format!("gi|{}", decode_integer(content))
            } else {
                String::new()
            }
        },
        // [11] ddbj
        0xab => parse_textseq_id(rest, "dbj"),
        // [12] prf
        0xac => parse_textseq_id(rest, "prf"),
        // [13] pdb
        // [14] tpg — genbank third party
        0xae => parse_textseq_id(rest, "tpg"),
        // [15] tpe
        0xaf => parse_textseq_id(rest, "tpe"),
        // [16] tpd
        0xb0 => parse_textseq_id(rest, "tpd"),
        // [17] gpipe
        0xb1 => parse_textseq_id(rest, "gpipe"),
        // [18] named-annot-track
        _ => {
            // Try treating it as a textseq-id if it looks like SEQUENCE
            if tag & 0x20 != 0 {
                parse_textseq_id(rest, "?")
            } else {
                format!("unknown_tag:{:#x}", tag)
            }
        }
    }
}

/// Parse Object-id ::= CHOICE { id INTEGER, str VisibleString }
fn parse_object_id(data: &[u8]) -> String {
    if data.is_empty() {
        return String::new();
    }
    match read_tlv(data) {
        Ok((0x02, content, _)) => format!("local:{}", decode_integer(content)),
        Ok((0x1a, content, _)) | Ok((0x13, content, _)) => {
            String::from_utf8_lossy(content).into_owned()
        }
        _ => String::new(),
    }
}

/// Parse Dbtag ::= SEQUENCE { db VisibleString, tag Object-id }
fn parse_dbtag(data: &[u8]) -> String {
    let mut cur = data;
    let mut db = String::new();
    let mut tag = String::new();
    while !cur.is_empty() {
        match read_tlv(cur) {
            Ok((0x1a, content, rest)) | Ok((0x13, content, rest)) => {
                if db.is_empty() {
                    db = String::from_utf8_lossy(content).into_owned();
                }
                cur = rest;
            }
            Ok((0xa0, content, rest)) => {
                tag = parse_object_id(content);
                cur = rest;
            }
            Ok((_, _, rest)) => { cur = rest; }
            Err(_) => break,
        }
    }
    if db.is_empty() {
        tag
    } else {
        format!("{}:{}", db, tag)
    }
}

/// Parse Text-seq-id ::= SEQUENCE {
///   name VisibleString OPTIONAL,
///   accession VisibleString OPTIONAL,
///   release VisibleString OPTIONAL,
///   version INTEGER OPTIONAL
/// }
fn parse_textseq_id(data: &[u8], prefix: &str) -> String {
    let mut cur = data;
    let mut name = String::new();
    let mut accession = String::new();
    let mut version: Option<i64> = None;
    while !cur.is_empty() {
        match read_tlv(cur) {
            Ok((tag, content, rest)) => {
                match tag {
                    0x1a | 0x13 => {  // VisibleString / PrintableString
                        let s = String::from_utf8_lossy(content).into_owned();
                        if name.is_empty() {
                            name = s;
                        } else if accession.is_empty() {
                            accession = s;
                        }
                    }
                    0x02 => {  // INTEGER = version
                        version = Some(decode_integer(content));
                    }
                    _ => {}
                }
                cur = rest;
            }
            Err(_) => break,
        }
    }
    // prefer accession; fall back to name
    let id = if !accession.is_empty() { accession } else { name };
    match version {
        Some(v) => format!("{}|{}.{}", prefix, id, v),
        None => format!("{}|{}", prefix, id),
    }
}

/// Parse one Blast-def-line SEQUENCE.
fn parse_def_line(data: &[u8]) -> Result<BlastDefLine> {
    let mut cur = data;
    let mut dl = BlastDefLine::default();

    // Fields are CONTEXT-tagged (implicit) in order:
    // [0] title VisibleString OPTIONAL
    // [1] seqid SEQUENCE OF Seq-id
    // [2] taxid INTEGER OPTIONAL
    // (higher tags: memberships, links, other-info — we skip)

    while !cur.is_empty() {
        let (tag, content, rest) = match read_tlv(cur) {
            Ok(x) => x,
            Err(_) => break,
        };
        cur = rest;

        match tag {
            // [0] IMPLICIT VisibleString — title
            // Context [0] constructed = 0xa0, primitive = 0x80
            0xa0 => {
                // constructed context [0] — title as VisibleString inside
                if let Ok((_, inner, _)) = read_tlv(content) {
                    dl.title = String::from_utf8_lossy(inner).into_owned();
                } else {
                    dl.title = String::from_utf8_lossy(content).into_owned();
                }
            }
            0x80 => {
                // primitive context [0] — direct string
                dl.title = String::from_utf8_lossy(content).into_owned();
            }
            // [1] IMPLICIT SEQUENCE OF Seq-id — seqid
            0xa1 => {
                let mut seqid_cur = content;
                let mut first = true;
                while !seqid_cur.is_empty() {
                    let (_, _, seqid_rest) = match read_tlv(seqid_cur) {
                        Ok(x) => x,
                        Err(_) => break,
                    };
                    if first {
                        dl.accession = parse_seq_id(seqid_cur);
                        first = false;
                    }
                    seqid_cur = seqid_rest;
                }
            }
            // [2] IMPLICIT INTEGER — taxid
            0xa2 => {
                // constructed [2] containing an INTEGER
                if let Ok((0x02, int_content, _)) = read_tlv(content) {
                    dl.taxid = decode_integer(int_content) as u32;
                } else {
                    dl.taxid = decode_integer(content) as u32;
                }
            }
            0x82 => {
                // primitive [2] — taxid
                dl.taxid = decode_integer(content) as u32;
            }
            // skip everything else
            _ => {}
        }
    }

    Ok(dl)
}

/// Parse a Blast-def-line-set from raw BER bytes.
/// Returns all deflines found; the first is the primary one.
pub fn parse_def_line_set(data: &[u8]) -> Result<Vec<BlastDefLine>> {
    // The outer wrapper is: SEQUENCE { ... } or SET { ... }
    // Actually Blast-def-line-set ::= SEQUENCE OF Blast-def-line
    // The encoding is: 0x30 (SEQUENCE) len [Blast-def-line ...]+
    if data.is_empty() {
        return Err(DbError::AsnParse("empty header data".into()));
    }

    let (tag, content, _) = read_tlv(data)?;
    if tag != 0x30 {
        return Err(DbError::AsnParse(format!("expected SEQUENCE (0x30), got {:#x}", tag)));
    }

    let mut cur = content;
    let mut result = Vec::new();

    while !cur.is_empty() {
        let (tag2, inner, rest) = match read_tlv(cur) {
            Ok(x) => x,
            Err(_) => break,
        };
        cur = rest;
        if tag2 == 0x30 {
            if let Ok(dl) = parse_def_line(inner) {
                result.push(dl);
            }
        }
    }

    if result.is_empty() {
        // Try parsing as a single defline without the outer set wrapper
        if let Ok(dl) = parse_def_line(content) {
            result.push(dl);
        }
    }

    Ok(result)
}
