use byteorder::{BigEndian, LittleEndian, ReadBytesExt};
use std::io::{Cursor, Read};
use crate::error::{DbError, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqType {
    Nucleotide,
    Protein,
}

#[derive(Debug)]
pub struct IndexFile {
    pub seq_type: SeqType,
    pub title: String,
    pub create_date: String,
    pub num_oids: u32,
    pub volume_length: u64,
    pub max_seq_length: u32,
    /// header-array: (num_oids+1) entries, byte offsets into .nhr/.phr
    pub header_array: Vec<u32>,
    /// sequence-array: (num_oids+1) entries, byte offsets into .nsq/.psq
    pub sequence_array: Vec<u32>,
    /// ambig-array: (num_oids+1) entries, only present for nucleotide
    pub ambig_array: Option<Vec<u32>>,
}

fn read_string(cur: &mut Cursor<&[u8]>) -> Result<String> {
    let len = cur.read_i32::<BigEndian>()? as usize;
    let mut buf = vec![0u8; len];
    cur.read_exact(&mut buf)?;
    // strip trailing NUL bytes
    while buf.last() == Some(&0) {
        buf.pop();
    }
    Ok(String::from_utf8_lossy(&buf).into_owned())
}

impl IndexFile {
    pub fn parse(data: &[u8]) -> Result<Self> {
        let mut cur = Cursor::new(data);

        let version = cur.read_i32::<BigEndian>()?;
        if version != 4 {
            return Err(DbError::UnsupportedVersion(version));
        }

        let seq_type_raw = cur.read_i32::<BigEndian>()?;
        let seq_type = match seq_type_raw {
            0 => SeqType::Nucleotide,
            1 => SeqType::Protein,
            _ => return Err(DbError::InvalidFormat(format!("Unknown seq-type {}", seq_type_raw))),
        };

        let title = read_string(&mut cur)?;
        let create_date = read_string(&mut cur)?;

        let num_oids = cur.read_u32::<BigEndian>()?;
        let volume_length = cur.read_u64::<LittleEndian>()?;
        let max_seq_length = cur.read_u32::<BigEndian>()?;

        let count = (num_oids + 1) as usize;

        let mut header_array = Vec::with_capacity(count);
        for _ in 0..count {
            header_array.push(cur.read_u32::<BigEndian>()?);
        }

        let mut sequence_array = Vec::with_capacity(count);
        for _ in 0..count {
            sequence_array.push(cur.read_u32::<BigEndian>()?);
        }

        let ambig_array = if seq_type == SeqType::Nucleotide {
            let mut arr = Vec::with_capacity(count);
            for _ in 0..count {
                arr.push(cur.read_u32::<BigEndian>()?);
            }
            Some(arr)
        } else {
            None
        };

        Ok(IndexFile {
            seq_type,
            title,
            create_date,
            num_oids,
            volume_length,
            max_seq_length,
            header_array,
            sequence_array,
            ambig_array,
        })
    }
}
