use std::path::{Path, PathBuf};
use std::fs;
use memmap2::Mmap;
use crate::error::{DbError, Result};
use crate::index::{IndexFile, SeqType};
use crate::header::{BlastDefLine, parse_def_line_set};
use crate::sequence::{get_protein_raw, get_nucleotide, decode_protein};

pub struct BlastDb {
    pub(crate) index: IndexFile,
    /// Memory-mapped sequence file (.psq or .nsq)
    seq_mmap: Mmap,
    /// Memory-mapped header file (.phr or .nhr)
    hdr_mmap: Mmap,
}

impl BlastDb {
    /// Open a BLAST database by base path (without extension).
    /// Detects whether it's protein or nucleotide from available files.
    pub fn open(path: &Path) -> Result<Self> {
        // Try protein first, then nucleotide
        let (index_ext, seq_ext, hdr_ext) = if path.with_extension("pin").exists() {
            ("pin", "psq", "phr")
        } else if path.with_extension("nin").exists() {
            ("nin", "nsq", "nhr")
        } else {
            return Err(DbError::InvalidFormat(format!(
                "No .pin or .nin file found at {}",
                path.display()
            )));
        };

        let index_path = path.with_extension(index_ext);
        let seq_path = path.with_extension(seq_ext);
        let hdr_path = path.with_extension(hdr_ext);

        let index_data = fs::read(&index_path)?;
        let index = IndexFile::parse(&index_data)?;

        let seq_file = fs::File::open(&seq_path)?;
        let seq_mmap = unsafe { Mmap::map(&seq_file)? };

        let hdr_file = fs::File::open(&hdr_path)?;
        let hdr_mmap = unsafe { Mmap::map(&hdr_file)? };

        Ok(BlastDb { index, seq_mmap, hdr_mmap })
    }

    pub fn num_sequences(&self) -> u32 {
        self.index.num_oids
    }

    pub fn seq_type(&self) -> SeqType {
        self.index.seq_type
    }

    pub fn volume_length(&self) -> u64 {
        self.index.volume_length
    }

    pub fn title(&self) -> &str {
        &self.index.title
    }

    /// Returns raw Ncbistdaa bytes for a protein sequence.
    pub fn get_sequence_protein_raw(&self, oid: u32) -> Result<&[u8]> {
        self.check_oid(oid)?;
        let start = self.index.sequence_array[oid as usize] as usize;
        let end = self.index.sequence_array[oid as usize + 1] as usize;
        Ok(get_protein_raw(&self.seq_mmap, start, end))
    }

    /// Returns ASCII amino acid sequence for a protein OID.
    pub fn get_sequence_protein(&self, oid: u32) -> Result<Vec<u8>> {
        Ok(decode_protein(self.get_sequence_protein_raw(oid)?))
    }

    /// Returns decoded ASCII nucleotide sequence for a nucleotide OID.
    pub fn get_sequence_nucleotide(&self, oid: u32) -> Result<Vec<u8>> {
        self.check_oid(oid)?;
        let ambig = self.index.ambig_array.as_ref()
            .ok_or_else(|| DbError::InvalidFormat("No ambig_array for nucleotide db".into()))?;
        let seq_start = self.index.sequence_array[oid as usize] as usize;
        let seq_end = ambig[oid as usize] as usize;
        let ambig_start = seq_end;
        let ambig_end = self.index.sequence_array[oid as usize + 1] as usize;
        Ok(get_nucleotide(&self.seq_mmap, seq_start, seq_end, ambig_start, ambig_end))
    }

    /// Returns the primary defline for an OID.
    pub fn get_header(&self, oid: u32) -> Result<BlastDefLine> {
        self.check_oid(oid)?;
        let start = self.index.header_array[oid as usize] as usize;
        let end = self.index.header_array[oid as usize + 1] as usize;
        let data = &self.hdr_mmap[start..end];
        let deflines = parse_def_line_set(data)?;
        Ok(deflines.into_iter().next().unwrap_or_default())
    }

    /// Returns all deflines for an OID (for non-redundant databases).
    pub fn get_headers(&self, oid: u32) -> Result<Vec<BlastDefLine>> {
        self.check_oid(oid)?;
        let start = self.index.header_array[oid as usize] as usize;
        let end = self.index.header_array[oid as usize + 1] as usize;
        let data = &self.hdr_mmap[start..end];
        parse_def_line_set(data)
    }

    fn check_oid(&self, oid: u32) -> Result<()> {
        if oid >= self.index.num_oids {
            Err(DbError::OidOutOfRange(oid))
        } else {
            Ok(())
        }
    }
}
