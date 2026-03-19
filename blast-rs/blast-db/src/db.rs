use std::path::Path;
use std::fs;
use memmap2::Mmap;
use crate::error::{DbError, Result};
use crate::index::{IndexFile, SeqType};
use crate::header::{BlastDefLine, parse_def_line_set};
use crate::sequence::{get_protein_raw, get_nucleotide, decode_protein};
use crate::lmdb_v5::LmdbV5;
use crate::oid_seqids::OidSeqIds;
use crate::oid_taxids::OidTaxIds;

pub struct BlastDb {
    pub(crate) index: IndexFile,
    /// Memory-mapped sequence file (.psq or .nsq)
    seq_mmap: Mmap,
    /// Memory-mapped header file (.phr or .nhr)
    hdr_mmap: Mmap,
    /// V5 LMDB accession index (.pdb or .ndb) — None for v4 databases.
    lmdb: Option<LmdbV5>,
    /// V5 OID→SeqIds file (.pos or .nos) — None for v4 databases.
    oid_seqids: Option<OidSeqIds>,
    /// V5 OID→TaxIds file (.pot or .not) — None for v4 databases.
    oid_taxids: Option<OidTaxIds>,
}

impl BlastDb {
    /// Open a BLAST database by base path (without extension).
    /// Auto-detects protein vs nucleotide and v4 vs v5.
    pub fn open(path: &Path) -> Result<Self> {
        // Detect type by which index file exists.
        let (is_protein, index_ext, seq_ext, hdr_ext) = if path.with_extension("pin").exists() {
            (true, "pin", "psq", "phr")
        } else if path.with_extension("nin").exists() {
            (false, "nin", "nsq", "nhr")
        } else {
            return Err(DbError::InvalidFormat(format!(
                "No .pin or .nin index file found at {}",
                path.display()
            )));
        };

        let index_data = fs::read(path.with_extension(index_ext))?;
        let index = IndexFile::parse(&index_data)?;

        let seq_file = fs::File::open(path.with_extension(seq_ext))?;
        let seq_mmap = unsafe { Mmap::map(&seq_file)? };

        let hdr_file = fs::File::open(path.with_extension(hdr_ext))?;
        let hdr_mmap = unsafe { Mmap::map(&hdr_file)? };

        // V5: try to open LMDB and auxiliary files (gracefully absent = v4).
        let (lmdb, oid_seqids, oid_taxids) = if index.format_version == 5 {
            let (lmdb_ext, seqids_ext, taxids_ext) =
                if is_protein { ("pdb", "pos", "pot") }
                else          { ("ndb", "nos", "not") };

            let lmdb_path = path.with_extension(lmdb_ext);
            let lmdb = if lmdb_path.exists() {
                Some(LmdbV5::open(&lmdb_path)?)
            } else {
                None
            };

            let oid_seqids_path = path.with_extension(seqids_ext);
            let oid_seqids = if oid_seqids_path.exists() {
                Some(OidSeqIds::open(&oid_seqids_path)?)
            } else {
                None
            };

            let oid_taxids_path = path.with_extension(taxids_ext);
            let oid_taxids = if oid_taxids_path.exists() {
                Some(OidTaxIds::open(&oid_taxids_path)?)
            } else {
                None
            };

            (lmdb, oid_seqids, oid_taxids)
        } else {
            (None, None, None)
        };

        Ok(BlastDb { index, seq_mmap, hdr_mmap, lmdb, oid_seqids, oid_taxids })
    }

    // ---- Metadata ----

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

    /// Returns 4 for classic BlastDB, 5 for LMDB-indexed BlastDB.
    pub fn format_version(&self) -> i32 {
        self.index.format_version
    }

    /// Returns true if this is a v5 database with an LMDB accession index.
    pub fn is_v5(&self) -> bool {
        self.index.format_version == 5
    }

    // ---- Sequence access ----

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

    // ---- Header access ----

    /// Returns the primary defline for an OID (from the BER-encoded .phr/.nhr file).
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

    // ---- V5 accession index ----

    /// Find OID(s) for an accession string using the LMDB index.
    /// Returns `None` if this is a v4 database (no LMDB index).
    /// Returns an empty `Vec` if the accession is not found.
    pub fn lookup_accession(&self, accession: &str) -> Option<Result<Vec<u32>>> {
        let lmdb = self.lmdb.as_ref()?;
        Some(lmdb.get_oids_for_accession(accession))
    }

    /// Iterate over all (accession, oid) pairs in the LMDB index.
    /// Returns `None` if this is a v4 database.
    pub fn iter_accessions<F>(&self, f: F) -> Option<Result<()>>
    where
        F: FnMut(&str, u32),
    {
        let lmdb = self.lmdb.as_ref()?;
        Some(lmdb.iter_accessions(f))
    }

    /// Get volume name → num_oids information from the LMDB.
    /// Returns `None` for v4 databases.
    pub fn get_volumes_info(&self) -> Option<Result<Vec<(String, u32)>>> {
        let lmdb = self.lmdb.as_ref()?;
        Some(lmdb.get_volumes_info())
    }

    // ---- V5 OID→SeqIDs ----

    /// Return all seq-id strings for an OID from the `.pos`/`.nos` file.
    /// Returns `None` if the file is not available (v4 or absent).
    pub fn get_seqids(&self, oid: u32) -> Option<Result<Vec<String>>> {
        let r = self.oid_seqids.as_ref()?;
        Some(r.get_seqids(oid))
    }

    // ---- V5 OID→TaxIDs ----

    /// Return all tax IDs for an OID from the `.pot`/`.not` file.
    /// Returns `None` if the file is not available (v4 or absent).
    pub fn get_taxids(&self, oid: u32) -> Option<Result<Vec<i32>>> {
        let r = self.oid_taxids.as_ref()?;
        Some(r.get_taxids(oid))
    }

    // ---- Internal ----

    fn check_oid(&self, oid: u32) -> Result<()> {
        if oid >= self.index.num_oids {
            Err(DbError::OidOutOfRange(oid))
        } else {
            Ok(())
        }
    }
}
