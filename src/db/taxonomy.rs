//! NCBI taxonomy database reader (taxdb.bti/btd files).
//!
//! Provides taxid → organism name resolution for tabular output columns
//! (ssciname, scomname, sblastname, sskingdom) and the organism report (format 18).

use std::path::{Path, PathBuf};
use std::fs;

/// Taxonomy information for one taxid.
#[derive(Debug, Clone, Default)]
pub struct TaxInfo {
    pub scientific_name: String,
    pub common_name: String,
    pub blast_name: String,
    pub kingdom: String,
}

/// Reader for NCBI taxdb.bti (index) + taxdb.btd (data) files.
/// Uses binary search on the index to resolve taxid → names.
pub struct TaxDb {
    /// Sorted array of (taxid, offset) pairs from taxdb.bti
    index: Vec<(u32, u32)>,
    /// Raw data from taxdb.btd (tab-separated name strings)
    data: Vec<u8>,
}

impl TaxDb {
    /// Try to open taxdb from standard search paths.
    /// Returns None if taxdb files are not found.
    pub fn open() -> Option<Self> {
        let paths = taxdb_search_paths();
        for dir in &paths {
            let bti = dir.join("taxdb.bti");
            let btd = dir.join("taxdb.btd");
            if bti.exists() && btd.exists() {
                if let Ok(db) = Self::open_from(&bti, &btd) {
                    return Some(db);
                }
            }
        }
        None
    }

    /// Open taxdb from specific file paths.
    pub fn open_from(bti_path: &Path, btd_path: &Path) -> Result<Self, String> {
        let bti_data = fs::read(bti_path)
            .map_err(|e| format!("Failed to read {}: {}", bti_path.display(), e))?;
        let data = fs::read(btd_path)
            .map_err(|e| format!("Failed to read {}: {}", btd_path.display(), e))?;

        // Parse bti header: magic (4 bytes LE) + num_entries (4 bytes LE) + 16 bytes reserved
        if bti_data.len() < 24 {
            return Err("taxdb.bti too small".into());
        }

        let magic = u32::from_le_bytes([bti_data[0], bti_data[1], bti_data[2], bti_data[3]]);
        if magic != 0x8739 {
            return Err(format!("taxdb.bti bad magic: {:#x}", magic));
        }

        let num_entries = u32::from_le_bytes([bti_data[4], bti_data[5], bti_data[6], bti_data[7]]) as usize;

        // Parse index entries: (taxid: u32, offset: u32) pairs starting at byte 24
        let entry_size = 8; // 4 bytes taxid + 4 bytes offset
        let expected_size = 24 + num_entries * entry_size;
        if bti_data.len() < expected_size {
            return Err(format!("taxdb.bti truncated: {} < {}", bti_data.len(), expected_size));
        }

        let mut index = Vec::with_capacity(num_entries);
        for i in 0..num_entries {
            let base = 24 + i * entry_size;
            let taxid = u32::from_le_bytes([
                bti_data[base], bti_data[base + 1], bti_data[base + 2], bti_data[base + 3],
            ]);
            let offset = u32::from_le_bytes([
                bti_data[base + 4], bti_data[base + 5], bti_data[base + 6], bti_data[base + 7],
            ]);
            index.push((taxid, offset));
        }

        Ok(TaxDb { index, data })
    }

    /// Look up taxonomy info for a taxid using binary search.
    pub fn lookup(&self, taxid: i32) -> Option<TaxInfo> {
        let taxid = taxid as u32;
        let pos = self.index.binary_search_by_key(&taxid, |&(t, _)| t).ok()?;
        let offset = self.index[pos].1 as usize;

        // Find end of this entry's data
        let end = if pos + 1 < self.index.len() {
            self.index[pos + 1].1 as usize
        } else {
            self.data.len()
        };

        if offset >= self.data.len() || end > self.data.len() {
            return None;
        }

        // Parse tab-separated: scientific_name\tcommon_name\tblast_name\tkingdom
        let record = std::str::from_utf8(&self.data[offset..end]).ok()?;
        let record = record.trim_end_matches('\0').trim();
        let fields: Vec<&str> = record.split('\t').collect();

        Some(TaxInfo {
            scientific_name: fields.first().unwrap_or(&"").to_string(),
            common_name: fields.get(1).unwrap_or(&"").to_string(),
            blast_name: fields.get(2).unwrap_or(&"").to_string(),
            kingdom: fields.get(3).unwrap_or(&"").to_string(),
        })
    }
}

/// Search paths for taxdb files, matching NCBI's lookup order.
fn taxdb_search_paths() -> Vec<PathBuf> {
    let mut paths = Vec::new();

    // $BLASTDB environment variable
    if let Ok(blastdb) = std::env::var("BLASTDB") {
        for dir in blastdb.split(':') {
            paths.push(PathBuf::from(dir));
        }
    }

    // ~/.ncbi/
    if let Some(home) = std::env::var_os("HOME") {
        paths.push(PathBuf::from(home).join(".ncbi"));
    }

    // System paths
    paths.push(PathBuf::from("/usr/local/ncbi/data"));
    paths.push(PathBuf::from("/usr/share/ncbi"));

    // Current directory
    paths.push(PathBuf::from("."));

    paths
}
