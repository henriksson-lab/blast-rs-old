//! BlastDB version 5 LMDB reader.
//!
//! V5 databases store accession→OID mappings in an LMDB file (`.pdb`/`.ndb`).
//!
//! Named databases inside the LMDB environment:
//!   "acc2oid"  — MDB_DUPSORT|MDB_DUPFIXED; key: accession bytes; value: LE u32 OID
//!   "volname"  — MDB_INTEGERKEY; key: native-endian u32 index; value: volume name bytes
//!   "volinfo"  — MDB_INTEGERKEY; key: native-endian u32 index; value: LE u32 num_oids

use std::path::Path;
use lmdb::{self, Cursor, Database, Environment, EnvironmentFlags, Transaction};
use crate::db::error::{DbError, Result};

pub struct LmdbV5 {
    env: Environment,
    db_acc2oid: Database,
    db_volname: Database,
    db_volinfo: Database,
}

impl LmdbV5 {
    /// Open a BlastDB v5 LMDB file path directly (not a directory; uses NO_SUB_DIR).
    pub fn open(path: &Path) -> Result<Self> {
        let env = Environment::new()
            .set_flags(
                EnvironmentFlags::NO_SUB_DIR
                    | EnvironmentFlags::NO_LOCK
                    | EnvironmentFlags::READ_ONLY,
            )
            .set_max_dbs(4)
            // Large virtual map; uses no physical RAM until pages are accessed.
            .set_map_size(1 << 40)
            .open(path)
            .map_err(|e| DbError::InvalidFormat(
                format!("LMDB open '{}': {}", path.display(), e)
            ))?;

        // Named databases must be opened inside a transaction.
        let txn = env.begin_ro_txn()
            .map_err(|e| DbError::InvalidFormat(format!("LMDB txn: {}", e)))?;

        let db_acc2oid = unsafe {
            txn.open_db(Some("acc2oid"))
                .map_err(|e| DbError::InvalidFormat(format!("LMDB open 'acc2oid': {}", e)))?
        };
        let db_volname = unsafe {
            txn.open_db(Some("volname"))
                .map_err(|e| DbError::InvalidFormat(format!("LMDB open 'volname': {}", e)))?
        };
        let db_volinfo = unsafe {
            txn.open_db(Some("volinfo"))
                .map_err(|e| DbError::InvalidFormat(format!("LMDB open 'volinfo': {}", e)))?
        };
        txn.abort();

        Ok(LmdbV5 { env, db_acc2oid, db_volname, db_volinfo })
    }

    /// Look up all OIDs for an accession string.
    /// Returns an empty Vec if not found.
    pub fn get_oids_for_accession(&self, accession: &str) -> Result<Vec<u32>> {
        let txn = self.env.begin_ro_txn()
            .map_err(|e| DbError::InvalidFormat(format!("LMDB txn: {}", e)))?;
        let mut cursor = txn.open_ro_cursor(self.db_acc2oid)
            .map_err(|e| DbError::InvalidFormat(format!("LMDB cursor: {}", e)))?;

        // iter_dup_of positions the cursor at the key and returns an Iter over all
        // duplicate values for that key. Item type is (&[u8], &[u8]) = (key, value).
        let key: &[u8] = accession.as_bytes();
        let oids = match cursor.iter_dup_of(&key) {
            Ok(iter) => iter
                .filter_map(|(_, val)| {
                    if val.len() >= 4 {
                        Some(u32::from_le_bytes([val[0], val[1], val[2], val[3]]))
                    } else {
                        None
                    }
                })
                .collect(),
            Err(lmdb::Error::NotFound) => vec![],
            Err(e) => return Err(DbError::InvalidFormat(format!("LMDB lookup: {}", e))),
        };

        drop(cursor);
        txn.abort();
        Ok(oids)
    }

    /// Return (volume_name, num_oids) for each volume stored in the LMDB.
    pub fn get_volumes_info(&self) -> Result<Vec<(String, u32)>> {
        let txn = self.env.begin_ro_txn()
            .map_err(|e| DbError::InvalidFormat(format!("LMDB txn: {}", e)))?;

        // iter_dup_start on a non-DUPSORT database yields one Iter per entry,
        // each containing exactly one item. flat_map(|i| i) flattens them.
        let mut cursor_name = txn.open_ro_cursor(self.db_volname)
            .map_err(|e| DbError::InvalidFormat(format!("LMDB cursor volname: {}", e)))?;
        let mut cursor_info = txn.open_ro_cursor(self.db_volinfo)
            .map_err(|e| DbError::InvalidFormat(format!("LMDB cursor volinfo: {}", e)))?;

        let names: Vec<String> = cursor_name
            .iter_start()
            .map(|(_, val)| String::from_utf8_lossy(val).into_owned())
            .collect();

        let counts: Vec<u32> = cursor_info
            .iter_start()
            .map(|(_, val)| {
                if val.len() >= 4 {
                    u32::from_le_bytes([val[0], val[1], val[2], val[3]])
                } else {
                    0
                }
            })
            .collect();

        drop(cursor_name);
        drop(cursor_info);
        txn.abort();

        Ok(names.into_iter().zip(counts).collect())
    }

    /// Iterate over all (accession, oid) pairs in the acc2oid database.
    /// Handles DUPSORT: calls `f(accession, oid)` for every (key, value) pair including
    /// multiple OIDs per accession.
    pub fn iter_accessions<F>(&self, mut f: F) -> Result<()>
    where
        F: FnMut(&str, u32),
    {
        let txn = self.env.begin_ro_txn()
            .map_err(|e| DbError::InvalidFormat(format!("LMDB txn: {}", e)))?;
        let mut cursor = txn.open_ro_cursor(self.db_acc2oid)
            .map_err(|e| DbError::InvalidFormat(format!("LMDB cursor: {}", e)))?;

        // iter_dup_start yields one Iter<'txn> per distinct key.
        // Each inner Iter yields all (key, val) duplicate pairs for that key.
        for inner_iter in cursor.iter_dup_start() {
            for (key_bytes, val) in inner_iter {
                let acc = std::str::from_utf8(key_bytes).unwrap_or("");
                if val.len() >= 4 {
                    let oid = u32::from_le_bytes([val[0], val[1], val[2], val[3]]);
                    f(acc, oid);
                }
            }
        }

        drop(cursor);
        txn.abort();
        Ok(())
    }
}
