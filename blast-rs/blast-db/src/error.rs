use thiserror::Error;

#[derive(Debug, Error)]
pub enum DbError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Invalid format: {0}")]
    InvalidFormat(String),
    #[error("Unsupported format version: {0}")]
    UnsupportedVersion(i32),
    #[error("OID out of range: {0}")]
    OidOutOfRange(u32),
    #[error("ASN.1 parse error: {0}")]
    AsnParse(String),
}

pub type Result<T> = std::result::Result<T, DbError>;
