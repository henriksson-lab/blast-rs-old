pub mod matrix;
pub mod stats;
mod tests;
pub mod lookup;
pub mod extend;
pub mod align;
pub mod search;
pub mod hsp;

pub use matrix::{ScoringMatrix, MatrixType};
pub use stats::KarlinAltschul;
pub use hsp::{Hsp, SearchResult};
pub use search::{SearchParams, blast_search};
