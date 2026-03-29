# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

blast-rs is a pure-Rust implementation of NCBI BLAST (Basic Local Alignment Search Tool). It reads/writes BLAST databases, implements blastp, blastn, blastx, tblastn, tblastx, and psiblast, and produces all standard BLAST output formats. No dependency on the NCBI C++ toolkit.

## Build and test commands

```sh
cargo build --release        # binary at target/release/blast-cli
cargo test                   # all tests
cargo test test_name         # single test by name
cargo check                  # type-check without building
```

Requires Rust 1.70+ and system LMDB library (for v5 database support).

## Architecture

Single crate (`blast-rs`) with both a library (`src/lib.rs`) and a binary (`src/main.rs`).

- **`src/db/`** — Database I/O layer. Reads v4/v5 BLAST databases (memory-mapped via `memmap2`), writes v4/v5. Handles binary formats for sequences, headers, index, and v5 LMDB accession indexes. Key types: `BlastDb`, `BlastDbBuilder`. Multi-volume and alias file support.

- **`src/`** (library modules) — Search algorithms and scoring. Contains the full BLAST pipeline: word lookup table construction (`lookup.rs`), ungapped/gapped extension (`extend.rs`), Karlin-Altschul statistics (`stats.rs`), all 8 scoring matrices (`matrix.rs`), low-complexity masking SEG/DUST/repeat (`mask.rs`), composition-based statistics modes 0-3 (`compo.rs`), PSSM construction for PSI-BLAST (`pssm.rs`), 6-frame translation with 23 genetic codes (`translate.rs`), and discontiguous megablast (`lookup.rs`). The `api.rs` module provides high-level entry points. Parallelism via Rayon.

- **`src/main.rs`** + **`src/output.rs`** — CLI binary. Subcommands: `blastp`, `blastn`, `blastx`, `tblastn`, `tblastx`, `psiblast`, `makeblastdb`, `dumpdb`, `blastdbcmd`, `blastdb-aliastool`. Output formats: pairwise, tabular, CSV, XML, JSON, XML2, SAM, FASTA, HTML. Uses `clap` derive for argument parsing.

From `main.rs`/`output.rs`, library items are accessed via `blast_rs::` (e.g., `blast_rs::db::BlastDb`). From library modules, database types are at `crate::db::`.

## Search pipeline flow

The core search (e.g. `blast_search()` in `search.rs`) follows this sequence:

1. **Masking** (`mask.rs`) — SEG for protein, DUST for nucleotide. Soft masking masks only during seeding; hard masking replaces with X/N permanently.
2. **Scoring setup** (`matrix.rs`, `stats.rs`) — Load scoring matrix (28×28 Ncbistdaa), look up pre-computed Karlin-Altschul parameters (λ, K, H) for the matrix+gap combination.
3. **Lookup table** (`lookup.rs`) — `ProteinLookup` builds a HashMap of neighboring words (base-28 encoding, scoring ≥ threshold). `NucleotideLookup` scans for exact seed matches. `DiscontiguousLookup` for megablast templates.
4. **Parallel scan** — Rayon `par_iter` over database OIDs. Each thread independently:
   - Retrieves subject sequence from memory-mapped DB
   - Finds seed hits; uses diagonal tracking (`diag = s_pos - q_pos + offset`) for 2-hit detection (protein default window=40)
   - Ungapped extension (`extend.rs`) with X-drop
   - Gapped extension (Smith-Waterman with affine gaps) for top ungapped hits
   - Covered-region tracking prevents duplicate overlapping HSPs
5. **Statistics** (`stats.rs`, `compo.rs`) — E-value = m·n·K·exp(−λ·S), effective length correction for edge effects, optional composition-based λ' adjustment (modes 0–3).
6. **Result filtering** — By E-value threshold and `max_target_seqs`, sorted by best E-value.

Entry points in `api.rs` (`blastp()`, `blastn()`, etc.) configure `SearchParams` with appropriate defaults and call into `search.rs`.

## Key types

- `SearchParams` — Builder-pattern struct with all search settings (~20 fields). Presets: `SearchParams::blastp()`, `::blastn()`, etc.
- `SearchResult` / `Hsp` (`hsp.rs`) — Per-subject results containing HSP list with scores, coordinates, alignment strings.
- `BlastDb` (`db/blastdb.rs`) — Database handle; holds memory-mapped `Volume` instances. `open()` auto-detects v4/v5 and alias files.
- `BlastDbBuilder` (`db/builder.rs`) — Accumulates `SequenceEntry` structs and writes v4/v5 binary format.
- `ScoringMatrix` (`matrix.rs`) — 28×28 table indexed by Ncbistdaa codes (gap=0, A=1, ..., *=25, O=26, J=27).
- `Pssm` (`pssm.rs`) — Position-specific scoring matrix for PSI-BLAST; scores\[position\]\[ncbistdaa_code\].

## Key design details

- Database files are memory-mapped, not read into heap buffers
- Nucleotide sequences are stored 4 bases per byte (NCBI NcbiNa2, MSB first; remainder count in low 2 bits of last byte) with ambiguity data in a separate section using IUPAC codes
- Protein sequences use NCBI stdaa (28-letter) encoding internally; `decode_protein()` in `db/sequence.rs` converts to ASCII
- 2-hit seeding for protein searches, single-hit for nucleotide (configurable)
- Gapped extension uses separate X-drop thresholds for preliminary vs final extension (final is higher, allowing longer alignments)
- PSSM pseudocount weighting uses β=10 with BLOSUM62 background frequencies
- The `ncbi-blast-2.17.0+-src/` directory contains the reference NCBI C++ source for comparison; it is not part of the build
