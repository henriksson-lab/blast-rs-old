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

## Key design details

- Database files are memory-mapped, not read into heap buffers
- Nucleotide sequences are stored 4 bases per byte (NCBI 2-bit encoding) with ambiguity data in a separate section
- Protein sequences use NCBI stdaa (28-letter) encoding internally
- 2-hit seeding for protein searches, single-hit for nucleotide (configurable)
- The `ncbi-blast-2.17.0+-src/` directory contains the reference NCBI C++ source for comparison; it is not part of the build
