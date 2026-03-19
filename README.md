# blast-rs

A pure-Rust implementation of BLAST (Basic Local Alignment Search Tool). Reads and writes databases created by NCBI `makeblastdb`, implements the core BLAST algorithm for protein and nucleotide search, and produces output in all standard BLAST formats.

No dependency on the NCBI C++ toolkit.

## Features

- Read BLAST databases v4 and v5 (protein and nucleotide)
- Build BLAST databases from FASTA input (`makeblastdb`)
- `blastp` — protein-protein search with BLOSUM/PAM scoring matrices
- `blastn` — nucleotide-nucleotide search
- Parallel search via Rayon
- All standard BLAST output formats (0–18)
- Custom tabular column selection

## Workspace layout

```
blast-rs/
  blast-db/     library: database reader and writer
  blast-core/   library: search algorithm, statistics, scoring matrices
  blast-cli/    binary: blast-cli (blastp, blastn, makeblastdb, dumpdb)
```

## Building

Requires Rust 1.70+ and a system LMDB library (for v5 database support).

```sh
cargo build --release
```

The binary is at `target/release/blast-cli`.

## Usage

### Build a database

```sh
# Protein database
blast-cli makeblastdb -i sequences.faa --dbtype prot -o mydb

# Nucleotide database
blast-cli makeblastdb -i sequences.fna --dbtype nucl -o mydb --title "My genome DB"
```

FASTA headers are parsed as `>accession description` (the first whitespace-delimited word becomes the accession). IUPAC ambiguity codes (N, R, Y, M, K, S, W, H, B, V, D) are preserved in nucleotide databases.

### Protein search (blastp)

```sh
blast-cli blastp -q query.faa -d mydb
```

### Nucleotide search (blastn)

```sh
blast-cli blastn -q query.fna -d mydb
```

### Common options

| Flag | Default | Description |
|------|---------|-------------|
| `-q / --query` | — | Query FASTA file |
| `-d / --db` | — | Database base path (no extension) |
| `-o / --out` | stdout | Output file |
| `--evalue` | 10 | E-value threshold |
| `--outfmt` | 0 | Output format (see below) |
| `--max-target-seqs` | 500 | Maximum hits returned |
| `--num_threads` | 0 (all) | Worker threads |
| `--matrix` | BLOSUM62 | Scoring matrix (`blastp` only) |
| `--gapopen` | matrix default | Gap open penalty |
| `--gapextend` | matrix default | Gap extend penalty |
| `--word-size` | 3 (prot) / 11 (nucl) | Seed word size |
| `--reward` | 2 | Match reward (`blastn` only) |
| `--penalty` | -3 | Mismatch penalty (`blastn` only) |

### Output formats

`--outfmt` accepts a format number, optionally followed by column names for tabular formats:

```sh
blast-cli blastp -q query.faa -d mydb --outfmt 6
blast-cli blastp -q query.faa -d mydb --outfmt "6 qseqid sseqid pident evalue bitscore"
blast-cli blastp -q query.faa -d mydb --outfmt 15   # JSON
blast-cli blastp -q query.faa -d mydb --outfmt 5    # BLAST XML
```

| Format | Description |
|--------|-------------|
| 0 | Pairwise (default) |
| 1 | Query-anchored with identities |
| 2 | Query-anchored, no identities |
| 3 | Flat query-anchored with identities |
| 4 | Flat query-anchored, no identities |
| 5 | BLAST XML |
| 6 | Tabular |
| 7 | Tabular with comment lines |
| 10 | Comma-separated (CSV) |
| 15 | Single-file BLAST JSON |
| 16 | Single-file BLAST XML2 |
| 17 | Subject sequences as FASTA |

Tabular columns available: `qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen nident positive gaps ppos qseq sseq btop staxid salltitles qcovs qcovhsp score`

### Inspect a database

```sh
# Dump all sequences as FASTA
blast-cli dumpdb -d mydb

# Show headers only
blast-cli dumpdb -d mydb --headers-only

# v5 databases: list all accessions
blast-cli dumpdb -d mydb --list-accessions

# v5 databases: look up OIDs for an accession
blast-cli dumpdb -d mydb --lookup NP_001234

# v5 databases: show volume info
blast-cli dumpdb -d mydb --volumes
```

## Database format support

| Feature | Supported |
|---------|-----------|
| v4 protein (`.pin / .psq / .phr`) | read + write |
| v4 nucleotide (`.nin / .nsq / .nhr`) | read + write |
| v5 LMDB accession index (`.pdb / .ndb`) | read |
| v5 OID→SeqIds (`.pos / .nos`) | read |
| v5 OID→TaxIds (`.pot / .not`) | read |
| Multi-volume databases | not yet |

## Scoring matrices

BLOSUM45, BLOSUM62, BLOSUM80, PAM30, PAM70, PAM250.
Karlin-Altschul parameters are pre-computed for standard gap penalty combinations.

## Dependencies

| Crate | Purpose |
|-------|---------|
| `byteorder` | Big/little-endian binary I/O |
| `memmap2` | Memory-mapped database files |
| `lmdb` | v5 LMDB accession index |
| `rayon` | Parallel search |
| `clap` | CLI argument parsing |
| `thiserror` | Error types |

## License

MIT
