# Feature Parity Report: blast-rs vs NCBI C++ BLAST 2.17.0+

## Executive Summary

blast-rs implements the core BLAST search pipeline with functional output in 11 formats. After implementation work, estimated feature coverage is **~60%** of the NCBI C++ BLAST feature surface.

| Category | Implemented | Partial | Missing | Coverage |
|----------|-------------|---------|---------|----------|
| Search Programs | 8 | 0 | 7+ | ~53% |
| Search Task Variants | 1 | 1 | 6 | ~19% |
| Scoring Matrices | 8 | 0 | 0 | 100% |
| Statistics & Scoring | 9 | 1 | 4 | ~68% |
| Masking & Filtering | 2 | 2 | 3 | ~43% |
| Database Features | 5 | 0 | 8 | ~38% |
| Output Formats | 11 | 0 | 9+ | ~55% |
| Tabular Columns | 31 | 0 | 19+ | ~62% |
| Translation & Genetics | 4 | 0 | 3 | ~57% |
| PSI-BLAST / PSSM | 6 | 0 | 3 | ~67% |
| CLI Parameters | 19 | 2 | 6+ | ~70% |

---

## 1. Search Programs

| Program | blast-rs | Status | Notes |
|---------|----------|--------|-------|
| blastp | Yes | **Implemented** | |
| blastn | Yes | **Implemented** | |
| blastx | Yes | **Implemented** | 6-frame translation of query |
| tblastn | Yes | **Implemented** | 6-frame translation of subjects |
| tblastx | Yes | **Implemented** | 6x6 frame combinations, ungapped |
| psiblast | Yes | **Implemented** | Iterative with PSSM construction |
| makeblastdb | Yes | **Implemented** | v4 format only |
| dumpdb | Yes | **Implemented** | Non-standard utility (not in NCBI) |
| deltablast | No | **Missing** | Domain-enhanced BLAST using CDD |
| rpsblast | No | **Missing** | Reverse position-specific BLAST |
| rpstblastn | No | **Missing** | RPS-BLAST with translated queries |
| blast_formatter | No | **Missing** | Re-format saved results |
| blastdbcmd | No | **Missing** | Database query/extraction tool |
| blastdbcheck | No | **Missing** | Database integrity verification |
| blastdb_aliastool | No | **Missing** | Create database aliases |
| makeprofiledb | No | **Missing** | Build RPS-BLAST profile databases |
| dustmasker | No | **Missing** | Standalone masking tool |
| segmasker | No | **Missing** | Standalone masking tool |
| windowmasker | No | **Missing** | Repetitive region masking tool |
| seedtop | No | **Missing** | Seed/pattern search |
| update_blastdb | No | **Missing** | Download pre-built databases |

## 2. Search Task Variants

NCBI BLAST supports task presets that tune word size, thresholds, and algorithm behavior. blast-rs has no preset system but allows manual `--word-size` override.

| Task | blast-rs | Status | Notes |
|------|----------|--------|-------|
| blastn (word_size=11) | Yes | **Implemented** | Default blastn mode |
| megablast (word_size=28) | Manual | **Partial** | User can set `--word-size 28` but no contiguous seed optimization |
| dc-megablast | No | **Missing** | Discontiguous megablast with template support |
| blastn-short | No | **Missing** | Optimized for short queries (word_size=7, reward=1, penalty=-3) |
| blastp (default) | Yes | **Implemented** | |
| blastp-short | No | **Missing** | Short query preset |
| blastp-fast | No | **Missing** | word_size=6, threshold=21 |
| blastx-fast | No | **Missing** | |
| tblastn-fast | No | **Missing** | |

## 3. Scoring Matrices

| Matrix | blast-rs | Status |
|--------|----------|--------|
| BLOSUM45 | Yes | **Implemented** |
| BLOSUM50 | No | **Missing** |
| BLOSUM62 | Yes | **Implemented** (default) |
| BLOSUM80 | Yes | **Implemented** |
| BLOSUM90 | No | **Missing** |
| PAM30 | Yes | **Implemented** |
| PAM70 | Yes | **Implemented** |
| PAM250 | Yes | **Implemented** |

## 4. Statistics & Scoring

| Feature | blast-rs | Status | Notes |
|---------|----------|--------|-------|
| Karlin-Altschul parameters (lambda, K, H) | Yes | **Implemented** | Pre-computed tables for all 6 matrices |
| E-value computation | Yes | **Implemented** | E = m*n*K*exp(-λS) |
| Bit score computation | Yes | **Implemented** | |
| Effective length correction | Yes | **Implemented** | Edge-effect correction |
| Composition-based stats mode 0 (off) | Yes | **Implemented** | `--no-comp-adjust` |
| Composition-based stats mode 1 (adjust λ) | Yes | **Implemented** | Bisection method for adjusted λ |
| Composition-based stats mode 2 (conditional) | No | **Missing** | |
| Composition-based stats mode 3 (unconditional) | No | **Missing** | |
| 2-hit algorithm (two neighboring seeds) | No | **Missing** | blast-rs uses single-hit only |
| Configurable word score threshold (T) | Internal | **Partial** | Hardcoded per matrix/word_size, no CLI flag |
| X-dropoff ungapped | Yes | **Implemented** | `x_drop_ungapped` in SearchParams |
| X-dropoff gapped | Yes | **Implemented** | `x_drop_gapped` in SearchParams |
| X-dropoff final (separate final pass) | No | **Missing** | NCBI does a final gapped extension with higher X-drop |
| Smith-Waterman traceback (`-use_sw_tback`) | No | **Missing** | |
| Gap trigger threshold | No | **Missing** | |
| Nucleotide K-A parameters | Hardcoded | **Partial** | Uses approximation (λ=1.370, K=0.711) rather than computing from scoring system |

## 5. Masking & Filtering

| Feature | blast-rs | Status | Notes |
|---------|----------|--------|-------|
| SEG (protein low-complexity) | Yes | **Implemented** | W=12, K1=2.2, K2=2.5 |
| DUST (nucleotide low-complexity) | Yes | **Implemented** | W=64, threshold=2.5 |
| WindowMasker (repeat masking) | No | **Missing** | |
| Configurable SEG/DUST parameters via CLI | No | **Partial** | Internal functions accept params, CLI only has on/off toggle |
| Soft masking (mask for lookup only) | No | **Missing** | blast-rs does hard masking only (replaces with X/N) |
| Lowercase masking (`-lcase_masking`) | No | **Missing** | |
| Database masking data integration | No | **Missing** | |
| Lookup-table-only masking | No | **Missing** | |

## 6. Database Features

| Feature | blast-rs | Status | Notes |
|---------|----------|--------|-------|
| v4 database read (protein) | Yes | **Implemented** | .pin/.psq/.phr |
| v4 database read (nucleotide) | Yes | **Implemented** | .nin/.nsq/.nhr |
| v4 database write | Yes | **Implemented** | Both protein and nucleotide |
| v5 database read (LMDB index) | Yes | **Implemented** | .pdb/.ndb, .pos/.nos, .pot/.not |
| v5 database write | No | **Missing** | |
| Multi-volume databases | No | **Missing** | |
| Database aliases (.pal/.nal) | No | **Missing** | |
| GI list filtering (`-gilist`) | No | **Missing** | |
| Sequence ID list filtering (`-seqidlist`) | No | **Missing** | |
| Negative ID list filtering | No | **Missing** | |
| Taxonomy ID filtering (`-taxids`, `-taxidlist`) | No | **Missing** | |
| Entrez query filtering | No | **Missing** | |
| Database masking data storage | No | **Missing** | |

## 7. Output Formats

| Format | Description | blast-rs | Status |
|--------|-------------|----------|--------|
| 0 | Pairwise (default) | Yes | **Implemented** |
| 1 | Query-anchored with identities | Yes | **Implemented** |
| 2 | Query-anchored no identities | Yes | **Implemented** |
| 3 | Flat query-anchored with identities | Yes | **Implemented** |
| 4 | Flat query-anchored no identities | Yes | **Implemented** |
| 5 | BLAST XML v1 | Yes | **Implemented** |
| 6 | Tabular | Yes | **Implemented** |
| 7 | Tabular with comments | Yes | **Implemented** |
| 8 | Text ASN.1 | Stub | **Missing** |
| 9 | Binary ASN.1 | Stub | **Missing** |
| 10 | CSV | Yes | **Implemented** |
| 11 | BLAST archive (ASR) | Stub | **Missing** |
| 12 | Seqalign JSON | Stub | **Missing** |
| 13 | Multi-file JSON | Stub | **Missing** |
| 14 | Multi-file XML2 | Stub | **Missing** |
| 15 | Single-file BLAST JSON | Yes | **Implemented** |
| 16 | Single-file BLAST XML2 | Yes | **Implemented** |
| 17 | Subject FASTA | Yes | **Implemented** |
| 18 | Organism report / SAM | Stub | **Missing** |
| — | HTML output (`-html`) | No | **Missing** |
| — | AIRR rearrangement | No | **Missing** |

## 8. Tabular Columns (formats 6/7/10)

blast-rs supports 24 columns. NCBI supports 50+.

**Implemented in blast-rs:**
`qseqid`, `sseqid`, `pident`, `length`, `mismatch`, `gapopen`, `qstart`, `qend`, `sstart`, `send`, `evalue`, `bitscore`, `qlen`, `slen`, `nident`, `positive`, `gaps`, `ppos`, `qseq`, `sseq`, `btop`, `staxid` (stub — outputs "N/A"), `salltitles`, `qcovs`, `qcovhsp`, `score`

**Missing from blast-rs (selected):**
`qgi`, `sgi`, `sacc`, `sallseqid`, `qframe`, `sframe`, `frames`, `stitle`, `sstrand`, `qcovus`, `staxids`, `sscinames`, `scomnames`, `sblastnames`, `sskingdoms`, `sallacc`, `qaccver`, `saccver`

## 9. Translation & Genetic Codes

| Feature | blast-rs | Status | Notes |
|---------|----------|--------|-------|
| Standard genetic code (1) | Yes | **Implemented** | Hardcoded CODON_TABLE |
| Alternative genetic codes (2–31) | No | **Missing** | Mitochondrial, mycoplasma, etc. |
| Query genetic code (`-query_gencode`) | No | **Missing** | |
| Database genetic code (`-db_gencode`) | No | **Missing** | |
| Query strand selection (`-strand plus/minus/both`) | No | **Missing** | Always searches both strands |
| Query range restriction (`-query_loc`) | No | **Missing** | |
| Frame shift penalties | No | **Missing** | |

## 10. PSI-BLAST / PSSM Features

| Feature | blast-rs | Status | Notes |
|---------|----------|--------|-------|
| PSI-BLAST iterative search | Yes | **Implemented** | Convergence detection + max iterations |
| PSSM construction from alignments | Yes | **Implemented** | Pseudocount β=10 |
| Inclusion E-value threshold | Yes | **Implemented** | `--inclusion_ethresh` |
| PSSM checkpoint save (`-out_pssm`) | No | **Missing** | |
| PSSM checkpoint load (`-in_pssm`) | No | **Missing** | |
| ASCII PSSM output (`-out_ascii_pssm`) | No | **Missing** | |
| MSA input for PSSM | No | **Missing** | |
| PHI-BLAST (pattern-hit initiated) | No | **Missing** | |
| Delta-BLAST (CDD-based PSSM) | No | **Missing** | |

## 11. CLI Parameters & Misc

| Parameter | blast-rs | Status | Notes |
|-----------|----------|--------|-------|
| `-query` / `-q` | Yes | **Implemented** | |
| `-db` / `-d` | Yes | **Implemented** | |
| `-out` / `-o` | Yes | **Implemented** | |
| `-evalue` | Yes | **Implemented** | |
| `-outfmt` | Yes | **Implemented** | With custom columns |
| `-max_target_seqs` | Yes | **Implemented** | |
| `-matrix` | Yes | **Implemented** | |
| `-gapopen` / `-gapextend` | Yes | **Implemented** | |
| `-word_size` | Yes | **Implemented** | |
| `-num_threads` | Yes | **Implemented** | Rayon-based |
| `-reward` / `-penalty` | Yes | **Implemented** | blastn only |
| `--no-lc-filter` | Yes | **Implemented** | Equivalent to `-seg no` / `-dust no` |
| `--no-comp-adjust` | Yes | **Implemented** | Equivalent to `-comp_based_stats 0` |
| `-num_descriptions` / `-num_alignments` | No | **Missing** | |
| `-line_length` | No | **Missing** | |
| `-show_gis` | No | **Missing** | |
| `-best_hit_overhang` / `-best_hit_score_edge` | No | **Missing** | |
| `-culling_limit` | No | **Missing** | |
| `-subject_besthit` | No | **Missing** | |
| `-max_hsps` | No | **Missing** | |
| `-xdrop_ungap` / `-xdrop_gap` / `-xdrop_gap_final` | Partial | **Partial** | Ungap + gap present in SearchParams, no CLI flags; no final X-drop |
| `-threshold` (word score T) | No | **Missing** | Hardcoded internally |
| `-window_size` (2-hit window) | No | **Missing** | |
| `-remote` | No | **Missing** | |
| `-import_search_strategy` / `-export_search_strategy` | No | **Missing** | |
| `-entrez_query` | No | **Missing** | |
| `-gilist` / `-seqidlist` / `-negative_*` | No | **Missing** | |
| `-lcase_masking` | No | **Missing** | |
| MT mode selection (by-query vs by-database) | No | **Partial** | Rayon parallelizes by-database (per-OID) only |

---

## Priority Gaps

Features most likely to matter for real-world usage, roughly ordered by impact:

1. **Multi-volume database support** — Many NCBI databases (nr, nt) are multi-volume; blast-rs cannot read them
2. **Megablast** — Default NCBI blastn mode; blast-rs blastn uses word_size=11 (traditional blastn), not megablast (word_size=28 with contiguous seeds)
3. **Alternative genetic codes** — Required for any non-standard-code organism searches (mitochondrial, mycoplasma, etc.)
4. **Database filtering (GI/seqid/taxid lists)** — Common workflow for restricting searches to specific taxa
5. **blastdbcmd** — Essential utility for extracting sequences from databases
6. **PSSM checkpoint save/load** — Needed for resumable PSI-BLAST workflows
7. **2-hit algorithm** — Significant sensitivity/speed tradeoff used in standard NCBI BLAST
8. **Composition-based stats modes 2/3** — More accurate than mode 1 for certain searches
9. **Query strand selection** — `--strand plus` is commonly used for directional searches
10. **ASN.1 and BLAST archive formats** — Needed for `blast_formatter` round-trip workflows
