# TODO

## Benchmark: blast-rs vs NCBI C++ BLAST

NCBI BLAST 2.17.0 is installed at `/usr/local/Cellar/blast/2.17.0_1/bin/`.

### Setup

1. Generate test data — create protein and nucleotide FASTA files of varying sizes (100, 1000, 10000 sequences)
2. Build databases with both tools:
   ```
   makeblastdb -in test.faa -dbtype prot -out testdb_ncbi
   ./target/release/blast-cli makeblastdb -i test.faa --dbtype prot -o testdb_rs
   ```
3. Pick representative queries (short ~50aa, medium ~300aa, long ~1000aa)

### Benchmarks to run

| Test | NCBI command | blast-rs command |
|------|-------------|-----------------|
| blastp small query vs small db | `blastp -query q.faa -db testdb -outfmt 6` | `blast-cli blastp -q q.faa -db testdb --outfmt 6` |
| blastp large db | same, larger db | same, larger db |
| blastn | `blastn -query q.fna -db testdb` | `blast-cli blastn -q q.fna -db testdb` |
| megablast | `blastn -task megablast ...` | `blast-cli blastn --task megablast ...` |
| blastx | `blastx -query q.fna -db protdb` | `blast-cli blastx -q q.fna -db protdb` |
| makeblastdb | `makeblastdb -in big.faa -dbtype prot -out db` | `blast-cli makeblastdb -i big.faa --dbtype prot -o db` |
| multi-thread scaling | vary `-num_threads 1,2,4,8` | vary `--num_threads 1,2,4,8` |

### Metrics to capture

- Wall-clock time (`time` or `/usr/bin/time -l`)
- Peak memory (RSS)
- Result correctness: compare output (outfmt 6) between both tools — diff hit counts, top hits, e-values

### Correctness validation

Compare tabular output (format 6) from both tools on the same query+db:
```sh
diff <(sort ncbi_results.tsv) <(sort rs_results.tsv)
```
Small e-value differences are expected; focus on same hits found and same ranking.

## Remaining features (~20%)

- deltablast / rpsblast / rpstblastn (need CDD/profile DB infrastructure)
- ASN.1 output formats (8, 9, 11)
- PHI-BLAST (pattern-hit initiated)
- Remote BLAST (NCBI server search)
- Full WindowMasker with taxonomy-based repeat databases
- Full LMDB v5 database write (currently index-only)
- SRA input support
- blastdbcheck utility
