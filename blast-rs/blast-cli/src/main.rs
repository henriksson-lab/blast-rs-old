use std::path::{Path, PathBuf};
use std::fs;
use std::io::{self, BufRead, Write};
use clap::{Parser, Subcommand};
use blast_db::BlastDb;
use blast_core::{
    SearchParams, SearchResult, Hsp,
    matrix::MatrixType,
    search::{blast_search, blastn_search},
};

#[derive(Parser)]
#[command(name = "blast-cli", about = "Pure-Rust BLAST implementation")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Protein-protein BLAST search
    Blastp(BlastArgs),
    /// Nucleotide-nucleotide BLAST search
    Blastn(BlastArgs),
    /// Dump sequences from a BLAST database (for testing)
    Dumpdb(DumpArgs),
}

#[derive(clap::Args)]
struct BlastArgs {
    /// Query FASTA file
    #[arg(short = 'q', long)]
    query: PathBuf,
    /// Database path (without extension)
    #[arg(short, long)]
    db: PathBuf,
    /// Output file (default: stdout)
    #[arg(short, long)]
    out: Option<PathBuf>,
    /// E-value threshold
    #[arg(long, default_value = "10")]
    evalue: f64,
    /// Output format: 0=pairwise, 6=tabular
    #[arg(long = "outfmt", default_value = "0")]
    outfmt: u32,
    /// Maximum target sequences
    #[arg(long, default_value = "500")]
    max_target_seqs: usize,
    /// Scoring matrix (blastp only)
    #[arg(long, default_value = "BLOSUM62")]
    matrix: String,
    /// Gap open penalty
    #[arg(long = "gapopen")]
    gap_open: Option<i32>,
    /// Gap extend penalty
    #[arg(long = "gapextend")]
    gap_extend: Option<i32>,
    /// Word size
    #[arg(long)]
    word_size: Option<usize>,
    /// Number of threads (0 = all available)
    #[arg(long = "num_threads", default_value = "0")]
    num_threads: usize,
    /// Match score (blastn only)
    #[arg(long = "reward", default_value = "2")]
    match_score: i32,
    /// Mismatch penalty (blastn only, positive value will be negated)
    #[arg(long = "penalty", default_value = "-3")]
    mismatch: i32,
}

#[derive(clap::Args)]
struct DumpArgs {
    /// Database path (without extension)
    #[arg(short, long)]
    db: PathBuf,
    /// Max sequences to dump (0 = all)
    #[arg(long, default_value = "0")]
    max: usize,
    /// Show headers only (no sequences)
    #[arg(long)]
    headers_only: bool,
    /// List all accessions from the v5 LMDB index
    #[arg(long)]
    list_accessions: bool,
    /// Look up OIDs for a specific accession (v5 only)
    #[arg(long)]
    lookup: Option<String>,
    /// Show v5 volume info from LMDB
    #[arg(long)]
    volumes: bool,
}

fn main() {
    let cli = Cli::parse();
    match &cli.command {
        Commands::Blastp(args) => run_blastp(args),
        Commands::Blastn(args) => run_blastn(args),
        Commands::Dumpdb(args) => run_dumpdb(args),
    }
}

fn run_blastp(args: &BlastArgs) {
    let db = BlastDb::open(&args.db).unwrap_or_else(|e| {
        eprintln!("Error opening database: {}", e);
        std::process::exit(1);
    });

    if args.num_threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.num_threads)
            .build_global()
            .ok();
    }

    let matrix: MatrixType = args.matrix.parse().unwrap_or_else(|e: String| {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    });

    let mut params = SearchParams::blastp_defaults();
    params.matrix = matrix;
    params.evalue_threshold = args.evalue;
    params.max_target_seqs = args.max_target_seqs;
    if let Some(go) = args.gap_open { params.gap_open = go; }
    if let Some(ge) = args.gap_extend { params.gap_extend = ge; }
    if let Some(ws) = args.word_size { params.word_size = ws; }

    let queries = read_fasta(&args.query).unwrap_or_else(|e| {
        eprintln!("Error reading query: {}", e);
        std::process::exit(1);
    });

    let out: Box<dyn Write> = match &args.out {
        Some(p) => Box::new(fs::File::create(p).unwrap()),
        None => Box::new(io::stdout()),
    };
    let mut out = io::BufWriter::new(out);

    for (query_title, query_seq) in &queries {
        // Convert ASCII amino acids to Ncbistdaa
        let query_ncbistdaa = aa_to_ncbistdaa(query_seq);
        let results = blast_search(&db, &query_ncbistdaa, &params);

        match args.outfmt {
            6 => write_tabular(&mut out, query_title, &results),
            _ => write_pairwise(&mut out, query_title, query_seq, &results, &db),
        }
    }
}

fn run_blastn(args: &BlastArgs) {
    let db = BlastDb::open(&args.db).unwrap_or_else(|e| {
        eprintln!("Error opening database: {}", e);
        std::process::exit(1);
    });

    if args.num_threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.num_threads)
            .build_global()
            .ok();
    }

    let mut params = SearchParams::blastn_defaults();
    params.evalue_threshold = args.evalue;
    params.max_target_seqs = args.max_target_seqs;
    params.match_score = args.match_score;
    params.mismatch = args.mismatch;
    if let Some(go) = args.gap_open { params.gap_open = go; }
    if let Some(ge) = args.gap_extend { params.gap_extend = ge; }
    if let Some(ws) = args.word_size { params.word_size = ws; }

    let queries = read_fasta(&args.query).unwrap_or_else(|e| {
        eprintln!("Error reading query: {}", e);
        std::process::exit(1);
    });

    let out: Box<dyn Write> = match &args.out {
        Some(p) => Box::new(fs::File::create(p).unwrap()),
        None => Box::new(io::stdout()),
    };
    let mut out = io::BufWriter::new(out);

    for (query_title, query_seq) in &queries {
        let results = blastn_search(&db, query_seq, &params);

        match args.outfmt {
            6 => write_tabular(&mut out, query_title, &results),
            _ => write_pairwise(&mut out, query_title, query_seq, &results, &db),
        }
    }
}

fn run_dumpdb(args: &DumpArgs) {
    let db = BlastDb::open(&args.db).unwrap_or_else(|e| {
        eprintln!("Error opening database: {}", e);
        std::process::exit(1);
    });

    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());

    // v5: accession lookup
    if let Some(acc) = &args.lookup {
        match db.lookup_accession(acc) {
            None => eprintln!("Database is v4 (no LMDB accession index)."),
            Some(Err(e)) => eprintln!("Lookup error: {}", e),
            Some(Ok(oids)) => {
                if oids.is_empty() {
                    writeln!(out, "Accession '{}' not found.", acc).unwrap();
                } else {
                    for oid in oids {
                        writeln!(out, "{}\t{}", acc, oid).unwrap();
                    }
                }
            }
        }
        return;
    }

    // v5: list volume info
    if args.volumes {
        match db.get_volumes_info() {
            None => eprintln!("Database is v4 (no LMDB volume info)."),
            Some(Err(e)) => eprintln!("Volume info error: {}", e),
            Some(Ok(vols)) => {
                writeln!(out, "DB version: {}", db.format_version()).unwrap();
                writeln!(out, "Volumes ({}):", vols.len()).unwrap();
                for (name, n_oids) in &vols {
                    writeln!(out, "  {} ({} OIDs)", name, n_oids).unwrap();
                }
            }
        }
        return;
    }

    // v5: list all accessions
    if args.list_accessions {
        match db.iter_accessions(|acc, oid| {
            writeln!(out, "{}\t{}", acc, oid).ok();
        }) {
            None => eprintln!("Database is v4 (no LMDB accession index)."),
            Some(Err(e)) => eprintln!("Iteration error: {}", e),
            Some(Ok(())) => {}
        }
        return;
    }

    // Standard sequence dump
    let n = if args.max == 0 { db.num_sequences() } else { args.max as u32 };

    for oid in 0..n {
        let header = db.get_header(oid).unwrap_or_default();

        // For v5, prefer the .pos/.nos seqid list for the defline.
        let accession = if let Some(Ok(ids)) = db.get_seqids(oid) {
            ids.into_iter().next().unwrap_or_else(|| header.accession.clone())
        } else {
            header.accession.clone()
        };

        let title = if header.title.is_empty() {
            format!("OID={}", oid)
        } else {
            header.title.clone()
        };

        if args.headers_only {
            writeln!(out, ">{} {} [taxid={}]", accession, title, header.taxid).unwrap();
            continue;
        }

        match db.seq_type() {
            blast_db::index::SeqType::Protein => {
                let seq = db.get_sequence_protein(oid).unwrap_or_default();
                writeln!(out, ">{} {}", accession, title).unwrap();
                for chunk in seq.chunks(60) {
                    writeln!(out, "{}", std::str::from_utf8(chunk).unwrap_or("?")).unwrap();
                }
            }
            blast_db::index::SeqType::Nucleotide => {
                let seq = db.get_sequence_nucleotide(oid).unwrap_or_default();
                writeln!(out, ">{} {}", accession, title).unwrap();
                for chunk in seq.chunks(60) {
                    writeln!(out, "{}", std::str::from_utf8(chunk).unwrap_or("?")).unwrap();
                }
            }
        }
    }
}

// ---- Output formatters ----

fn write_tabular(out: &mut impl Write, query_title: &str, results: &[SearchResult]) {
    // Standard 12-column tabular (outfmt 6)
    // qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    for r in results {
        for hsp in &r.hsps {
            let qid = query_title.split_whitespace().next().unwrap_or(query_title);
            let sid = if !r.subject_accession.is_empty() { &r.subject_accession } else { &r.subject_title };
            let pident = hsp.percent_identity();
            let length = hsp.alignment_length;
            let mismatches = length - hsp.num_identities - hsp.num_gaps;
            let gapopen = count_gap_opens(&hsp.query_aln) + count_gap_opens(&hsp.subject_aln);
            writeln!(
                out,
                "{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2e}\t{:.1}",
                qid,
                sid,
                pident,
                length,
                mismatches,
                gapopen,
                hsp.query_start + 1,  // 1-based
                hsp.query_end,
                hsp.subject_start + 1,
                hsp.subject_end,
                hsp.evalue,
                hsp.bit_score,
            ).unwrap();
        }
    }
}

fn count_gap_opens(aln: &[u8]) -> usize {
    let mut opens = 0usize;
    let mut in_gap = false;
    for &b in aln {
        if b == b'-' {
            if !in_gap {
                opens += 1;
                in_gap = true;
            }
        } else {
            in_gap = false;
        }
    }
    opens
}

fn write_pairwise(
    out: &mut impl Write,
    query_title: &str,
    _query_seq: &[u8],
    results: &[SearchResult],
    _db: &BlastDb,
) {
    writeln!(out, "Query: {}", query_title).unwrap();
    writeln!(out, "").unwrap();

    if results.is_empty() {
        writeln!(out, "***** No significant similarity found. *****").unwrap();
        writeln!(out, "").unwrap();
        return;
    }

    // Summary section
    writeln!(out, "Sequences producing significant alignments:").unwrap();
    writeln!(out, "").unwrap();
    for r in results {
        let title = if !r.subject_accession.is_empty() {
            format!("{} {}", r.subject_accession, r.subject_title)
        } else {
            r.subject_title.clone()
        };
        let best_hsp = &r.hsps[0];
        writeln!(
            out,
            "{:<60}  {:.1}  {:.2e}",
            &title[..title.len().min(60)],
            best_hsp.bit_score,
            best_hsp.evalue
        ).unwrap();
    }
    writeln!(out, "").unwrap();

    // Alignment section
    for r in results {
        let title = if !r.subject_accession.is_empty() {
            format!("{} {}", r.subject_accession, r.subject_title)
        } else {
            r.subject_title.clone()
        };
        writeln!(out, "> {}", title).unwrap();
        writeln!(out, "   Length = {}", r.subject_len).unwrap();
        writeln!(out, "").unwrap();

        for hsp in &r.hsps {
            writeln!(out, " Score = {:.1} bits ({}),  Expect = {:.2e}",
                hsp.bit_score, hsp.score, hsp.evalue).unwrap();
            writeln!(out, " Identities = {}/{} ({:.0}%), Gaps = {}/{}",
                hsp.num_identities,
                hsp.alignment_length,
                hsp.percent_identity(),
                hsp.num_gaps,
                hsp.alignment_length,
            ).unwrap();
            writeln!(out, "").unwrap();

            // Print alignment in 60-char blocks
            let aln_len = hsp.query_aln.len();
            let mut pos = 0;
            let mut q_pos = hsp.query_start;
            let mut s_pos = hsp.subject_start;
            while pos < aln_len {
                let end = (pos + 60).min(aln_len);
                let q_chunk = &hsp.query_aln[pos..end];
                let m_chunk = &hsp.midline[pos..end];
                let s_chunk = &hsp.subject_aln[pos..end];

                let q_chars: usize = q_chunk.iter().filter(|&&b| b != b'-').count();
                let s_chars: usize = s_chunk.iter().filter(|&&b| b != b'-').count();

                writeln!(out, "Query  {:>5}  {}  {}",
                    q_pos + 1,
                    std::str::from_utf8(q_chunk).unwrap_or("?"),
                    q_pos + q_chars,
                ).unwrap();
                writeln!(out, "             {}",
                    std::str::from_utf8(m_chunk).unwrap_or("?"),
                ).unwrap();
                writeln!(out, "Sbjct  {:>5}  {}  {}",
                    s_pos + 1,
                    std::str::from_utf8(s_chunk).unwrap_or("?"),
                    s_pos + s_chars,
                ).unwrap();
                writeln!(out, "").unwrap();

                q_pos += q_chars;
                s_pos += s_chars;
                pos = end;
            }
        }
    }
}

// ---- FASTA reader ----

fn read_fasta(path: &Path) -> io::Result<Vec<(String, Vec<u8>)>> {
    let file = fs::File::open(path)?;
    let reader = io::BufReader::new(file);
    let mut sequences = Vec::new();
    let mut current_title = String::new();
    let mut current_seq: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim().to_string();
        if line.is_empty() { continue; }
        if line.starts_with('>') {
            if !current_title.is_empty() {
                sequences.push((current_title.clone(), current_seq.clone()));
                current_seq.clear();
            }
            current_title = line[1..].to_string();
        } else {
            current_seq.extend_from_slice(line.as_bytes());
        }
    }
    if !current_title.is_empty() {
        sequences.push((current_title, current_seq));
    }
    Ok(sequences)
}

// ---- Encoding converters ----

/// Convert ASCII amino acid sequence to Ncbistdaa encoding.
/// Inverse of NCBISTDAA_TO_AA table.
fn aa_to_ncbistdaa(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&c| ascii_to_ncbistdaa(c)).collect()
}

/// Convert a single ASCII amino acid character to its Ncbistdaa code.
fn ascii_to_ncbistdaa(c: u8) -> u8 {
    // Ncbistdaa: - A B C D E F G H I K L M N P Q R S T V W X Y Z U * O J
    //             0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
    match c.to_ascii_uppercase() {
        b'A' => 1,
        b'B' => 2,
        b'C' => 3,
        b'D' => 4,
        b'E' => 5,
        b'F' => 6,
        b'G' => 7,
        b'H' => 8,
        b'I' => 9,
        b'K' => 10,
        b'L' => 11,
        b'M' => 12,
        b'N' => 13,
        b'P' => 14,
        b'Q' => 15,
        b'R' => 16,
        b'S' => 17,
        b'T' => 18,
        b'V' => 19,
        b'W' => 20,
        b'X' => 21,
        b'Y' => 22,
        b'Z' => 23,
        b'U' => 24,
        b'*' => 25,
        b'O' => 26,
        b'J' => 27,
        _ => 21, // X = unknown
    }
}
