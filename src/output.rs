//! All BLAST output format writers.
//!
//! Formats implemented:
//!   0  – Pairwise
//!   1  – Query-anchored with identities
//!   2  – Query-anchored no identities (positive marks)
//!   3  – Flat query-anchored with identities
//!   4  – Flat query-anchored no identities
//!   5  – BLAST XML (version 1)
//!   6  – Tabular (12 standard columns, or user-specified)
//!   7  – Tabular with comment lines
//!   8  – Text ASN.1 (Seq-annot with Dense-seg alignments)
//!   9  – Binary ASN.1 (BER-encoded Seq-annot)
//!  10  – Comma-separated values (CSV)
//!  11  – BLAST archive (Blast4-archive with search params + results)
//!  12  – Seqalign JSON
//!  13  – Multiple-file JSON (delegates to format 15)
//!  14  – Multiple-file XML2 (delegates to format 16)
//!  15  – Single-file BLAST JSON
//!  16  – Single-file BLAST XML2
//!  17  – Subject sequences (FASTA)
//!  18  – Organism report (stub)

use std::io::{self, Write};
use blast_rs::hsp::{Hsp, SearchResult};
use blast_rs::db::BlastDb;

// ─── Public API ────────────────────────────────────────────────────────────

/// Parsed `--outfmt` argument.
#[derive(Debug, Clone)]
pub struct OutputFormat {
    pub fmt_id: u32,
    /// For tabular formats 6/7/10: selected columns (empty = default 12).
    pub columns: Vec<TabularColumn>,
}

impl OutputFormat {
    /// Parse `"6"` or `"6 qseqid sseqid pident"` style strings.
    pub fn parse(s: &str) -> Result<Self, String> {
        let mut tokens = s.split_whitespace();
        let fmt_id: u32 = tokens
            .next()
            .unwrap_or("0")
            .parse()
            .map_err(|_| format!("invalid outfmt number: '{}'", s))?;

        let mut columns: Vec<TabularColumn> = tokens
            .map(|tok| TabularColumn::parse(tok).ok_or_else(|| format!("unknown column: '{}'", tok)))
            .collect::<Result<Vec<_>, _>>()?;

        if columns.is_empty() && matches!(fmt_id, 6 | 7 | 10) {
            columns = TabularColumn::default_columns();
        }

        Ok(OutputFormat { fmt_id, columns })
    }
}

/// Context passed to every formatter.
pub struct SearchContext<'a> {
    pub program: &'a str,
    pub db_path: &'a str,
    pub db_title: &'a str,
    pub db_num_seqs: u64,
    pub db_len: u64,
    pub query_title: &'a str,
    #[allow(dead_code)]
    pub query_seq: &'a [u8],
    pub query_len: usize,
    pub matrix: &'a str,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub evalue_threshold: f64,
    pub iter_num: usize, // 1-based, for multi-query XML/JSON
    /// Maximum number of one-line descriptions to show (None = all)
    pub num_descriptions: Option<usize>,
    /// Maximum number of alignments to show (None = all)
    pub num_alignments: Option<usize>,
}

/// Write all results for one query in the requested format.
pub fn write_results(
    out: &mut impl Write,
    fmt: &OutputFormat,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
    db: Option<&BlastDb>,
) -> io::Result<()> {
    match fmt.fmt_id {
        0 => fmt0_pairwise(out, ctx, results, false, false),
        1 => fmt0_pairwise(out, ctx, results, true, false),
        2 => fmt0_pairwise(out, ctx, results, true, true),
        3 => fmt34_flat(out, ctx, results, false),
        4 => fmt34_flat(out, ctx, results, true),
        5 => fmt5_xml(out, ctx, results),
        6 | 7 | 10 => fmt6_tabular(out, ctx, results, fmt),
        8 => fmt8_asn_text(out, ctx, results),
        9 => fmt9_asn_binary(out, ctx, results),
        11 => fmt11_archive(out, ctx, results),
        12 => fmt12_seqalign_json(out, ctx, results),
        13 => fmt13_json_multi(out, ctx, results),
        14 => fmt14_xml2_multi(out, ctx, results),
        18 => fmt18_sam(out, ctx, results),
        15 => fmt15_json(out, ctx, results),
        16 => fmt16_xml2(out, ctx, results),
        17 => fmt17_fasta(out, ctx, results, db),
        19 => write_results_html(out, ctx, results),
        _ => writeln!(out, "# Unknown output format {}.", fmt.fmt_id),
    }
}

/// Write the XML document-open boilerplate (call once before first query).
pub fn write_xml_header(out: &mut impl Write, program: &str, db: &str, version: &str) -> io::Result<()> {
    writeln!(out, r#"<?xml version="1.0"?>"#)?;
    writeln!(out, r#"<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">"#)?;
    writeln!(out, "<BlastOutput>")?;
    writeln!(out, "  <BlastOutput_program>{}</BlastOutput_program>", xml_escape(program))?;
    writeln!(out, "  <BlastOutput_version>{}</BlastOutput_version>", xml_escape(version))?;
    writeln!(out, "  <BlastOutput_db>{}</BlastOutput_db>", xml_escape(db))?;
    writeln!(out, "  <BlastOutput_iterations>")
}

pub fn write_xml_footer(out: &mut impl Write) -> io::Result<()> {
    writeln!(out, "  </BlastOutput_iterations>")?;
    writeln!(out, "</BlastOutput>")
}

/// Write the JSON array open `[` (call once).
pub fn write_json_header(out: &mut impl Write) -> io::Result<()> {
    writeln!(out, "[")
}
pub fn write_json_footer(out: &mut impl Write) -> io::Result<()> {
    writeln!(out, "]")
}

// ─── Format 0/1/2: Pairwise ────────────────────────────────────────────────

/// Format 0: classic pairwise.
/// Format 1: midline shows '|' for identical, ' ' for mismatch/conservative (no '+').
/// Format 2: midline shows '+' for identical OR conservative, ' ' for mismatch.
fn fmt0_pairwise(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
    query_anchored: bool,    // true → formats 1/2
    positive_marks: bool,    // true → format 2 (use '+' for all conserved, no '|')
) -> io::Result<()> {
    let _qid = first_word(ctx.query_title);

    writeln!(out, "BLAST{} {}  Reference: Rust BLAST implementation", if query_anchored {"X"} else {""}, ctx.program.to_uppercase())?;
    writeln!(out)?;
    writeln!(out, "Database: {}", ctx.db_title)?;
    writeln!(out, "           {} sequences; {} total letters", ctx.db_num_seqs, ctx.db_len)?;
    writeln!(out)?;
    writeln!(out, "Query= {}", ctx.query_title)?;
    writeln!(out, "         ({} letters)", ctx.query_len)?;
    writeln!(out)?;

    if results.is_empty() {
        writeln!(out, "                                                                 Score     E")?;
        writeln!(out, "Sequences producing significant alignments:                      (Bits)  Value")?;
        writeln!(out)?;
        writeln!(out, "***** No significant similarity found. *****")?;
        writeln!(out)?;
        return Ok(());
    }

    writeln!(out, "                                                                 Score     E")?;
    writeln!(out, "Sequences producing significant alignments:                      (Bits)  Value")?;
    writeln!(out)?;
    let desc_limit = ctx.num_descriptions.unwrap_or(results.len());
    for r in results.iter().take(desc_limit) {
        let title = subject_title(r);
        let best = &r.hsps[0];
        writeln!(out, "{:<67} {:.0}  {:.2e}", &title[..title.len().min(67)], best.bit_score, best.evalue)?;
    }
    writeln!(out)?;

    let aln_limit = ctx.num_alignments.unwrap_or(results.len());
    writeln!(out, "ALIGNMENTS")?;
    for r in results.iter().take(aln_limit) {
        let title = subject_title(r);
        writeln!(out, "> {}", title)?;
        writeln!(out, "   Length = {}", r.subject_len)?;
        writeln!(out)?;

        for (hi, hsp) in r.hsps.iter().enumerate() {
            write_hsp_header(out, hsp, hi + 1)?;
            write_alignment_blocks(out, hsp, query_anchored, positive_marks, 60)?;
        }
    }
    Ok(())
}

// ─── Format 3/4: Flat query-anchored ───────────────────────────────────────

fn fmt34_flat(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
    positive_marks: bool,
) -> io::Result<()> {
    let _qid = first_word(ctx.query_title);
    writeln!(out, "Query= {} ({} letters)", ctx.query_title, ctx.query_len)?;
    writeln!(out)?;

    if results.is_empty() {
        writeln!(out, "***** No significant similarity found. *****")?;
        return Ok(());
    }

    for r in results {
        let title = subject_title(r);
        writeln!(out, "> {}", title)?;
        writeln!(out, "   Length = {}", r.subject_len)?;
        for (hi, hsp) in r.hsps.iter().enumerate() {
            write_hsp_header(out, hsp, hi + 1)?;
            write_alignment_blocks(out, hsp, true, positive_marks, 60)?;
        }
    }
    Ok(())
}

// ─── Format 5: BLAST XML 1 ─────────────────────────────────────────────────

fn fmt5_xml(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
) -> io::Result<()> {
    writeln!(out, "    <Iteration>")?;
    writeln!(out, "      <Iteration_iter-num>{}</Iteration_iter-num>", ctx.iter_num)?;
    writeln!(out, "      <Iteration_query-ID>Query_{}</Iteration_query-ID>", ctx.iter_num)?;
    writeln!(out, "      <Iteration_query-def>{}</Iteration_query-def>", xml_escape(ctx.query_title))?;
    writeln!(out, "      <Iteration_query-len>{}</Iteration_query-len>", ctx.query_len)?;
    writeln!(out, "      <Iteration_hits>")?;

    for (hi, r) in results.iter().enumerate() {
        let sid = if !r.subject_accession.is_empty() { &r.subject_accession } else { &r.subject_title };
        writeln!(out, "        <Hit>")?;
        writeln!(out, "          <Hit_num>{}</Hit_num>", hi + 1)?;
        writeln!(out, "          <Hit_id>{}</Hit_id>", xml_escape(sid))?;
        writeln!(out, "          <Hit_def>{}</Hit_def>", xml_escape(&r.subject_title))?;
        writeln!(out, "          <Hit_len>{}</Hit_len>", r.subject_len)?;
        writeln!(out, "          <Hit_hsps>")?;

        for (hj, hsp) in r.hsps.iter().enumerate() {
            let positives = count_positive_chars(&hsp.midline);
            writeln!(out, "            <Hsp>")?;
            writeln!(out, "              <Hsp_num>{}</Hsp_num>", hj + 1)?;
            writeln!(out, "              <Hsp_bit-score>{:.4}</Hsp_bit-score>", hsp.bit_score)?;
            writeln!(out, "              <Hsp_score>{}</Hsp_score>", hsp.score)?;
            writeln!(out, "              <Hsp_evalue>{:.6e}</Hsp_evalue>", hsp.evalue)?;
            writeln!(out, "              <Hsp_query-from>{}</Hsp_query-from>", hsp.query_start + 1)?;
            writeln!(out, "              <Hsp_query-to>{}</Hsp_query-to>", hsp.query_end)?;
            writeln!(out, "              <Hsp_hit-from>{}</Hsp_hit-from>", hsp.subject_start + 1)?;
            writeln!(out, "              <Hsp_hit-to>{}</Hsp_hit-to>", hsp.subject_end)?;
            writeln!(out, "              <Hsp_query-frame>0</Hsp_query-frame>")?;
            writeln!(out, "              <Hsp_hit-frame>0</Hsp_hit-frame>")?;
            writeln!(out, "              <Hsp_identity>{}</Hsp_identity>", hsp.num_identities)?;
            writeln!(out, "              <Hsp_positive>{}</Hsp_positive>", positives)?;
            writeln!(out, "              <Hsp_gaps>{}</Hsp_gaps>", hsp.num_gaps)?;
            writeln!(out, "              <Hsp_align-len>{}</Hsp_align-len>", hsp.alignment_length)?;
            writeln!(out, "              <Hsp_qseq>{}</Hsp_qseq>", xml_escape(str_from_aln(&hsp.query_aln)))?;
            writeln!(out, "              <Hsp_hseq>{}</Hsp_hseq>", xml_escape(str_from_aln(&hsp.subject_aln)))?;
            writeln!(out, "              <Hsp_midline>{}</Hsp_midline>", xml_escape(str_from_aln(&hsp.midline)))?;
            writeln!(out, "            </Hsp>")?;
        }

        writeln!(out, "          </Hit_hsps>")?;
        writeln!(out, "        </Hit>")?;
    }

    writeln!(out, "      </Iteration_hits>")?;
    writeln!(out, "      <Iteration_stat>")?;
    writeln!(out, "        <Statistics>")?;
    writeln!(out, "          <Statistics_db-num>{}</Statistics_db-num>", ctx.db_num_seqs)?;
    writeln!(out, "          <Statistics_db-len>{}</Statistics_db-len>", ctx.db_len)?;
    writeln!(out, "          <Statistics_hsp-len>0</Statistics_hsp-len>")?;
    writeln!(out, "          <Statistics_eff-space>0</Statistics_eff-space>")?;
    writeln!(out, "          <Statistics_kappa>0</Statistics_kappa>")?;
    writeln!(out, "          <Statistics_lambda>0</Statistics_lambda>")?;
    writeln!(out, "          <Statistics_entropy>0</Statistics_entropy>")?;
    writeln!(out, "        </Statistics>")?;
    writeln!(out, "      </Iteration_stat>")?;
    writeln!(out, "    </Iteration>")
}

// ─── Formats 6/7/10: Tabular / CSV ─────────────────────────────────────────

/// All recognised tabular column identifiers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TabularColumn {
    Qseqid,
    Sseqid,
    Pident,
    Length,
    Mismatch,
    Gapopen,
    Qstart,
    Qend,
    Sstart,
    Send,
    Evalue,
    Bitscore,
    // Extended columns
    Qlen,
    Slen,
    Nident,
    Positive,
    Gaps,
    Ppos,
    Qseq,
    Sseq,
    Btop,
    Staxid,
    Salltitles,
    Qcovs,
    QcovHsp,
    Score,
    // Additional columns
    Qframe,
    Sframe,
    Frames,
    Stitle,
    Sacc,
    Sstrand,
    Qcovus,
}

impl TabularColumn {
    fn parse(s: &str) -> Option<Self> {
        Some(match s {
            "qseqid"     => Self::Qseqid,
            "sseqid"     => Self::Sseqid,
            "pident"     => Self::Pident,
            "length"     => Self::Length,
            "mismatch"   => Self::Mismatch,
            "gapopen"    => Self::Gapopen,
            "qstart"     => Self::Qstart,
            "qend"       => Self::Qend,
            "sstart"     => Self::Sstart,
            "send"       => Self::Send,
            "evalue"     => Self::Evalue,
            "bitscore"   => Self::Bitscore,
            "qlen"       => Self::Qlen,
            "slen"       => Self::Slen,
            "nident"     => Self::Nident,
            "positive"   => Self::Positive,
            "gaps"       => Self::Gaps,
            "ppos"       => Self::Ppos,
            "qseq"       => Self::Qseq,
            "sseq"       => Self::Sseq,
            "btop"       => Self::Btop,
            "staxid"     => Self::Staxid,
            "salltitles" => Self::Salltitles,
            "qcovs"      => Self::Qcovs,
            "qcovhsp"    => Self::QcovHsp,
            "score"      => Self::Score,
            "qframe"     => Self::Qframe,
            "sframe"     => Self::Sframe,
            "frames"     => Self::Frames,
            "stitle"     => Self::Stitle,
            "sacc"       => Self::Sacc,
            "sstrand"    => Self::Sstrand,
            "qcovus"     => Self::Qcovus,
            _ => return None,
        })
    }

    fn default_columns() -> Vec<Self> {
        vec![
            Self::Qseqid, Self::Sseqid, Self::Pident, Self::Length,
            Self::Mismatch, Self::Gapopen,
            Self::Qstart, Self::Qend, Self::Sstart, Self::Send,
            Self::Evalue, Self::Bitscore,
        ]
    }

    fn header_name(&self) -> &'static str {
        match self {
            Self::Qseqid     => "qseqid",
            Self::Sseqid     => "sseqid",
            Self::Pident     => "pident",
            Self::Length     => "length",
            Self::Mismatch   => "mismatch",
            Self::Gapopen    => "gapopen",
            Self::Qstart     => "qstart",
            Self::Qend       => "qend",
            Self::Sstart     => "sstart",
            Self::Send       => "send",
            Self::Evalue     => "evalue",
            Self::Bitscore   => "bitscore",
            Self::Qlen       => "qlen",
            Self::Slen       => "slen",
            Self::Nident     => "nident",
            Self::Positive   => "positive",
            Self::Gaps       => "gaps",
            Self::Ppos       => "ppos",
            Self::Qseq       => "qseq",
            Self::Sseq       => "sseq",
            Self::Btop       => "btop",
            Self::Staxid     => "staxid",
            Self::Salltitles => "salltitles",
            Self::Qcovs      => "qcovs",
            Self::QcovHsp    => "qcovhsp",
            Self::Score      => "score",
            Self::Qframe     => "qframe",
            Self::Sframe     => "sframe",
            Self::Frames     => "frames",
            Self::Stitle     => "stitle",
            Self::Sacc       => "sacc",
            Self::Sstrand    => "sstrand",
            Self::Qcovus     => "qcovus",
        }
    }
}

fn fmt6_tabular(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
    fmt: &OutputFormat,
) -> io::Result<()> {
    let sep = if fmt.fmt_id == 10 { "," } else { "\t" };
    let with_comments = fmt.fmt_id == 7;

    if with_comments {
        writeln!(out, "# BLAST{}", ctx.program.to_uppercase())?;
        writeln!(out, "# Query: {}", ctx.query_title)?;
        writeln!(out, "# Database: {}", ctx.db_path)?;
        let col_names: Vec<&str> = fmt.columns.iter().map(|c| c.header_name()).collect();
        writeln!(out, "# Fields: {}", col_names.join(", "))?;
        writeln!(out, "# {} hits found", results.iter().map(|r| r.hsps.len()).sum::<usize>())?;
    }

    let qid = first_word(ctx.query_title);

    for r in results {
        let sid = if !r.subject_accession.is_empty() { &r.subject_accession } else { &r.subject_title };
        for hsp in &r.hsps {
            let fields: Vec<String> = fmt.columns.iter().map(|col| {
                tabular_field(col, qid, sid, hsp, r, ctx)
            }).collect();
            writeln!(out, "{}", fields.join(sep))?;
        }
    }

    if with_comments {
        writeln!(out, "# BLAST processed {} queries", 1)?;
    }

    Ok(())
}

fn tabular_field(
    col: &TabularColumn,
    qid: &str,
    sid: &str,
    hsp: &Hsp,
    r: &SearchResult,
    ctx: &SearchContext<'_>,
) -> String {
    let aln_len = hsp.alignment_length;
    let mismatches = aln_len.saturating_sub(hsp.num_identities + hsp.num_gaps);
    let gap_opens = count_gap_opens(&hsp.query_aln) + count_gap_opens(&hsp.subject_aln);
    let positives = count_positive_chars(&hsp.midline);

    match col {
        TabularColumn::Qseqid     => qid.to_string(),
        TabularColumn::Sseqid     => sid.to_string(),
        TabularColumn::Pident     => format!("{:.2}", hsp.percent_identity()),
        TabularColumn::Length     => aln_len.to_string(),
        TabularColumn::Mismatch   => mismatches.to_string(),
        TabularColumn::Gapopen    => gap_opens.to_string(),
        TabularColumn::Qstart     => (hsp.query_start + 1).to_string(),
        TabularColumn::Qend       => hsp.query_end.to_string(),
        TabularColumn::Sstart     => (hsp.subject_start + 1).to_string(),
        TabularColumn::Send       => hsp.subject_end.to_string(),
        TabularColumn::Evalue     => format!("{:.2e}", hsp.evalue),
        TabularColumn::Bitscore   => format!("{:.1}", hsp.bit_score),
        TabularColumn::Score      => hsp.score.to_string(),
        TabularColumn::Qlen       => ctx.query_len.to_string(),
        TabularColumn::Slen       => r.subject_len.to_string(),
        TabularColumn::Nident     => hsp.num_identities.to_string(),
        TabularColumn::Positive   => positives.to_string(),
        TabularColumn::Gaps       => hsp.num_gaps.to_string(),
        TabularColumn::Ppos       => {
            if aln_len == 0 { "0.00".to_string() }
            else { format!("{:.2}", 100.0 * positives as f64 / aln_len as f64) }
        }
        TabularColumn::Qseq       => str_from_aln(&hsp.query_aln).to_string(),
        TabularColumn::Sseq       => str_from_aln(&hsp.subject_aln).to_string(),
        TabularColumn::Btop       => build_btop(&hsp.query_aln, &hsp.subject_aln),
        TabularColumn::Staxid     => "N/A".to_string(),
        TabularColumn::Salltitles => r.subject_title.clone(),
        TabularColumn::Qcovs      => {
            if ctx.query_len == 0 { "0".to_string() }
            else {
                let cov_len = hsp.query_end - hsp.query_start;
                format!("{}", (100 * cov_len / ctx.query_len).min(100))
            }
        }
        TabularColumn::QcovHsp    => {
            if ctx.query_len == 0 { "0".to_string() }
            else {
                let cov_len = hsp.query_end - hsp.query_start;
                format!("{}", (100 * cov_len / ctx.query_len).min(100))
            }
        }
        TabularColumn::Qframe     => hsp.query_frame.to_string(),
        TabularColumn::Sframe     => hsp.subject_frame.to_string(),
        TabularColumn::Frames     => format!("{}/{}", hsp.query_frame, hsp.subject_frame),
        TabularColumn::Stitle     => r.subject_title.clone(),
        TabularColumn::Sacc       => r.subject_accession.clone(),
        TabularColumn::Sstrand    => {
            if hsp.subject_frame < 0 { "minus".to_string() }
            else { "plus".to_string() }
        }
        TabularColumn::Qcovus     => {
            if ctx.query_len == 0 { "0".to_string() }
            else {
                let cov_len = hsp.query_end - hsp.query_start;
                format!("{}", (100 * cov_len / ctx.query_len).min(100))
            }
        }
    }
}

// ─── Format 8: Text ASN.1 (Seq-annot) ─────────────────────────────────────

/// Build Dense-seg segment data from gapped alignment strings.
/// Returns (starts, lens, numseg) where starts is interleaved [q_start, s_start] per segment.
/// A start of -1 means a gap in that sequence for that segment.
fn build_denseg(hsp: &Hsp) -> (Vec<i64>, Vec<usize>, usize) {
    let mut starts: Vec<i64> = Vec::new();
    let mut lens: Vec<usize> = Vec::new();

    let mut q_pos = hsp.query_start;
    let mut s_pos = hsp.subject_start;
    let mut seg_len = 0usize;
    let mut seg_type: u8 = 0; // 0=match, 1=gap-in-query, 2=gap-in-subject
    let mut seg_q_start: i64 = q_pos as i64;
    let mut seg_s_start: i64 = s_pos as i64;

    for i in 0..hsp.query_aln.len() {
        let q = hsp.query_aln[i];
        let s = hsp.subject_aln[i];
        let cur_type = if q == b'-' { 1u8 } else if s == b'-' { 2u8 } else { 0u8 };

        if i == 0 {
            seg_type = cur_type;
            seg_q_start = if cur_type == 1 { -1 } else { q_pos as i64 };
            seg_s_start = if cur_type == 2 { -1 } else { s_pos as i64 };
            seg_len = 1;
        } else if cur_type != seg_type {
            // flush previous segment
            starts.push(seg_q_start);
            starts.push(seg_s_start);
            lens.push(seg_len);
            // start new segment
            seg_type = cur_type;
            seg_q_start = if cur_type == 1 { -1 } else { q_pos as i64 };
            seg_s_start = if cur_type == 2 { -1 } else { s_pos as i64 };
            seg_len = 1;
        } else {
            seg_len += 1;
        }

        if q != b'-' { q_pos += 1; }
        if s != b'-' { s_pos += 1; }
    }

    // flush last segment
    if seg_len > 0 {
        starts.push(seg_q_start);
        starts.push(seg_s_start);
        lens.push(seg_len);
    }

    let numseg = lens.len();
    (starts, lens, numseg)
}

/// Write a single Seq-align in text ASN.1 notation.
fn write_seqalign_text(
    out: &mut impl Write,
    indent: &str,
    ctx: &SearchContext<'_>,
    r: &SearchResult,
    hsp: &Hsp,
) -> io::Result<()> {
    let (starts, lens, numseg) = build_denseg(hsp);
    let qid = first_word(ctx.query_title);
    let sid = if !r.subject_accession.is_empty() { &r.subject_accession } else { first_word(&r.subject_title) };

    writeln!(out, "{}{{", indent)?;
    writeln!(out, "{}  type partial ,", indent)?;
    writeln!(out, "{}  dim 2 ,", indent)?;
    writeln!(out, "{}  score {{", indent)?;
    writeln!(out, "{}    {{", indent)?;
    writeln!(out, "{}      id str \"score\" ,", indent)?;
    writeln!(out, "{}      value int {}", indent, hsp.score)?;
    writeln!(out, "{}    }} ,", indent)?;
    writeln!(out, "{}    {{", indent)?;
    writeln!(out, "{}      id str \"e_value\" ,", indent)?;
    write_asn_real(out, &format!("{}      ", indent), hsp.evalue)?;
    writeln!(out, "{}    }} ,", indent)?;
    writeln!(out, "{}    {{", indent)?;
    writeln!(out, "{}      id str \"bit_score\" ,", indent)?;
    write_asn_real(out, &format!("{}      ", indent), hsp.bit_score)?;
    writeln!(out, "{}    }} ,", indent)?;
    writeln!(out, "{}    {{", indent)?;
    writeln!(out, "{}      id str \"num_ident\" ,", indent)?;
    writeln!(out, "{}      value int {}", indent, hsp.num_identities)?;
    writeln!(out, "{}    }}", indent)?;
    writeln!(out, "{}  }} ,", indent)?;
    writeln!(out, "{}  segs denseg {{", indent)?;
    writeln!(out, "{}    dim 2 ,", indent)?;
    writeln!(out, "{}    numseg {} ,", indent, numseg)?;
    writeln!(out, "{}    ids {{", indent)?;
    writeln!(out, "{}      local str \"{}\" ,", indent, asn_escape(qid))?;
    writeln!(out, "{}      local str \"{}\"", indent, asn_escape(sid))?;
    writeln!(out, "{}    }} ,", indent)?;
    // starts
    writeln!(out, "{}    starts {{", indent)?;
    for (i, &s) in starts.iter().enumerate() {
        let comma = if i + 1 < starts.len() { " ," } else { "" };
        writeln!(out, "{}      {}{}", indent, s, comma)?;
    }
    writeln!(out, "{}    }} ,", indent)?;
    // lens
    writeln!(out, "{}    lens {{", indent)?;
    for (i, &l) in lens.iter().enumerate() {
        let comma = if i + 1 < lens.len() { " ," } else { "" };
        writeln!(out, "{}      {}{}", indent, l, comma)?;
    }
    writeln!(out, "{}    }}", indent)?;
    writeln!(out, "{}  }}", indent)?;
    writeln!(out, "{}}}", indent)
}

/// Write an ASN.1 REAL value in NCBI's mantissa/base/exponent notation.
fn write_asn_real(out: &mut impl Write, indent: &str, val: f64) -> io::Result<()> {
    if val == 0.0 {
        writeln!(out, "{}value real {{ 0, 10, 0 }}", indent)
    } else {
        // Represent as mantissa * 10^exponent where mantissa is an integer
        let s = format!("{:.15e}", val);
        let parts: Vec<&str> = s.split('e').collect();
        let mantissa_str = parts[0]; // e.g. "1.234000000000000"
        let exp: i32 = parts[1].parse().unwrap_or(0);

        // Remove decimal point and trailing zeros to get integer mantissa
        let cleaned = mantissa_str.replace('.', "");
        let trimmed = cleaned.trim_end_matches('0');
        let trimmed = if trimmed.is_empty() { "0" } else { trimmed };

        // Number of digits after decimal point that we kept
        let frac_digits = mantissa_str.len() - mantissa_str.find('.').unwrap_or(mantissa_str.len()) - 1;
        let actual_exp = exp - frac_digits as i32 + (cleaned.len() - trimmed.len()) as i32;

        // mantissa as integer, exponent adjusted
        let mantissa_int: i64 = trimmed.parse().unwrap_or(0);
        let sign = if val < 0.0 && mantissa_int > 0 { -1i64 } else { 1 };
        writeln!(out, "{}value real {{ {}, 10, {} }}", indent, mantissa_int * sign, actual_exp)
    }
}

fn fmt8_asn_text(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
) -> io::Result<()> {
    writeln!(out, "Seq-annot ::= {{")?;
    writeln!(out, "  desc {{")?;
    writeln!(out, "    name \"BLAST {}\" ,", asn_escape(ctx.program))?;
    writeln!(out, "    title \"Database: {}\"", asn_escape(ctx.db_title))?;
    writeln!(out, "  }} ,")?;
    writeln!(out, "  data align {{")?;

    let total_hsps: usize = results.iter().map(|r| r.hsps.len()).sum();
    let mut hsp_idx = 0;
    for r in results {
        for hsp in &r.hsps {
            write_seqalign_text(out, "    ", ctx, r, hsp)?;
            hsp_idx += 1;
            if hsp_idx < total_hsps {
                writeln!(out, "    ,")?;
            }
        }
    }

    writeln!(out, "  }}")?;
    writeln!(out, "}}")
}

// ─── Format 9: Binary ASN.1 (Seq-annot, BER) ─────────────────────────────

/// Encode an ASN.1 tag in BER.
fn ber_tag(class: u8, constructed: bool, tag_num: u32) -> Vec<u8> {
    let class_bits = class << 6;
    let constructed_bit = if constructed { 0x20 } else { 0x00 };

    if tag_num < 31 {
        vec![class_bits | constructed_bit | (tag_num as u8)]
    } else {
        let mut bytes = vec![class_bits | constructed_bit | 0x1F];
        let mut t = tag_num;
        let mut stack = Vec::new();
        while t > 0 {
            stack.push((t & 0x7F) as u8);
            t >>= 7;
        }
        for (i, &b) in stack.iter().rev().enumerate() {
            if i + 1 < stack.len() {
                bytes.push(b | 0x80);
            } else {
                bytes.push(b);
            }
        }
        bytes
    }
}

/// Encode BER length.
fn ber_length(len: usize) -> Vec<u8> {
    if len < 128 {
        vec![len as u8]
    } else if len < 256 {
        vec![0x81, len as u8]
    } else if len < 65536 {
        vec![0x82, (len >> 8) as u8, (len & 0xFF) as u8]
    } else if len < 16777216 {
        vec![0x83, (len >> 16) as u8, ((len >> 8) & 0xFF) as u8, (len & 0xFF) as u8]
    } else {
        vec![0x84, (len >> 24) as u8, ((len >> 16) & 0xFF) as u8, ((len >> 8) & 0xFF) as u8, (len & 0xFF) as u8]
    }
}

/// Encode a BER TLV.
fn ber_tlv(class: u8, constructed: bool, tag_num: u32, value: &[u8]) -> Vec<u8> {
    let mut out = ber_tag(class, constructed, tag_num);
    out.extend_from_slice(&ber_length(value.len()));
    out.extend_from_slice(value);
    out
}

/// BER INTEGER encoding.
fn ber_integer(val: i64) -> Vec<u8> {
    let mut bytes = Vec::new();
    if val == 0 {
        bytes.push(0);
    } else {
        let mut v = val;
        let mut stack = Vec::new();
        if val > 0 {
            while v > 0 {
                stack.push((v & 0xFF) as u8);
                v >>= 8;
            }
            // ensure high bit clear for positive
            if stack.last().is_some_and(|&b| b & 0x80 != 0) {
                stack.push(0);
            }
        } else {
            while v < -1 || (stack.last().is_none_or(|&b| b & 0x80 == 0)) {
                stack.push((v & 0xFF) as u8);
                v >>= 8;
                if stack.len() > 8 { break; }
            }
        }
        stack.reverse();
        bytes = stack;
    }
    ber_tlv(0, false, 2, &bytes) // UNIVERSAL INTEGER tag=2
}

/// BER REAL encoding (base-10 NR3 form).
fn ber_real(val: f64) -> Vec<u8> {
    if val == 0.0 {
        return ber_tlv(0, false, 9, &[]); // zero-length for 0.0
    }
    // Use NR3 form: mantissa.Eexponent as ASCII
    let s = format!("{:.15E}", val);
    let nr3 = format!(" {}", s); // leading space = NR3 form indicator
    let mut content = vec![0x03u8]; // NR3 form
    content.extend_from_slice(nr3.as_bytes());
    ber_tlv(0, false, 9, &content) // UNIVERSAL REAL tag=9
}

/// BER VisibleString.
fn ber_visiblestring(s: &str) -> Vec<u8> {
    ber_tlv(0, false, 26, s.as_bytes()) // UNIVERSAL VisibleString tag=26
}

/// BER ENUMERATED.
fn ber_enumerated(val: i32) -> Vec<u8> {
    let mut bytes = Vec::new();
    if (0..128).contains(&val) {
        bytes.push(val as u8);
    } else {
        let mut v = val as i64;
        let mut stack = Vec::new();
        if val >= 0 {
            while v > 0 { stack.push((v & 0xFF) as u8); v >>= 8; }
            if stack.last().is_some_and(|&b| b & 0x80 != 0) { stack.push(0); }
        } else {
            while v < -1 || stack.last().is_none_or(|&b| b & 0x80 == 0) {
                stack.push((v & 0xFF) as u8); v >>= 8;
                if stack.len() > 4 { break; }
            }
        }
        stack.reverse();
        bytes = stack;
    }
    ber_tlv(0, false, 10, &bytes) // UNIVERSAL ENUMERATED tag=10
}

/// Wrap content with a context-specific constructed tag.
fn ber_context(tag_num: u32, content: &[u8]) -> Vec<u8> {
    ber_tlv(2, true, tag_num, content) // CONTEXT-SPECIFIC, CONSTRUCTED
}

/// BER SEQUENCE (UNIVERSAL 16 constructed).
fn ber_sequence(content: &[u8]) -> Vec<u8> {
    ber_tlv(0, true, 16, content)
}

/// BER SET OF (UNIVERSAL 17 constructed).
fn ber_set_of(content: &[u8]) -> Vec<u8> {
    ber_tlv(0, true, 17, content)
}

/// Build a BER-encoded Score structure.
fn ber_score_int(name: &str, val: i64) -> Vec<u8> {
    // Score ::= SEQUENCE { id [0] Object-id, value [1] CHOICE { int INTEGER } }
    let id = ber_context(0, &ber_visiblestring(name)); // Object-id as VisibleString via [0]
    let value_inner = ber_integer(val);
    let value = ber_context(1, &value_inner);
    let mut content = Vec::new();
    content.extend_from_slice(&id);
    content.extend_from_slice(&value);
    ber_sequence(&content)
}

fn ber_score_real(name: &str, val: f64) -> Vec<u8> {
    let id = ber_context(0, &ber_visiblestring(name));
    let value_inner = ber_real(val);
    let value = ber_context(1, &value_inner);
    let mut content = Vec::new();
    content.extend_from_slice(&id);
    content.extend_from_slice(&value);
    ber_sequence(&content)
}

/// Build a BER Seq-align for one HSP.
fn ber_seqalign(ctx: &SearchContext<'_>, r: &SearchResult, hsp: &Hsp) -> Vec<u8> {
    let (starts_data, lens_data, numseg) = build_denseg(hsp);
    let qid = first_word(ctx.query_title);
    let sid = if !r.subject_accession.is_empty() { &r.subject_accession } else { first_word(&r.subject_title) };

    // type [0] ENUMERATED: partial=3
    let sa_type = ber_context(0, &ber_enumerated(3));
    // dim [1] INTEGER
    let sa_dim = ber_context(1, &ber_integer(2));

    // score [2] SET OF Score
    let mut scores_content = Vec::new();
    scores_content.extend_from_slice(&ber_score_int("score", hsp.score as i64));
    scores_content.extend_from_slice(&ber_score_real("e_value", hsp.evalue));
    scores_content.extend_from_slice(&ber_score_real("bit_score", hsp.bit_score));
    scores_content.extend_from_slice(&ber_score_int("num_ident", hsp.num_identities as i64));
    let sa_scores = ber_context(2, &ber_set_of(&scores_content));

    // segs [3] CHOICE denseg [1]
    // Dense-seg ::= SEQUENCE { dim, numseg, ids, starts, lens }
    let ds_dim = ber_context(0, &ber_integer(2));
    let ds_numseg = ber_context(1, &ber_integer(numseg as i64));

    // ids [2] SEQUENCE OF Seq-id; using local str
    let qid_ber = ber_context(0, &ber_visiblestring(qid)); // Seq-id local [0]
    let sid_ber = ber_context(0, &ber_visiblestring(sid));
    let mut ids_content = Vec::new();
    // Wrap each as a CHOICE: local [0] Object-id str VisibleString
    ids_content.extend_from_slice(&ber_sequence(&qid_ber));
    ids_content.extend_from_slice(&ber_sequence(&sid_ber));
    let ds_ids = ber_context(2, &ber_sequence(&ids_content));

    // starts [3] SEQUENCE OF INTEGER
    let mut starts_enc = Vec::new();
    for &s in &starts_data {
        starts_enc.extend_from_slice(&ber_integer(s));
    }
    let ds_starts = ber_context(3, &ber_sequence(&starts_enc));

    // lens [4] SEQUENCE OF INTEGER
    let mut lens_enc = Vec::new();
    for &l in &lens_data {
        lens_enc.extend_from_slice(&ber_integer(l as i64));
    }
    let ds_lens = ber_context(4, &ber_sequence(&lens_enc));

    let mut denseg_content = Vec::new();
    denseg_content.extend_from_slice(&ds_dim);
    denseg_content.extend_from_slice(&ds_numseg);
    denseg_content.extend_from_slice(&ds_ids);
    denseg_content.extend_from_slice(&ds_starts);
    denseg_content.extend_from_slice(&ds_lens);

    let sa_segs = ber_context(3, &ber_context(1, &ber_sequence(&denseg_content))); // denseg is CHOICE [1]

    let mut sa_content = Vec::new();
    sa_content.extend_from_slice(&sa_type);
    sa_content.extend_from_slice(&sa_dim);
    sa_content.extend_from_slice(&sa_scores);
    sa_content.extend_from_slice(&sa_segs);
    ber_sequence(&sa_content)
}

fn fmt9_asn_binary(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
) -> io::Result<()> {
    // Seq-annot ::= SEQUENCE { desc [0], data [1] CHOICE align [2] SET OF Seq-align }
    let db_title = ctx.db_title;
    let prog_name = format!("BLAST {}", ctx.program);

    // desc [0]: Annot-descr ::= SET OF Annotdesc
    // name [0] VisibleString, title [1] VisibleString
    let name_desc = ber_context(0, &ber_visiblestring(&prog_name));
    let title_desc = ber_context(1, &ber_visiblestring(&format!("Database: {}", db_title)));
    let mut desc_content = Vec::new();
    desc_content.extend_from_slice(&name_desc);
    desc_content.extend_from_slice(&title_desc);
    let sa_desc = ber_context(0, &ber_set_of(&desc_content));

    // data [1] CHOICE align [2] SET OF Seq-align
    let mut aligns = Vec::new();
    for r in results {
        for hsp in &r.hsps {
            aligns.extend_from_slice(&ber_seqalign(ctx, r, hsp));
        }
    }
    let sa_data = ber_context(1, &ber_context(2, &ber_set_of(&aligns)));

    let mut content = Vec::new();
    content.extend_from_slice(&sa_desc);
    content.extend_from_slice(&sa_data);
    let encoded = ber_sequence(&content);

    out.write_all(&encoded)
}

// ─── Format 11: BLAST Archive (Blast4-archive, text ASN.1) ───────────────

fn fmt11_archive(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
) -> io::Result<()> {
    let qid = first_word(ctx.query_title);

    writeln!(out, "Blast4-archive ::= {{")?;

    // request
    writeln!(out, "  request {{")?;
    writeln!(out, "    ident \"blast-rs\" ,")?;
    writeln!(out, "    body queue-search {{")?;
    writeln!(out, "      program \"{}\" ,", asn_escape(ctx.program))?;
    writeln!(out, "      service \"plain\" ,")?;
    writeln!(out, "      queries bioseq-set {{")?;
    writeln!(out, "        seq-set {{")?;
    writeln!(out, "          seq {{")?;
    writeln!(out, "            id {{")?;
    writeln!(out, "              local str \"{}\"", asn_escape(qid))?;
    writeln!(out, "            }} ,")?;
    writeln!(out, "            descr {{")?;
    writeln!(out, "              title \"{}\"", asn_escape(ctx.query_title))?;
    writeln!(out, "            }} ,")?;
    writeln!(out, "            inst {{")?;
    writeln!(out, "              repr raw ,")?;
    // Determine molecule type from program
    let mol = if matches!(ctx.program, "blastn" | "tblastx" | "tblastn") { "na" } else { "aa" };
    writeln!(out, "              mol {} ,", mol)?;
    writeln!(out, "              length {}", ctx.query_len)?;
    writeln!(out, "            }}")?;
    writeln!(out, "          }}")?;
    writeln!(out, "        }}")?;
    writeln!(out, "      }} ,")?;
    writeln!(out, "      subject database \"{}\" ,", asn_escape(ctx.db_path))?;
    writeln!(out, "      algorithm-options {{")?;
    writeln!(out, "        {{")?;
    writeln!(out, "          name \"EvalueThreshold\" ,")?;
    write!(out, "          value cutoff e-value ")?;
    write_asn_real_inline(out, ctx.evalue_threshold)?;
    writeln!(out)?;
    writeln!(out, "        }} ,")?;
    writeln!(out, "        {{")?;
    writeln!(out, "          name \"GapOpeningCost\" ,")?;
    writeln!(out, "          value integer {}", ctx.gap_open)?;
    writeln!(out, "        }} ,")?;
    writeln!(out, "        {{")?;
    writeln!(out, "          name \"GapExtensionCost\" ,")?;
    writeln!(out, "          value integer {}", ctx.gap_extend)?;
    writeln!(out, "        }} ,")?;
    writeln!(out, "        {{")?;
    writeln!(out, "          name \"MatrixName\" ,")?;
    writeln!(out, "          value string \"{}\"", asn_escape(ctx.matrix))?;
    writeln!(out, "        }}")?;
    writeln!(out, "      }}")?;
    writeln!(out, "    }}")?;
    writeln!(out, "  }} ,")?;

    // results
    writeln!(out, "  results {{")?;

    // alignments
    if !results.is_empty() {
        writeln!(out, "    alignments {{")?;
        let total_hsps: usize = results.iter().map(|r| r.hsps.len()).sum();
        let mut hsp_idx = 0;
        for r in results {
            for hsp in &r.hsps {
                write_seqalign_text(out, "      ", ctx, r, hsp)?;
                hsp_idx += 1;
                if hsp_idx < total_hsps {
                    writeln!(out, "      ,")?;
                }
            }
        }
        writeln!(out, "    }} ,")?;
    }

    // search-stats
    writeln!(out, "    search-stats {{")?;
    writeln!(out, "      \"db-num={}\" ,", ctx.db_num_seqs)?;
    writeln!(out, "      \"db-len={}\"", ctx.db_len)?;
    writeln!(out, "    }}")?;

    writeln!(out, "  }}")?;
    writeln!(out, "}}")
}

/// Write ASN.1 REAL value inline as {{ mantissa, 10, exponent }}.
fn write_asn_real_inline(out: &mut impl Write, val: f64) -> io::Result<()> {
    if val == 0.0 {
        write!(out, "{{ 0, 10, 0 }}")
    } else {
        // Decompose into integer mantissa * 10^exp
        let s = format!("{:.15e}", val);
        let parts: Vec<&str> = s.split('e').collect();
        let mantissa_str = parts[0];
        let exp: i32 = parts[1].parse().unwrap_or(0);
        let cleaned = mantissa_str.replace('.', "");
        let trimmed = cleaned.trim_end_matches('0');
        let trimmed = if trimmed.is_empty() { "0" } else { trimmed };
        let frac_digits = mantissa_str.len() - mantissa_str.find('.').unwrap_or(mantissa_str.len()) - 1;
        let actual_exp = exp - frac_digits as i32 + (cleaned.len() - trimmed.len()) as i32;
        let mantissa_int: i64 = trimmed.parse().unwrap_or(0);
        let sign = if val < 0.0 && mantissa_int > 0 { -1i64 } else { 1 };
        write!(out, "{{ {}, 10, {} }}", mantissa_int * sign, actual_exp)
    }
}

/// Escape a string for text ASN.1 (double-quote delimited).
fn asn_escape(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}

// ─── Format 12: Seqalign JSON ─────────────────────────────────────────────

fn fmt12_seqalign_json(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
) -> io::Result<()> {
    writeln!(out, "{{")?;
    writeln!(out, "  \"query\": {},", json_str(ctx.query_title))?;
    writeln!(out, "  \"query_len\": {},", ctx.query_len)?;
    writeln!(out, "  \"program\": {},", json_str(ctx.program))?;
    writeln!(out, "  \"database\": {},", json_str(ctx.db_path))?;
    writeln!(out, "  \"alignments\": [")?;

    for (hi, r) in results.iter().enumerate() {
        let sid = if !r.subject_accession.is_empty() { &r.subject_accession } else { &r.subject_title };
        let comma_r = if hi + 1 < results.len() { "," } else { "" };
        writeln!(out, "    {{")?;
        writeln!(out, "      \"subject_id\": {},", json_str(sid))?;
        writeln!(out, "      \"subject_title\": {},", json_str(&r.subject_title))?;
        writeln!(out, "      \"subject_len\": {},", r.subject_len)?;
        writeln!(out, "      \"hsps\": [")?;

        for (hj, hsp) in r.hsps.iter().enumerate() {
            let comma_h = if hj + 1 < r.hsps.len() { "," } else { "" };
            writeln!(out, "        {{")?;
            writeln!(out, "          \"score\": {},", hsp.score)?;
            writeln!(out, "          \"bit_score\": {:.1},", hsp.bit_score)?;
            writeln!(out, "          \"evalue\": {:.2e},", hsp.evalue)?;
            writeln!(out, "          \"identity\": {},", hsp.num_identities)?;
            writeln!(out, "          \"align_len\": {},", hsp.alignment_length)?;
            writeln!(out, "          \"gaps\": {},", hsp.num_gaps)?;
            writeln!(out, "          \"query_from\": {},", hsp.query_start + 1)?;
            writeln!(out, "          \"query_to\": {},", hsp.query_end)?;
            writeln!(out, "          \"subject_from\": {},", hsp.subject_start + 1)?;
            writeln!(out, "          \"subject_to\": {},", hsp.subject_end)?;
            writeln!(out, "          \"query_strand\": {},", json_str(if hsp.query_frame < 0 { "Minus" } else { "Plus" }))?;
            writeln!(out, "          \"subject_strand\": {},", json_str(if hsp.subject_frame < 0 { "Minus" } else { "Plus" }))?;
            writeln!(out, "          \"qseq\": {},", json_str(str_from_aln(&hsp.query_aln)))?;
            writeln!(out, "          \"sseq\": {},", json_str(str_from_aln(&hsp.subject_aln)))?;
            writeln!(out, "          \"midline\": {}", json_str(str_from_aln(&hsp.midline)))?;
            writeln!(out, "        }}{}", comma_h)?;
        }
        writeln!(out, "      ]")?;
        writeln!(out, "    }}{}", comma_r)?;
    }

    writeln!(out, "  ]")?;
    writeln!(out, "}}")
}

// ─── Format 13: Multi-file BLAST JSON ─────────────────────────────────────

fn fmt13_json_multi(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
) -> io::Result<()> {
    // Output same as format 15 (single-file JSON) since we don't do actual multi-file
    fmt15_json(out, ctx, results)
}

// ─── Format 14: Multi-file BLAST XML2 ─────────────────────────────────────

fn fmt14_xml2_multi(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
) -> io::Result<()> {
    // Output same as format 16 (single-file XML2) since we don't do actual multi-file
    fmt16_xml2(out, ctx, results)
}

// ─── Format 15: Single-file BLAST JSON ─────────────────────────────────────

fn fmt15_json(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
) -> io::Result<()> {
    writeln!(out, "  {{")?;
    writeln!(out, "    \"BlastOutput2\": [{{")?;
    writeln!(out, "      \"report\": {{")?;
    writeln!(out, "        \"program\": {},", json_str(ctx.program))?;
    writeln!(out, "        \"version\": \"blast-rs 0.1.0\",")?;
    writeln!(out, "        \"reference\": \"Rust BLAST implementation\",")?;
    writeln!(out, "        \"search_target\": {{ \"db\": {} }},", json_str(ctx.db_path))?;
    writeln!(out, "        \"params\": {{")?;
    writeln!(out, "          \"matrix\": {},", json_str(ctx.matrix))?;
    writeln!(out, "          \"expect\": {},", ctx.evalue_threshold)?;
    writeln!(out, "          \"gap_open\": {},", ctx.gap_open)?;
    writeln!(out, "          \"gap_extend\": {}", ctx.gap_extend)?;
    writeln!(out, "        }},")?;
    writeln!(out, "        \"results\": {{")?;
    writeln!(out, "          \"search\": {{")?;
    writeln!(out, "            \"query_id\": {},", json_str(&format!("Query_{}", ctx.iter_num)))?;
    writeln!(out, "            \"query_title\": {},", json_str(ctx.query_title))?;
    writeln!(out, "            \"query_len\": {},", ctx.query_len)?;
    writeln!(out, "            \"hits\": [")?;

    for (hi, r) in results.iter().enumerate() {
        let sid = if !r.subject_accession.is_empty() { &r.subject_accession } else { &r.subject_title };
        let comma_r = if hi + 1 < results.len() { "," } else { "" };
        writeln!(out, "              {{")?;
        writeln!(out, "                \"num\": {},", hi + 1)?;
        writeln!(out, "                \"description\": [{{")?;
        writeln!(out, "                  \"id\": {},", json_str(sid))?;
        writeln!(out, "                  \"title\": {}", json_str(&r.subject_title))?;
        writeln!(out, "                }}],")?;
        writeln!(out, "                \"len\": {},", r.subject_len)?;
        writeln!(out, "                \"hsps\": [")?;

        for (hj, hsp) in r.hsps.iter().enumerate() {
            let positives = count_positive_chars(&hsp.midline);
            let _mismatches = hsp.alignment_length.saturating_sub(hsp.num_identities + hsp.num_gaps);
            let comma_h = if hj + 1 < r.hsps.len() { "," } else { "" };
            writeln!(out, "                  {{")?;
            writeln!(out, "                    \"num\": {},", hj + 1)?;
            writeln!(out, "                    \"bit_score\": {:.4},", hsp.bit_score)?;
            writeln!(out, "                    \"score\": {},", hsp.score)?;
            writeln!(out, "                    \"evalue\": {:.6e},", hsp.evalue)?;
            writeln!(out, "                    \"identity\": {},", hsp.num_identities)?;
            writeln!(out, "                    \"positive\": {},", positives)?;
            writeln!(out, "                    \"gaps\": {},", hsp.num_gaps)?;
            writeln!(out, "                    \"align_len\": {},", hsp.alignment_length)?;
            writeln!(out, "                    \"query_from\": {},", hsp.query_start + 1)?;
            writeln!(out, "                    \"query_to\": {},", hsp.query_end)?;
            writeln!(out, "                    \"hit_from\": {},", hsp.subject_start + 1)?;
            writeln!(out, "                    \"hit_to\": {},", hsp.subject_end)?;
            writeln!(out, "                    \"qseq\": {},", json_str(str_from_aln(&hsp.query_aln)))?;
            writeln!(out, "                    \"hseq\": {},", json_str(str_from_aln(&hsp.subject_aln)))?;
            writeln!(out, "                    \"midline\": {}", json_str(str_from_aln(&hsp.midline)))?;
            writeln!(out, "                  }}{}", comma_h)?;
        }

        writeln!(out, "                ]")?;
        writeln!(out, "              }}{}", comma_r)?;
    }

    writeln!(out, "            ],")?;
    writeln!(out, "            \"stat\": {{")?;
    writeln!(out, "              \"db_num\": {},", ctx.db_num_seqs)?;
    writeln!(out, "              \"db_len\": {}", ctx.db_len)?;
    writeln!(out, "            }}")?;
    writeln!(out, "          }}")?;
    writeln!(out, "        }}")?;
    writeln!(out, "      }}")?;
    writeln!(out, "    }}]")?;
    writeln!(out, "  }}")
}

// ─── Format 16: Single-file BLAST XML2 ─────────────────────────────────────

fn fmt16_xml2(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
) -> io::Result<()> {
    writeln!(out, r#"<?xml version="1.0"?>"#)?;
    writeln!(out, r#"<BlastXML2>"#)?;
    writeln!(out, "  <BlastXML2_report>")?;
    writeln!(out, "    <Report>")?;
    writeln!(out, "      <Report_program>{}</Report_program>", xml_escape(ctx.program))?;
    writeln!(out, "      <Report_version>blast-rs 0.1.0</Report_version>")?;
    writeln!(out, "      <Report_search-target>")?;
    writeln!(out, "        <Target>")?;
    writeln!(out, "          <Target_db>{}</Target_db>", xml_escape(ctx.db_path))?;
    writeln!(out, "        </Target>")?;
    writeln!(out, "      </Report_search-target>")?;
    writeln!(out, "      <Report_params>")?;
    writeln!(out, "        <Parameters>")?;
    writeln!(out, "          <Parameters_matrix>{}</Parameters_matrix>", xml_escape(ctx.matrix))?;
    writeln!(out, "          <Parameters_expect>{}</Parameters_expect>", ctx.evalue_threshold)?;
    writeln!(out, "          <Parameters_gap-open>{}</Parameters_gap-open>", ctx.gap_open)?;
    writeln!(out, "          <Parameters_gap-extend>{}</Parameters_gap-extend>", ctx.gap_extend)?;
    writeln!(out, "        </Parameters>")?;
    writeln!(out, "      </Report_params>")?;
    writeln!(out, "      <Report_results>")?;
    writeln!(out, "        <Results>")?;
    writeln!(out, "          <Results_search>")?;
    writeln!(out, "            <Search>")?;
    writeln!(out, "              <Search_query-id>Query_{}</Search_query-id>", ctx.iter_num)?;
    writeln!(out, "              <Search_query-title>{}</Search_query-title>", xml_escape(ctx.query_title))?;
    writeln!(out, "              <Search_query-len>{}</Search_query-len>", ctx.query_len)?;
    writeln!(out, "              <Search_hits>")?;

    for (hi, r) in results.iter().enumerate() {
        let sid = if !r.subject_accession.is_empty() { &r.subject_accession } else { &r.subject_title };
        writeln!(out, "                <Hit>")?;
        writeln!(out, "                  <Hit_num>{}</Hit_num>", hi + 1)?;
        writeln!(out, "                  <Hit_description>")?;
        writeln!(out, "                    <HitDescr>")?;
        writeln!(out, "                      <HitDescr_id>{}</HitDescr_id>", xml_escape(sid))?;
        writeln!(out, "                      <HitDescr_title>{}</HitDescr_title>", xml_escape(&r.subject_title))?;
        writeln!(out, "                    </HitDescr>")?;
        writeln!(out, "                  </Hit_description>")?;
        writeln!(out, "                  <Hit_len>{}</Hit_len>", r.subject_len)?;
        writeln!(out, "                  <Hit_hsps>")?;

        for (hj, hsp) in r.hsps.iter().enumerate() {
            let positives = count_positive_chars(&hsp.midline);
            writeln!(out, "                    <Hsp>")?;
            writeln!(out, "                      <Hsp_num>{}</Hsp_num>", hj + 1)?;
            writeln!(out, "                      <Hsp_bit-score>{:.4}</Hsp_bit-score>", hsp.bit_score)?;
            writeln!(out, "                      <Hsp_score>{}</Hsp_score>", hsp.score)?;
            writeln!(out, "                      <Hsp_evalue>{:.6e}</Hsp_evalue>", hsp.evalue)?;
            writeln!(out, "                      <Hsp_query-from>{}</Hsp_query-from>", hsp.query_start + 1)?;
            writeln!(out, "                      <Hsp_query-to>{}</Hsp_query-to>", hsp.query_end)?;
            writeln!(out, "                      <Hsp_hit-from>{}</Hsp_hit-from>", hsp.subject_start + 1)?;
            writeln!(out, "                      <Hsp_hit-to>{}</Hsp_hit-to>", hsp.subject_end)?;
            writeln!(out, "                      <Hsp_identity>{}</Hsp_identity>", hsp.num_identities)?;
            writeln!(out, "                      <Hsp_positive>{}</Hsp_positive>", positives)?;
            writeln!(out, "                      <Hsp_gaps>{}</Hsp_gaps>", hsp.num_gaps)?;
            writeln!(out, "                      <Hsp_align-len>{}</Hsp_align-len>", hsp.alignment_length)?;
            writeln!(out, "                      <Hsp_qseq>{}</Hsp_qseq>", xml_escape(str_from_aln(&hsp.query_aln)))?;
            writeln!(out, "                      <Hsp_hseq>{}</Hsp_hseq>", xml_escape(str_from_aln(&hsp.subject_aln)))?;
            writeln!(out, "                      <Hsp_midline>{}</Hsp_midline>", xml_escape(str_from_aln(&hsp.midline)))?;
            writeln!(out, "                    </Hsp>")?;
        }

        writeln!(out, "                  </Hit_hsps>")?;
        writeln!(out, "                </Hit>")?;
    }

    writeln!(out, "              </Search_hits>")?;
    writeln!(out, "              <Search_stat>")?;
    writeln!(out, "                <Statistics>")?;
    writeln!(out, "                  <Statistics_db-num>{}</Statistics_db-num>", ctx.db_num_seqs)?;
    writeln!(out, "                  <Statistics_db-len>{}</Statistics_db-len>", ctx.db_len)?;
    writeln!(out, "                </Statistics>")?;
    writeln!(out, "              </Search_stat>")?;
    writeln!(out, "            </Search>")?;
    writeln!(out, "          </Results_search>")?;
    writeln!(out, "        </Results>")?;
    writeln!(out, "      </Report_results>")?;
    writeln!(out, "    </Report>")?;
    writeln!(out, "  </BlastXML2_report>")?;
    writeln!(out, "</BlastXML2>")
}

// ─── Format 17: Subject FASTA sequences ────────────────────────────────────

fn fmt17_fasta(
    out: &mut impl Write,
    _ctx: &SearchContext<'_>,
    results: &[SearchResult],
    db: Option<&BlastDb>,
) -> io::Result<()> {
    let db = match db {
        Some(d) => d,
        None => {
            writeln!(out, "# Error: database not available for format 17")?;
            return Ok(());
        }
    };

    for r in results {
        let sid = subject_title(r);
        writeln!(out, ">{}", sid)?;

        let seq: Vec<u8> = match db.seq_type() {
            blast_rs::db::index::SeqType::Protein => {
                db.get_sequence_protein(r.subject_oid).unwrap_or_default()
            }
            blast_rs::db::index::SeqType::Nucleotide => {
                db.get_sequence_nucleotide(r.subject_oid).unwrap_or_default()
            }
        };

        for chunk in seq.chunks(60) {
            out.write_all(chunk)?;
            writeln!(out)?;
        }
    }
    Ok(())
}

// ─── Format 18: SAM ───────────────────────────────────────────────────────

fn fmt18_sam(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
) -> io::Result<()> {
    // SAM header
    writeln!(out, "@HD\tVN:1.6\tSO:queryname")?;
    for r in results {
        writeln!(out, "@SQ\tSN:{}\tLN:{}", r.subject_accession, r.subject_len)?;
    }
    writeln!(out, "@PG\tID:blast-cli\tPN:{}\tVN:0.1.0", ctx.program)?;

    let qname = first_word(ctx.query_title);

    for r in results {
        let rname = if !r.subject_accession.is_empty() { &r.subject_accession } else { &r.subject_title };

        for hsp in &r.hsps {
            // FLAG: 0 for forward, 16 for reverse
            let flag = if hsp.subject_frame < 0 || hsp.query_frame < 0 { 16 } else { 0 };

            // POS: 1-based subject start
            let pos = hsp.subject_start + 1;

            // MAPQ: derived from e-value (approximate)
            let mapq = if hsp.evalue <= 0.0 { 255 }
                else { ((-10.0 * hsp.evalue.log10()).clamp(0.0, 255.0)) as u8 };

            // Build CIGAR from alignment
            let cigar = build_cigar(&hsp.query_aln, &hsp.subject_aln);

            // SEQ: query sequence from alignment (gaps removed)
            let seq: String = hsp.query_aln.iter()
                .filter(|&&b| b != b'-')
                .map(|&b| b as char)
                .collect();

            // QUAL: unavailable
            let qual = "*";

            // Optional fields
            let as_tag = format!("AS:i:{}", hsp.score);
            let bs_tag = format!("BS:f:{:.1}", hsp.bit_score);
            let ev_tag = format!("EV:f:{:.2e}", hsp.evalue);
            let nm_tag = format!("NM:i:{}", hsp.alignment_length - hsp.num_identities);

            writeln!(out, "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}\t{}\t{}\t{}\t{}",
                qname, flag, rname, pos, mapq, cigar, seq, qual,
                as_tag, bs_tag, ev_tag, nm_tag)?;
        }
    }
    Ok(())
}

/// Build a CIGAR string from query and subject alignment strings.
fn build_cigar(query_aln: &[u8], subject_aln: &[u8]) -> String {
    let mut cigar = String::new();
    let mut last_op = ' ';
    let mut run = 0usize;

    for (&q, &s) in query_aln.iter().zip(subject_aln.iter()) {
        let op = if q == b'-' {
            'D' // deletion from reference (gap in query)
        } else if s == b'-' {
            'I' // insertion to reference (gap in subject)
        } else {
            'M' // alignment match (or mismatch)
        };

        if op == last_op {
            run += 1;
        } else {
            if run > 0 {
                cigar.push_str(&format!("{}{}", run, last_op));
            }
            last_op = op;
            run = 1;
        }
    }
    if run > 0 {
        cigar.push_str(&format!("{}{}", run, last_op));
    }

    if cigar.is_empty() { "*".to_string() } else { cigar }
}

// ─── HTML output ───────────────────────────────────────────────────────────

/// Write results in HTML format (wraps pairwise output in HTML).
pub fn write_results_html(
    out: &mut impl Write,
    ctx: &SearchContext<'_>,
    results: &[SearchResult],
) -> io::Result<()> {
    writeln!(out, "<!DOCTYPE html>")?;
    writeln!(out, "<html><head><title>BLAST Results</title>")?;
    writeln!(out, "<style>")?;
    writeln!(out, "body {{ font-family: monospace; white-space: pre; }}")?;
    writeln!(out, ".header {{ color: #333; font-weight: bold; }}")?;
    writeln!(out, ".hit-title {{ color: #006; }}")?;
    writeln!(out, ".score {{ color: #060; }}")?;
    writeln!(out, "a {{ color: #009; }}")?;
    writeln!(out, "</style></head><body>")?;

    // Use pairwise format inside the HTML body
    fmt0_pairwise(out, ctx, results, false, false)?;

    writeln!(out, "</body></html>")
}

// ─── Shared alignment helpers ───────────────────────────────────────────────

/// Write the Score / Expect / Identity / Gaps header for one HSP.
fn write_hsp_header(out: &mut impl Write, hsp: &Hsp, _hsp_num: usize) -> io::Result<()> {
    writeln!(out, " Score = {:.1} bits ({}),  Expect = {:.2e}", hsp.bit_score, hsp.score, hsp.evalue)?;
    let positives = count_positive_chars(&hsp.midline);
    writeln!(out,
        " Identities = {}/{} ({:.0}%), Positives = {}/{} ({:.0}%), Gaps = {}/{}",
        hsp.num_identities, hsp.alignment_length, hsp.percent_identity(),
        positives, hsp.alignment_length, 100.0 * positives as f64 / hsp.alignment_length.max(1) as f64,
        hsp.num_gaps, hsp.alignment_length,
    )?;
    writeln!(out)
}

/// Print the alignment in 60-character blocks.
///
/// `query_anchored`: formats 1-4 — replace '+' midline chars with ' ' or use '+' for all.
/// `positive_marks`: formats 2/4 — use '+' for everything non-gap, ' ' for mismatch.
fn write_alignment_blocks(
    out: &mut impl Write,
    hsp: &Hsp,
    query_anchored: bool,
    positive_marks: bool,
    line_width: usize,
) -> io::Result<()> {
    let aln_len = hsp.query_aln.len();
    let mut pos = 0;
    let mut q_pos = hsp.query_start;
    let mut s_pos = hsp.subject_start;

    while pos < aln_len {
        let end = (pos + line_width).min(aln_len);
        let q_chunk = &hsp.query_aln[pos..end];
        let m_chunk = &hsp.midline[pos..end];
        let s_chunk = &hsp.subject_aln[pos..end];

        let q_chars = q_chunk.iter().filter(|&&b| b != b'-').count();
        let s_chars = s_chunk.iter().filter(|&&b| b != b'-').count();

        // Possibly remap midline for formats 1-4
        let midline_display: Vec<u8> = if query_anchored || positive_marks {
            m_chunk.iter().map(|&c| {
                if positive_marks {
                    // Format 2/4: '+' for identical, space for mismatch/gap
                    if c == b'|' || c == b'+' { b'+' }
                    else { b' ' }
                } else {
                    // Format 1/3: '|' for identical, space for everything else
                    if c == b'|' { b'|' } else { b' ' }
                }
            }).collect()
        } else {
            m_chunk.to_vec()
        };

        writeln!(out, "Query  {:>6}  {}  {}",
            q_pos + 1,
            std::str::from_utf8(q_chunk).unwrap_or("?"),
            q_pos + q_chars,
        )?;
        writeln!(out, "              {}",
            std::str::from_utf8(&midline_display).unwrap_or("?"),
        )?;
        writeln!(out, "Sbjct  {:>6}  {}  {}",
            s_pos + 1,
            std::str::from_utf8(s_chunk).unwrap_or("?"),
            s_pos + s_chars,
        )?;
        writeln!(out)?;

        q_pos += q_chars;
        s_pos += s_chars;
        pos = end;
    }
    Ok(())
}

/// Build BTOP (BLAST Traceback Operations) string.
/// Format: alternating run-lengths of matches and substitution pairs.
/// Matches: an integer (e.g., "5")
/// Substitutions: two characters (e.g., "TG" = query T aligned to subject G)
/// Gaps in query: "-" + subject char (e.g., "-A")
/// Gaps in subject: query char + "-" (e.g., "A-")
fn build_btop(query_aln: &[u8], subject_aln: &[u8]) -> String {
    let mut btop = String::new();
    let mut match_run = 0usize;

    let flush_matches = |s: &mut String, n: &mut usize| {
        if *n > 0 {
            s.push_str(&n.to_string());
            *n = 0;
        }
    };

    for (&q, &s) in query_aln.iter().zip(subject_aln.iter()) {
        if q == b'-' {
            flush_matches(&mut btop, &mut match_run);
            btop.push('-');
            btop.push(s as char);
        } else if s == b'-' {
            flush_matches(&mut btop, &mut match_run);
            btop.push(q as char);
            btop.push('-');
        } else if q == s {
            match_run += 1;
        } else {
            flush_matches(&mut btop, &mut match_run);
            btop.push(q as char);
            btop.push(s as char);
        }
    }
    flush_matches(&mut btop, &mut match_run);
    btop
}

// ─── Utility functions ──────────────────────────────────────────────────────

fn first_word(s: &str) -> &str {
    s.split_whitespace().next().unwrap_or(s)
}

fn subject_title(r: &SearchResult) -> String {
    if !r.subject_accession.is_empty() {
        format!("{} {}", r.subject_accession, r.subject_title)
    } else {
        r.subject_title.clone()
    }
}

fn str_from_aln(aln: &[u8]) -> &str {
    std::str::from_utf8(aln).unwrap_or("?")
}

fn count_gap_opens(aln: &[u8]) -> usize {
    let mut opens = 0;
    let mut in_gap = false;
    for &b in aln {
        if b == b'-' {
            if !in_gap { opens += 1; in_gap = true; }
        } else {
            in_gap = false;
        }
    }
    opens
}

/// Count '|' and '+' characters in the midline (identities + positives).
fn count_positive_chars(midline: &[u8]) -> usize {
    midline.iter().filter(|&&c| c == b'|' || c == b'+').count()
}

fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;")
     .replace('<', "&lt;")
     .replace('>', "&gt;")
     .replace('"', "&quot;")
     .replace('\'', "&apos;")
}

fn json_str(s: &str) -> String {
    let escaped = s
        .replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
        .replace('\r', "\\r")
        .replace('\t', "\\t");
    format!("\"{}\"", escaped)
}
