from __future__ import annotations
import argparse
import html as html_mod
import logging
import os
import re
import sys
import multiprocessing as mp
from collections import Counter
from typing import Dict, Iterator, List, Tuple

from telomerediff import __version__
from telomerediff.io import read_fasta_fastq, buffer_records
from telomerediff.repeats import parse_k
from telomerediff.engine import (
    RepeatBlock,
    run_telomerediff_pass1,
    run_telomerediff_pass2,
    summarize_blocks_incremental,
    block_to_bed_name,
)

logger = logging.getLogger("telomerediff")


def _setup_logging(verbose: bool = False) -> None:
    """Configure logging for telomerediff."""
    level = logging.DEBUG if verbose else logging.INFO
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    root = logging.getLogger("telomerediff")
    root.setLevel(level)
    root.addHandler(handler)


def _wrap_progress(records: Iterator[Tuple[str, str]], desc: str, disable: bool) -> Iterator[Tuple[str, str]]:
    """Wrap an iterator with a tqdm progress bar if available."""
    if disable:
        yield from records
        return
    try:
        from tqdm import tqdm
        for item in tqdm(records, desc=desc, unit=" reads", dynamic_ncols=True):
            yield item
    except ImportError:
        yield from records


def _sanitize_field(s: str) -> str:
    """Sanitize a string for use in tab-delimited fields (GFF3/VCF/BED)."""
    return re.sub(r'[\t\n\r;=]', '_', s)


def _write_gff3(blocks: List[RepeatBlock], path: str) -> None:
    """Write blocks in GFF3 format."""
    with open(path, "w") as f:
        f.write("##gff-version 3\n")
        for i, b in enumerate(blocks):
            read_id = _sanitize_field(b.read)
            attrs = (
                f"ID=telrep_{i};canonical={b.canonical};mode={b.mode};"
                f"copies={b.copies};confidence={b.confidence:.4f};entropy={b.entropy:.4f}"
            )
            f.write(
                f"{read_id}\ttelomerediff\trepeat_region\t{b.start + 1}\t{b.end}\t"
                f"{b.confidence:.4f}\t{b.strand}\t.\t{attrs}\n"
            )


def _write_vcf(blocks: List[RepeatBlock], path: str) -> None:
    """Write blocks in VCF-like format for telomeric repeat annotations."""
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.3\n")
        f.write('##INFO=<ID=MOTIF,Number=1,Type=String,Description="Canonical repeat motif">\n')
        f.write('##INFO=<ID=COPIES,Number=1,Type=Integer,Description="Number of repeat copies">\n')
        f.write('##INFO=<ID=RUNBP,Number=1,Type=Integer,Description="Total run length in bp">\n')
        f.write('##INFO=<ID=MODE,Number=1,Type=String,Description="Detection mode">\n')
        f.write('##INFO=<ID=CONF,Number=1,Type=Float,Description="Confidence score">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for i, b in enumerate(blocks):
            read_id = _sanitize_field(b.read)
            info = f"MOTIF={b.canonical};COPIES={b.copies};RUNBP={b.run_bp};MODE={b.mode};CONF={b.confidence:.4f}"
            f.write(
                f"{read_id}\t{b.start + 1}\tTELREP_{i}\t.\t<TELREPEAT>\t.\tPASS\t{info}\n"
            )


def _validate_output_path(path: str, label: str) -> None:
    """Check that the output directory exists and is writable before processing."""
    parent = os.path.dirname(os.path.abspath(path))
    if not os.path.isdir(parent):
        raise SystemExit(f"Error: output directory does not exist for {label}: {parent}")
    if not os.access(parent, os.W_OK):
        raise SystemExit(f"Error: output directory is not writable for {label}: {parent}")


def _validate_args(args: argparse.Namespace) -> None:
    """Validate CLI arguments early, before any processing."""
    if args.min_run_bp < 1:
        raise SystemExit("Error: --min-run-bp must be >= 1")
    if args.top_motifs_per_k < 1:
        raise SystemExit("Error: --top-motifs-per-k must be >= 1")
    if args.max_mismatch < 0:
        raise SystemExit("Error: --max-mismatch must be >= 0")
    if args.threads < 1:
        raise SystemExit("Error: --threads must be >= 1")
    if args.chunk_size < 1:
        raise SystemExit("Error: --chunk-size must be >= 1")
    if args.min_quality is not None and args.min_quality < 0:
        raise SystemExit("Error: --min-quality must be >= 0")

    # Validate input exists
    if args.input != "-" and not os.path.isfile(args.input):
        raise SystemExit(f"Error: input file not found: {args.input}")

    # Validate output paths
    _validate_output_path(args.out_per_read, "--out-per-read")
    _validate_output_path(args.out_bed, "--out-bed")
    _validate_output_path(args.out_summary, "--out-summary")
    if args.out_gff3:
        _validate_output_path(args.out_gff3, "--out-gff3")
    if args.out_vcf:
        _validate_output_path(args.out_vcf, "--out-vcf")
    if args.out_html:
        _validate_output_path(args.out_html, "--out-html")


def main() -> None:
    ap = argparse.ArgumentParser(
        prog="telomerediff",
        description="Alignment-free de novo discovery of telomeric repeats from long-read sequencing (ONT / PacBio).",
    )
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    # Input
    ap.add_argument("-i", "--input", required=True,
                    help="FASTA/FASTQ (.gz), BAM, CRAM, or '-' for stdin")
    ap.add_argument("--reference", default=None,
                    help="Reference FASTA for CRAM decoding")
    ap.add_argument("--min-quality", type=float, default=None,
                    help="Minimum average base quality to keep a read (FASTQ/BAM only)")
    ap.add_argument("--format", choices=["fasta", "fastq", "bam", "cram"], default=None,
                    help="Force input format (useful for stdin, default: auto-detect)")

    # Algorithm
    ap.add_argument("-k", "--k", default="4-15",
                    help="k-mer sizes: e.g. '6', '4-15', or '5,6,7'")
    ap.add_argument("--min-run-bp", type=int, default=150,
                    help="Minimum repeat block length (bp)")
    ap.add_argument("--top-motifs-per-k", type=int, default=10,
                    help="Top motifs per k for refine pass")
    ap.add_argument("--refine", action="store_true",
                    help="Enable error-tolerant refine pass")
    ap.add_argument("--max-mismatch", type=int, default=1,
                    help="Max mismatches per k-mer copy")

    # Performance
    ap.add_argument("-t", "--threads", type=int, default=max(1, mp.cpu_count() // 2))
    ap.add_argument("--chunk-size", type=int, default=2000)
    ap.add_argument("--no-progress", action="store_true",
                    help="Disable progress bar")
    ap.add_argument("--verbose", action="store_true",
                    help="Enable debug logging")

    # Output paths
    ap.add_argument("--out-per-read", default="telomerediff.per_read.tsv")
    ap.add_argument("--out-bed", default="telomerediff.blocks.bed")
    ap.add_argument("--out-summary", default="telomerediff.summary.tsv")
    ap.add_argument("--out-gff3", default=None,
                    help="Output GFF3 file (optional)")
    ap.add_argument("--out-vcf", default=None,
                    help="Output VCF file (optional)")
    ap.add_argument("--out-html", default=None,
                    help="Output self-contained HTML report (optional)")
    ap.add_argument("--out-plots", default=None,
                    help="Directory for plot PNG files (optional)")

    args = ap.parse_args()

    _setup_logging(args.verbose)
    _validate_args(args)

    k_values: List[int] = parse_k(args.k)
    threads = max(1, args.threads)
    chunk_size = max(200, args.chunk_size)
    io_kwargs = dict(
        reference=args.reference,
        min_quality=args.min_quality,
        format_hint=args.format,
    )

    # Stdin requires buffering since we need two passes
    is_stdin = args.input == "-"
    if is_stdin:
        logger.info("Reading from stdin — buffering all records into memory for two-pass analysis...")
        buffered = buffer_records(args.input, **io_kwargs)
        logger.info("Buffered %d records from stdin.", len(buffered))
        if not buffered:
            logger.warning("No records read from stdin. Exiting.")
            raise SystemExit(1)

    # Pass 1: discover top motifs (perfect)
    logger.info("Pass 1: Discovering motifs (k=%s)...", args.k)
    if is_stdin:
        records_p1 = iter(buffered)
    else:
        records_p1 = read_fasta_fastq(args.input, **io_kwargs)
    records_p1 = _wrap_progress(records_p1, "Pass 1", args.no_progress)

    motif_counts_by_k = run_telomerediff_pass1(
        records=records_p1,
        k_values=k_values,
        min_run_bp=args.min_run_bp,
        threads=threads,
        chunk_size=chunk_size,
    )

    # Check if pass 1 found anything
    total_motifs = sum(len(c) for c in motif_counts_by_k.values())
    if total_motifs == 0:
        logger.warning("Pass 1 found no telomeric motifs. Input may not contain telomeric repeats "
                        "at the current settings (k=%s, min-run-bp=%d).", args.k, args.min_run_bp)

    fuzzy_seeds: Dict[int, List[str]] = {}
    if args.refine:
        for k in k_values:
            tops = [m for m, _ in motif_counts_by_k[k].most_common(args.top_motifs_per_k)]
            if tops:
                fuzzy_seeds[k] = tops
        logger.debug("Fuzzy seeds: %s", fuzzy_seeds)

    # Pass 2: emit blocks + write outputs incrementally
    logger.info("Pass 2: Detecting blocks...")
    total_bp_by_motif: Counter = Counter()
    runs_by_motif: Counter = Counter()
    by_k: Dict[int, Counter] = {}
    by_strand: Dict[str, Counter] = {}
    all_blocks: List[RepeatBlock] = []

    if is_stdin:
        records_p2 = iter(buffered)
    else:
        records_p2 = read_fasta_fastq(args.input, **io_kwargs)
    records_p2 = _wrap_progress(records_p2, "Pass 2", args.no_progress)

    with open(args.out_per_read, "w") as tsv, open(args.out_bed, "w") as bed:
        tsv.write("read\tstart\tend\tk\tcanonical\tmode\tstrand\tcopies\trun_bp\tconfidence\tentropy\n")

        for block_batch in run_telomerediff_pass2(
            records=records_p2,
            k_values=k_values,
            min_run_bp=args.min_run_bp,
            fuzzy_seeds=fuzzy_seeds if args.refine else None,
            max_mismatch=args.max_mismatch,
            threads=threads,
            chunk_size=chunk_size,
        ):
            for b in block_batch:
                read_field = _sanitize_field(b.read)
                tsv.write(
                    f"{read_field}\t{b.start}\t{b.end}\t{b.k}\t{b.canonical}\t"
                    f"{b.mode}\t{b.strand}\t{b.copies}\t{b.run_bp}\t"
                    f"{b.confidence:.4f}\t{b.entropy:.4f}\n"
                )
                name = block_to_bed_name(b.canonical, b.k)
                score = min(1000, b.run_bp)
                bed.write(f"{read_field}\t{b.start}\t{b.end}\t{name}\t{score}\t{b.strand}\n")

            all_blocks.extend(block_batch)
            summarize_blocks_incremental(block_batch, total_bp_by_motif, runs_by_motif, by_k, by_strand)

    if not all_blocks:
        logger.warning("No telomeric blocks detected. Try adjusting --min-run-bp or -k values.")

    # Summary TSV (strand-aware)
    with open(args.out_summary, "w") as out:
        out.write("canonical\tk\tobserved_runs\ttotal_bp\tstrand_plus_runs\tstrand_minus_runs\n")
        for canon, bp in total_bp_by_motif.most_common():
            for k, cts in sorted(by_k.items()):
                if canon in cts:
                    sp = by_strand.get(f"{canon}:+", Counter()).get("runs", 0)
                    sm = by_strand.get(f"{canon}:-", Counter()).get("runs", 0)
                    out.write(f"{canon}\t{k}\t{cts[canon]}\t{bp}\t{sp}\t{sm}\n")

    outputs = [args.out_per_read, args.out_bed, args.out_summary]

    # Optional GFF3
    if args.out_gff3:
        _write_gff3(all_blocks, args.out_gff3)
        outputs.append(args.out_gff3)

    # Optional VCF
    if args.out_vcf:
        _write_vcf(all_blocks, args.out_vcf)
        outputs.append(args.out_vcf)

    # Optional plots and HTML report
    plot_images: Dict[str, str] = {}
    if args.out_plots or args.out_html:
        try:
            from telomerediff.visualize import (
                plot_length_distribution,
                plot_motif_abundance,
                plot_per_read_heatmap,
                plot_strand_distribution,
            )

            if args.out_plots:
                os.makedirs(args.out_plots, exist_ok=True)

            plot_funcs = [
                ("Block Length Distribution", plot_length_distribution, {"blocks": all_blocks}),
                ("Motif Abundance", plot_motif_abundance, {"total_bp_by_motif": total_bp_by_motif}),
                ("Block Locations", plot_per_read_heatmap, {"blocks": all_blocks}),
                ("Strand Distribution", plot_strand_distribution, {"blocks": all_blocks}),
            ]

            for title, func, kwargs in plot_funcs:
                if args.out_plots:
                    fname = title.lower().replace(" ", "_") + ".png"
                    fpath = os.path.join(args.out_plots, fname)
                    func(output_path=fpath, **kwargs)
                    outputs.append(fpath)

                if args.out_html:
                    b64 = func(**kwargs)
                    if b64:
                        plot_images[title] = b64

        except ImportError:
            logger.warning("matplotlib not installed, skipping plots.")

    # Optional HTML report
    if args.out_html:
        from telomerediff.report import generate_html_report
        report_html = generate_html_report(
            blocks=all_blocks,
            total_bp_by_motif=total_bp_by_motif,
            runs_by_motif=runs_by_motif,
            by_k=by_k,
            by_strand=by_strand,
            input_file=args.input,
            plot_images=plot_images,
        )
        with open(args.out_html, "w") as f:
            f.write(report_html)
        outputs.append(args.out_html)

    logger.info("Outputs written:\n%s", "\n".join(f"  - {o}" for o in outputs))


if __name__ == "__main__":
    main()
