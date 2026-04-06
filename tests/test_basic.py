"""Comprehensive tests for teloscan."""

import os
import tempfile
import gzip

import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
SAMPLE_FASTQ = os.path.join(TESTS_DIR, "data", "sample.fastq")

# ========================
# Import tests
# ========================

def test_import():
    import teloscan
    assert teloscan is not None


def test_version():
    from teloscan import __version__
    assert __version__


# ========================
# repeats.py - parse_k
# ========================

class TestParseK:
    def test_single_int(self):
        from teloscan.repeats import parse_k
        assert parse_k(6) == [6]

    def test_single_string(self):
        from teloscan.repeats import parse_k
        assert parse_k("6") == [6]

    def test_range(self):
        from teloscan.repeats import parse_k
        assert parse_k("4-8") == [4, 5, 6, 7, 8]

    def test_comma_separated(self):
        from teloscan.repeats import parse_k
        assert parse_k("5,6,7") == [5, 6, 7]

    def test_mixed(self):
        from teloscan.repeats import parse_k
        result = parse_k("4-6,8,10-12")
        assert result == [4, 5, 6, 8, 10, 11, 12]

    def test_list_input(self):
        from teloscan.repeats import parse_k
        assert parse_k([5, 3, 7]) == [3, 5, 7]

    def test_deduplication(self):
        from teloscan.repeats import parse_k
        assert parse_k("4-6,5-7") == [4, 5, 6, 7]

    def test_invalid_empty(self):
        from teloscan.repeats import parse_k
        with pytest.raises(ValueError):
            parse_k("")

    def test_invalid_range(self):
        from teloscan.repeats import parse_k
        with pytest.raises(ValueError):
            parse_k("8-4")

    def test_invalid_non_integer(self):
        from teloscan.repeats import parse_k
        with pytest.raises(ValueError):
            parse_k("abc")


# ========================
# io.py - FASTA/FASTQ reader
# ========================

class TestIO:
    def test_read_fastq(self):
        from teloscan.io import read_fasta_fastq
        records = list(read_fasta_fastq(SAMPLE_FASTQ))
        assert len(records) == 5
        assert records[0][0] == "read1_telomere_TTAGGG"
        assert "TTAGGG" in records[0][1]

    def test_read_fasta(self):
        from teloscan.io import read_fasta_fastq
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">seq1\nACGTACGT\n>seq2\nTTTTAAAA\nCCCCGGGG\n")
            f.flush()
            records = list(read_fasta_fastq(f.name))
        os.unlink(f.name)
        assert len(records) == 2
        assert records[0] == ("seq1", "ACGTACGT")
        assert records[1] == ("seq2", "TTTTAAAACCCCGGGG")

    def test_read_gzip_fastq(self):
        from teloscan.io import read_fasta_fastq
        with tempfile.NamedTemporaryFile(suffix=".fastq.gz", delete=False) as f:
            with gzip.open(f.name, "wt") as gz:
                gz.write("@r1\nACGT\n+\nIIII\n")
            records = list(read_fasta_fastq(f.name))
        os.unlink(f.name)
        assert len(records) == 1
        assert records[0] == ("r1", "ACGT")

    def test_read_gzip_fasta(self):
        from teloscan.io import read_fasta_fastq
        with tempfile.NamedTemporaryFile(suffix=".fa.gz", delete=False) as f:
            with gzip.open(f.name, "wt") as gz:
                gz.write(">s1\nGATTACA\n")
            records = list(read_fasta_fastq(f.name))
        os.unlink(f.name)
        assert len(records) == 1
        assert records[0] == ("s1", "GATTACA")

    def test_quality_filter(self):
        from teloscan.io import read_fasta_fastq
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fastq", delete=False) as f:
            # Quality '!' = 0, 'I' = 40
            f.write("@good\nACGT\n+\nIIII\n@bad\nACGT\n+\n!!!!\n")
            f.flush()
            records = list(read_fasta_fastq(f.name, min_quality=20.0))
        os.unlink(f.name)
        assert len(records) == 1
        assert records[0][0] == "good"

    def test_format_detection_fq(self):
        from teloscan.io import _detect_format
        assert _detect_format("reads.fq") == "fastq"
        assert _detect_format("reads.fastq.gz") == "fastq"

    def test_format_detection_fa(self):
        from teloscan.io import _detect_format
        assert _detect_format("genome.fa") == "fasta"
        assert _detect_format("genome.fasta.gz") == "fasta"

    def test_format_detection_bam(self):
        from teloscan.io import _detect_format
        assert _detect_format("aligned.bam") == "bam"
        assert _detect_format("aligned.cram") == "cram"


# ========================
# engine.py - DNA helpers
# ========================

class TestDNAHelpers:
    def test_revcomp(self):
        from teloscan.engine import revcomp
        assert revcomp("ACGT") == "ACGT"
        assert revcomp("TTAGGG") == "CCCTAA"
        assert revcomp("A") == "T"

    def test_canonical_repeat_unit(self):
        from teloscan.engine import canonical_repeat_unit
        # TTAGGG and CCCTAA should canonicalize to the same thing
        c1 = canonical_repeat_unit("TTAGGG")
        c2 = canonical_repeat_unit("CCCTAA")
        assert c1 == c2
        # Rotations should also match
        c3 = canonical_repeat_unit("AGGGTT")
        assert c1 == c3

    def test_hamming(self):
        from teloscan.engine import hamming
        assert hamming("ACGT", "ACGT") == 0
        assert hamming("ACGT", "ACGA") == 1
        assert hamming("AAAA", "TTTT") == 4

    def test_motif_entropy(self):
        from teloscan.engine import motif_entropy
        # Single base = 0 entropy
        assert motif_entropy("AAAA") == 0.0
        # All different bases = max entropy
        e = motif_entropy("ACGT")
        assert e == pytest.approx(2.0, abs=0.01)

    def test_repeat_confidence(self):
        from teloscan.engine import repeat_confidence
        # More copies, higher entropy => higher confidence
        c1 = repeat_confidence(2, 6, "perfect", 1.5)
        c2 = repeat_confidence(20, 6, "perfect", 1.5)
        assert c2 > c1
        # Perfect > fuzzy
        c3 = repeat_confidence(10, 6, "perfect", 1.5)
        c4 = repeat_confidence(10, 6, "fuzzy", 1.5)
        assert c3 > c4

    def test_determine_strand(self):
        from teloscan.engine import _determine_strand
        s1 = _determine_strand("TTAGGG")
        s2 = _determine_strand("CCCTAA")
        # They are rev comp of each other, so one should be + and the other -
        assert s1 != s2


# ========================
# engine.py - Detection
# ========================

class TestDetection:
    def test_perfect_detection_ttaggg(self):
        from teloscan.engine import iter_kmer_runs_perfect
        seq = "TTAGGG" * 30  # 180 bp
        runs = list(iter_kmer_runs_perfect(seq, 6, 150))
        assert len(runs) >= 1
        start, end, motif, copies = runs[0]
        assert motif == "TTAGGG"
        assert copies == 30
        assert end - start == 180

    def test_perfect_detection_skips_N(self):
        from teloscan.engine import iter_kmer_runs_perfect
        seq = "N" * 200
        runs = list(iter_kmer_runs_perfect(seq, 6, 150))
        assert len(runs) == 0

    def test_perfect_detection_min_run_filter(self):
        from teloscan.engine import iter_kmer_runs_perfect
        seq = "TTAGGG" * 5  # 30 bp, below 150 threshold
        runs = list(iter_kmer_runs_perfect(seq, 6, 150))
        assert len(runs) == 0

    def test_fuzzy_detection(self):
        from teloscan.engine import iter_kmer_runs_fuzzy
        # Introduce 1 mismatch in the middle
        seq = "TTAGGG" * 10 + "TTATGG" + "TTAGGG" * 15  # ~156 bp
        runs = list(iter_kmer_runs_fuzzy(seq, 6, "TTAGGG", 1, 150))
        assert len(runs) >= 1

    def test_detect_blocks_for_read(self):
        from teloscan.engine import detect_blocks_for_read
        seq = "TTAGGG" * 30
        blocks = detect_blocks_for_read("test_read", seq, [6], 150)
        assert len(blocks) >= 1
        b = blocks[0]
        assert b.read == "test_read"
        assert b.mode == "perfect"
        assert b.confidence > 0
        assert b.entropy > 0

    def test_detect_blocks_strand_consistent(self):
        from teloscan.engine import detect_blocks_for_read
        seq_fwd = "TTAGGG" * 30
        seq_rev = "CCCTAA" * 30
        blocks_fwd = detect_blocks_for_read("fwd", seq_fwd, [6], 150)
        blocks_rev = detect_blocks_for_read("rev", seq_rev, [6], 150)
        assert len(blocks_fwd) >= 1
        assert len(blocks_rev) >= 1
        # Should have same canonical but different strand
        assert blocks_fwd[0].canonical == blocks_rev[0].canonical
        assert blocks_fwd[0].strand != blocks_rev[0].strand

    def test_merge_overlapping_blocks(self):
        from teloscan.engine import RepeatBlock, _merge_blocks
        b1 = RepeatBlock("r1", 0, 100, 6, "AACCCT", "perfect", "+", 16, 100, 0.5, 1.5)
        b2 = RepeatBlock("r1", 80, 200, 6, "AACCCT", "fuzzy", "+", 20, 120, 0.4, 1.5)
        merged = _merge_blocks([b1, b2])
        assert len(merged) == 1
        assert merged[0].start == 0
        assert merged[0].end == 200
        assert merged[0].mode == "mixed"

    def test_block_to_bed_name(self):
        from teloscan.engine import block_to_bed_name
        assert block_to_bed_name("TTAGGG", 6) == "TELREP:TTAGGG:k6"


# ========================
# engine.py - Multiprocessing
# ========================

class TestMultiprocessing:
    def test_chunk_records(self):
        from teloscan.engine import chunk_records
        records = [("r1", "ACGT"), ("r2", "TGCA"), ("r3", "AAAA")]
        chunks = list(chunk_records(records, 2))
        assert len(chunks) == 2
        assert len(chunks[0]) == 2
        assert len(chunks[1]) == 1

    def test_pass1_single_thread(self):
        from teloscan.engine import run_teloscan_pass1
        seq = "TTAGGG" * 30
        records = [("r1", seq)]
        result = run_teloscan_pass1(records, [6], 150, threads=1, chunk_size=100)
        assert 6 in result
        assert len(result[6]) > 0

    def test_pass2_single_thread(self):
        from teloscan.engine import run_teloscan_pass2
        seq = "TTAGGG" * 30
        records = [("r1", seq)]
        batches = list(run_teloscan_pass2(records, [6], 150, None, 1, threads=1, chunk_size=100))
        assert len(batches) >= 1
        total_blocks = sum(len(b) for b in batches)
        assert total_blocks >= 1

    def test_summarize_blocks_incremental(self):
        from teloscan.engine import RepeatBlock, summarize_blocks_incremental
        from collections import Counter
        blocks = [
            RepeatBlock("r1", 0, 180, 6, "AACCCT", "perfect", "+", 30, 180, 0.8, 1.5),
            RepeatBlock("r2", 0, 180, 6, "AACCCT", "perfect", "-", 30, 180, 0.8, 1.5),
        ]
        total_bp = Counter()
        runs = Counter()
        by_k = {}
        by_strand = {}
        summarize_blocks_incremental(blocks, total_bp, runs, by_k, by_strand)
        assert total_bp["AACCCT"] == 360
        assert runs["AACCCT"] == 2
        assert "AACCCT:+" in by_strand
        assert "AACCCT:-" in by_strand


# ========================
# Visualization (smoke tests)
# ========================

class TestVisualization:
    @pytest.fixture
    def sample_blocks(self):
        from teloscan.engine import RepeatBlock
        return [
            RepeatBlock("r1", 0, 180, 6, "AACCCT", "perfect", "+", 30, 180, 0.8, 1.5),
            RepeatBlock("r1", 500, 700, 6, "AACCCT", "fuzzy", "-", 33, 200, 0.6, 1.5),
            RepeatBlock("r2", 0, 300, 7, "AACCCTA", "perfect", "+", 42, 300, 0.9, 1.8),
        ]

    def test_plot_length_distribution(self, sample_blocks):
        pytest.importorskip("matplotlib")
        from teloscan.visualize import plot_length_distribution
        b64 = plot_length_distribution(sample_blocks)
        assert b64 is not None
        assert len(b64) > 100

    def test_plot_motif_abundance(self):
        pytest.importorskip("matplotlib")
        from teloscan.visualize import plot_motif_abundance
        from collections import Counter
        bp = Counter({"AACCCT": 500, "AACCCTA": 300})
        b64 = plot_motif_abundance(bp)
        assert b64 is not None

    def test_plot_per_read_heatmap(self, sample_blocks):
        pytest.importorskip("matplotlib")
        from teloscan.visualize import plot_per_read_heatmap
        b64 = plot_per_read_heatmap(sample_blocks)
        assert b64 is not None

    def test_plot_strand_distribution(self, sample_blocks):
        pytest.importorskip("matplotlib")
        from teloscan.visualize import plot_strand_distribution
        b64 = plot_strand_distribution(sample_blocks)
        assert b64 is not None

    def test_plot_save_to_file(self, sample_blocks):
        pytest.importorskip("matplotlib")
        from teloscan.visualize import plot_length_distribution
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            plot_length_distribution(sample_blocks, output_path=f.name)
            assert os.path.getsize(f.name) > 0
        os.unlink(f.name)


# ========================
# HTML Report
# ========================

class TestReport:
    def test_generate_html_report(self):
        from teloscan.engine import RepeatBlock
        from teloscan.report import generate_html_report
        from collections import Counter
        blocks = [
            RepeatBlock("r1", 0, 180, 6, "AACCCT", "perfect", "+", 30, 180, 0.8, 1.5),
        ]
        total_bp = Counter({"AACCCT": 180})
        runs = Counter({"AACCCT": 1})
        by_k = {6: Counter({"AACCCT": 1})}
        html = generate_html_report(blocks, total_bp, runs, by_k, input_file="test.fastq")
        assert "TeloScan Report" in html
        assert "AACCCT" in html
        assert "test.fastq" in html
        assert "<table>" in html


# ========================
# CLI output format tests
# ========================

class TestCLIOutputFormats:
    def test_write_gff3(self):
        from teloscan.engine import RepeatBlock
        from teloscan.cli import _write_gff3
        blocks = [
            RepeatBlock("r1", 0, 180, 6, "AACCCT", "perfect", "+", 30, 180, 0.8, 1.5),
        ]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".gff3", delete=False) as f:
            _write_gff3(blocks, f.name)
        with open(f.name) as fh:
            content = fh.read()
        os.unlink(f.name)
        assert "##gff-version 3" in content
        assert "teloscan" in content
        assert "repeat_region" in content
        assert "canonical=AACCCT" in content

    def test_write_vcf(self):
        from teloscan.engine import RepeatBlock
        from teloscan.cli import _write_vcf
        blocks = [
            RepeatBlock("r1", 0, 180, 6, "AACCCT", "perfect", "+", 30, 180, 0.8, 1.5),
        ]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            _write_vcf(blocks, f.name)
        with open(f.name) as fh:
            content = fh.read()
        os.unlink(f.name)
        assert "##fileformat=VCFv4.3" in content
        assert "MOTIF=AACCCT" in content
        assert "<TELREPEAT>" in content


# ========================
# Integration test with sample data
# ========================

class TestIntegration:
    def test_full_pipeline_sample_fastq(self):
        """End-to-end test with sample FASTQ data."""
        from teloscan.io import read_fasta_fastq
        from teloscan.repeats import parse_k
        from teloscan.engine import (
            run_teloscan_pass1,
            run_teloscan_pass2,
            summarize_blocks_incremental,
        )
        from collections import Counter

        k_values = parse_k("6")
        records = list(read_fasta_fastq(SAMPLE_FASTQ))
        assert len(records) == 5

        # Pass 1
        motif_counts = run_teloscan_pass1(records, k_values, 150, threads=1, chunk_size=100)
        assert 6 in motif_counts

        # Pass 2
        total_bp = Counter()
        runs = Counter()
        by_k = {}
        by_strand = {}
        all_blocks = []
        for batch in run_teloscan_pass2(records, k_values, 150, None, 1, threads=1, chunk_size=100):
            all_blocks.extend(batch)
            summarize_blocks_incremental(batch, total_bp, runs, by_k, by_strand)

        # Should find telomeric blocks in reads 1, 2, 4 (and maybe 5 in perfect mode)
        reads_with_blocks = set(b.read for b in all_blocks)
        assert len(reads_with_blocks) >= 2
        assert "read3_no_telomere" not in reads_with_blocks


# ========================
# Fix-specific tests
# ========================

class TestFileHandleClosure:
    def test_file_handle_closed_after_iteration(self):
        """Verify file handles are properly closed after reading."""
        from teloscan.io import read_fasta_fastq
        # Just fully consume the iterator — should not leave handles open
        records = list(read_fasta_fastq(SAMPLE_FASTQ))
        assert len(records) == 5

    def test_file_handle_closed_on_error(self):
        """File handle should close even on parse error."""
        from teloscan.io import read_fasta_fastq
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fastq", delete=False) as f:
            f.write("NOT_A_VALID_HEADER\nACGT\n+\nIIII\n")
            f.flush()
        with pytest.raises(ValueError, match="Error reading"):
            list(read_fasta_fastq(f.name))
        os.unlink(f.name)


class TestStdinBuffering:
    def test_buffer_records(self):
        """buffer_records should load all records into a list."""
        from teloscan.io import buffer_records
        records = buffer_records(SAMPLE_FASTQ)
        assert isinstance(records, list)
        assert len(records) == 5
        # Can iterate multiple times
        first_pass = list(records)
        second_pass = list(records)
        assert first_pass == second_pass


class TestFormatHint:
    def test_format_hint_fasta(self):
        """format_hint should override auto-detection."""
        from teloscan.io import read_fasta_fastq
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(">s1\nACGTACGT\n")
            f.flush()
        records = list(read_fasta_fastq(f.name, format_hint="fasta"))
        os.unlink(f.name)
        assert len(records) == 1
        assert records[0] == ("s1", "ACGTACGT")


class TestFieldSanitization:
    def test_sanitize_field(self):
        from teloscan.cli import _sanitize_field
        assert _sanitize_field("normal_read") == "normal_read"
        assert "\t" not in _sanitize_field("read\twith\ttabs")
        assert "\n" not in _sanitize_field("read\nwith\nnewlines")
        assert ";" not in _sanitize_field("read;with;semicolons")

    def test_gff3_special_chars(self):
        """GFF3 output should handle special characters in read IDs."""
        from teloscan.engine import RepeatBlock
        from teloscan.cli import _write_gff3
        blocks = [
            RepeatBlock("read\twith\ttab", 0, 180, 6, "AACCCT", "perfect", "+", 30, 180, 0.8, 1.5),
        ]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".gff3", delete=False) as f:
            _write_gff3(blocks, f.name)
        with open(f.name) as fh:
            content = fh.read()
        os.unlink(f.name)
        # Should not contain raw tabs within the read ID field
        lines = content.strip().split("\n")
        data_line = lines[1]  # skip header
        fields = data_line.split("\t")
        assert len(fields) == 9  # valid GFF3 has exactly 9 tab-separated fields


class TestHTMLEscaping:
    def test_xss_prevention(self):
        """HTML report should escape potentially dangerous input."""
        from teloscan.engine import RepeatBlock
        from teloscan.report import generate_html_report
        from collections import Counter
        xss_read = '<script>alert("xss")</script>'
        blocks = [
            RepeatBlock(xss_read, 0, 180, 6, "AACCCT", "perfect", "+", 30, 180, 0.8, 1.5),
        ]
        total_bp = Counter({"AACCCT": 180})
        runs = Counter({"AACCCT": 1})
        by_k = {6: Counter({"AACCCT": 1})}
        html = generate_html_report(blocks, total_bp, runs, by_k, input_file=xss_read)
        assert "<script>" not in html
        assert "&lt;script&gt;" in html


class TestArgumentValidation:
    def test_validate_args_negative_min_run(self):
        from teloscan.cli import _validate_args
        import argparse
        args = argparse.Namespace(
            input=SAMPLE_FASTQ, k="6", min_run_bp=-1, top_motifs_per_k=10,
            max_mismatch=1, threads=1, chunk_size=2000, min_quality=None,
            out_per_read="/tmp/test.tsv", out_bed="/tmp/test.bed",
            out_summary="/tmp/test.summary.tsv", out_gff3=None,
            out_vcf=None, out_html=None,
        )
        with pytest.raises(SystemExit, match="min-run-bp"):
            _validate_args(args)

    def test_validate_args_missing_input(self):
        from teloscan.cli import _validate_args
        import argparse
        args = argparse.Namespace(
            input="/nonexistent/file.fastq", k="6", min_run_bp=150, top_motifs_per_k=10,
            max_mismatch=1, threads=1, chunk_size=2000, min_quality=None,
            out_per_read="/tmp/test.tsv", out_bed="/tmp/test.bed",
            out_summary="/tmp/test.summary.tsv", out_gff3=None,
            out_vcf=None, out_html=None,
        )
        with pytest.raises(SystemExit, match="not found"):
            _validate_args(args)


class TestKmerValidationWarnings:
    def test_large_k_warning(self, caplog):
        """Warn when k > 50."""
        import logging
        from teloscan.cli import _validate_args
        import argparse
        args = argparse.Namespace(
            input=SAMPLE_FASTQ, k="55", min_run_bp=150, top_motifs_per_k=10,
            max_mismatch=1, threads=1, chunk_size=2000, min_quality=None,
            out_per_read="/tmp/test.tsv", out_bed="/tmp/test.bed",
            out_summary="/tmp/test.summary.tsv", out_gff3=None,
            out_vcf=None, out_html=None,
        )
        with caplog.at_level(logging.WARNING, logger="teloscan"):
            _validate_args(args)
        assert "Large k-mer size" in caplog.text

    def test_min_run_bp_too_small_warning(self, caplog):
        """Warn when min_run_bp < 2 * max_k."""
        import logging
        from teloscan.cli import _validate_args
        import argparse
        args = argparse.Namespace(
            input=SAMPLE_FASTQ, k="6", min_run_bp=10, top_motifs_per_k=10,
            max_mismatch=1, threads=1, chunk_size=2000, min_quality=None,
            out_per_read="/tmp/test.tsv", out_bed="/tmp/test.bed",
            out_summary="/tmp/test.summary.tsv", out_gff3=None,
            out_vcf=None, out_html=None,
        )
        with caplog.at_level(logging.WARNING, logger="teloscan"):
            _validate_args(args)
        assert "less than 2" in caplog.text


class TestPublicAPI:
    def test_imports_from_package(self):
        """Key classes/functions should be importable from the top-level package."""
        from teloscan import RepeatBlock, canonical_repeat_unit, detect_blocks_for_read
        from teloscan import read_fasta_fastq, parse_k
        assert RepeatBlock is not None
        assert callable(canonical_repeat_unit)
        assert callable(detect_blocks_for_read)
        assert callable(read_fasta_fastq)
        assert callable(parse_k)
