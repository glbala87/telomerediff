"""FASTA / FASTQ / BAM / CRAM reader with gzip support."""

from __future__ import annotations

import gzip
import logging
import os
import sys
from typing import IO, Iterator, List, Optional, Tuple

logger = logging.getLogger("teloscan")


def _open_text(path: str) -> IO[str]:
    """Open a plain or gzip-compressed text file."""
    if path == "-":
        return sys.stdin
    if path.endswith(".gz") or path.endswith(".gzip"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def _read_fastq(fh: IO[str]) -> Iterator[Tuple[str, str, Optional[str]]]:
    """Yield (read_id, sequence, quality) from a FASTQ stream."""
    while True:
        header = fh.readline()
        if not header:
            return
        header = header.rstrip("\n\r")
        if not header:
            continue
        if not header.startswith("@"):
            raise ValueError(f"Expected FASTQ header starting with '@', got: {header!r}")
        read_id = header[1:].split()[0]
        seq = fh.readline().rstrip("\n\r")
        plus = fh.readline().rstrip("\n\r")
        if not plus.startswith("+"):
            raise ValueError(f"Expected '+' line, got: {plus!r}")
        qual = fh.readline().rstrip("\n\r")
        yield read_id, seq, qual


def _read_fasta(fh: IO[str]) -> Iterator[Tuple[str, str, Optional[str]]]:
    """Yield (read_id, sequence, None) from a FASTA stream."""
    read_id = None
    seq_parts: List[str] = []
    for line in fh:
        line = line.rstrip("\n\r")
        if not line:
            continue
        if line.startswith(">"):
            if read_id is not None:
                yield read_id, "".join(seq_parts), None
            read_id = line[1:].split()[0]
            seq_parts = []
        else:
            seq_parts.append(line)
    if read_id is not None:
        yield read_id, "".join(seq_parts), None


def _detect_format(path: str) -> str:
    """Detect file format from extension or first character."""
    lower = path.lower().replace(".gz", "").replace(".gzip", "")
    if lower.endswith((".fq", ".fastq")):
        return "fastq"
    if lower.endswith((".fa", ".fasta", ".fna")):
        return "fasta"
    if lower.endswith(".bam"):
        return "bam"
    if lower.endswith(".cram"):
        return "cram"
    # Peek at first byte for regular files
    if path != "-" and os.path.isfile(path):
        opener = gzip.open if path.endswith((".gz", ".gzip")) else open
        mode = "rt" if path.endswith((".gz", ".gzip")) else "r"
        try:
            with opener(path, mode) as fh:
                first = fh.read(1)
                if first == "@":
                    return "fastq"
                if first == ">":
                    return "fasta"
        except Exception:
            pass
    return "fastq"  # default


def _read_bam_cram(path: str, reference: Optional[str] = None) -> Iterator[Tuple[str, str, Optional[str]]]:
    """Read sequences from BAM/CRAM files using pysam."""
    try:
        import pysam
    except ImportError:
        raise ImportError(
            "pysam is required for BAM/CRAM input. Install with: pip install pysam"
        )
    kwargs = {}
    if reference:
        kwargs["reference_filename"] = reference
    with pysam.AlignmentFile(path, "rb", **kwargs) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            seq = read.query_sequence
            if seq:
                qual_str = None
                if read.query_qualities is not None:
                    qual_str = "".join(chr(q + 33) for q in read.query_qualities)
                yield read.query_name, seq, qual_str


def read_fasta_fastq(
    path: str,
    reference: Optional[str] = None,
    min_quality: Optional[float] = None,
    format_hint: Optional[str] = None,
) -> Iterator[Tuple[str, str]]:
    """
    Universal sequence reader.

    Supports FASTA, FASTQ (plain or gzip), BAM, and CRAM.
    Yields (read_id, sequence) tuples.

    Parameters
    ----------
    path : str
        Input file path, or '-' for stdin.
    reference : str, optional
        Reference FASTA for CRAM decoding.
    min_quality : float, optional
        Minimum average base quality to keep a read (FASTQ/BAM only).
    format_hint : str, optional
        Force format: 'fasta', 'fastq', 'bam', 'cram'. Useful for stdin.
    """
    fmt = format_hint if format_hint else _detect_format(path)

    if fmt in ("bam", "cram"):
        if path == "-":
            raise ValueError("BAM/CRAM input cannot be read from stdin. Provide a file path.")
        for item in _read_bam_cram(path, reference=reference):
            read_id, seq, qual = item
            if min_quality is not None and qual and len(qual) > 0:
                avg_q = sum(ord(c) - 33 for c in qual) / len(qual)
                if avg_q < min_quality:
                    continue
            yield read_id, seq
        return

    fh = _open_text(path)
    is_stdin = (path == "-")
    try:
        if fmt == "fasta":
            records = _read_fasta(fh)
        else:
            records = _read_fastq(fh)

        for item in records:
            read_id, seq, qual = item
            if min_quality is not None and qual and len(qual) > 0:
                avg_q = sum(ord(c) - 33 for c in qual) / len(qual)
                if avg_q < min_quality:
                    continue
            yield read_id, seq
    except ValueError as e:
        raise ValueError(f"Error reading {path!r} as {fmt}: {e}") from e
    finally:
        if not is_stdin:
            fh.close()


def buffer_records(path: str, **kwargs) -> List[Tuple[str, str]]:
    """
    Read all records into memory. Required for stdin or when multiple passes are needed.

    Parameters are the same as read_fasta_fastq().
    """
    return list(read_fasta_fastq(path, **kwargs))
