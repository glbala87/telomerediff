from __future__ import annotations

import math
from dataclasses import dataclass, field
from collections import Counter
from typing import Dict, Iterable, Iterator, List, Optional, Tuple
import multiprocessing as mp

# -------------------------
# Basic DNA helpers
# -------------------------
_DNA_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_DNA_COMP)[::-1]


def rotations(s: str) -> Iterator[str]:
    for i in range(len(s)):
        yield s[i:] + s[:i]


def canonical_repeat_unit(motif: str) -> str:
    """
    Canonicalize a repeat unit by selecting the lexicographically smallest string among:
      - all rotations of motif
      - all rotations of reverse-complement(motif)
    """
    motif = motif.upper()
    rc = revcomp(motif).upper()
    return min(min(rotations(motif)), min(rotations(rc)))


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def motif_entropy(motif: str) -> float:
    """Shannon entropy of base composition (bits). Higher = more complex."""
    n = len(motif)
    if n == 0:
        return 0.0
    counts = Counter(motif.upper())
    return -sum((c / n) * math.log2(c / n) for c in counts.values() if c > 0)


def repeat_confidence(copies: int, k: int, mode: str, entropy: float) -> float:
    """
    Confidence score [0, 1] for a repeat block.

    Combines:
      - copy count (more copies = higher confidence)
      - sequence complexity (entropy; low-complexity single-base repeats score lower)
      - detection mode (perfect > fuzzy)
    """
    # Copy contribution: asymptotic to 1
    copy_score = 1.0 - 1.0 / (1.0 + copies * 0.2)

    # Entropy: max for DNA is 2.0 bits; normalize
    max_entropy = 2.0
    entropy_score = min(entropy / max_entropy, 1.0)

    # Mode: perfect is more reliable
    mode_weight = 1.0 if mode == "perfect" else 0.85 if mode == "fuzzy" else 0.9

    return round(copy_score * entropy_score * mode_weight, 4)


def _determine_strand(motif: str) -> str:
    """Determine strand by comparing motif to its reverse complement canonically."""
    motif_u = motif.upper()
    rc = revcomp(motif_u)
    fwd_min = min(rotations(motif_u))
    rev_min = min(rotations(rc))
    if fwd_min <= rev_min:
        return "+"
    return "-"


# -------------------------
# Repeat block model
# -------------------------
@dataclass
class RepeatBlock:
    read: str
    start: int
    end: int
    k: int
    canonical: str
    mode: str        # perfect / fuzzy / mixed
    strand: str      # + or -
    copies: int
    run_bp: int
    confidence: float = 0.0
    entropy: float = 0.0


# -------------------------
# Core detection routines
# -------------------------
def iter_kmer_runs_perfect(seq: str, k: int, min_run_bp: int) -> Iterator[Tuple[int, int, str, int]]:
    """
    Find in-frame consecutive runs of identical k-mers.
    Returns tuples: (start, end, motif, copies)
    """
    seq = seq.upper()
    n = len(seq)
    if n < 2 * k:
        return

    i = 0
    while i + k <= n:
        motif = seq[i:i + k]
        if "N" in motif:
            i += 1
            continue

        j = i + k
        copies = 1
        while j + k <= n and seq[j:j + k] == motif:
            copies += 1
            j += k

        run_len = copies * k
        if copies >= 2 and run_len >= min_run_bp:
            yield i, j, motif, copies

        i = j if copies >= 2 else i + 1


def iter_kmer_runs_fuzzy(seq: str, k: int, seed: str, max_mismatch: int, min_run_bp: int) -> Iterator[Tuple[int, int, int]]:
    """
    Find in-frame runs matching a seed motif allowing up to max_mismatch per k-mer copy.
    Returns tuples: (start, end, copies)
    """
    seq = seq.upper()
    seed = seed.upper()
    n = len(seq)
    if n < 2 * k:
        return

    i = 0
    while i + k <= n:
        chunk = seq[i:i + k]
        if "N" in chunk:
            i += 1
            continue

        if hamming(chunk, seed) <= max_mismatch:
            j = i + k
            copies = 1
            while j + k <= n:
                nxt = seq[j:j + k]
                if "N" in nxt:
                    break
                if hamming(nxt, seed) <= max_mismatch:
                    copies += 1
                    j += k
                else:
                    break

            run_bp = copies * k
            if copies >= 2 and run_bp >= min_run_bp:
                yield i, j, copies

            i = j
        else:
            i += 1


def _merge_blocks(blocks: List[RepeatBlock]) -> List[RepeatBlock]:
    """Merge overlapping/adjacent blocks with same (read, canonical, k)."""
    if not blocks:
        return blocks

    blocks.sort(key=lambda b: (b.read, b.k, b.canonical, b.start, b.end))
    merged: List[RepeatBlock] = [blocks[0]]

    for b in blocks[1:]:
        prev = merged[-1]
        same_key = (b.read == prev.read and b.k == prev.k and b.canonical == prev.canonical)
        if same_key and b.start <= prev.end:
            prev.end = max(prev.end, b.end)
            prev.run_bp = prev.end - prev.start
            prev.copies = prev.run_bp // prev.k
            if prev.mode != b.mode:
                prev.mode = "mixed"
            # Recompute confidence after merge
            prev.entropy = motif_entropy(prev.canonical)
            prev.confidence = repeat_confidence(prev.copies, prev.k, prev.mode, prev.entropy)
        else:
            merged.append(b)

    return merged


def detect_blocks_for_read(
    read_id: str,
    seq: str,
    k_values: List[int],
    min_run_bp: int,
    fuzzy_seeds: Optional[Dict[int, List[str]]] = None,
    max_mismatch: int = 1,
) -> List[RepeatBlock]:
    """
    Detect repeat blocks for a single read.
    - perfect pass for all k (with proper strand detection)
    - optional fuzzy pass for provided seeds per k
    """
    seq_u = seq.upper()
    blocks: List[RepeatBlock] = []

    # Perfect discovery (with strand detection)
    for k in k_values:
        for start, end, motif, copies in iter_kmer_runs_perfect(seq_u, k, min_run_bp):
            canon = canonical_repeat_unit(motif)
            strand = _determine_strand(motif)
            ent = motif_entropy(canon)
            conf = repeat_confidence(copies, k, "perfect", ent)
            blocks.append(
                RepeatBlock(
                    read=read_id,
                    start=start,
                    end=end,
                    k=k,
                    canonical=canon,
                    mode="perfect",
                    strand=strand,
                    copies=copies,
                    run_bp=copies * k,
                    confidence=conf,
                    entropy=ent,
                )
            )

    # Fuzzy refine for top seeds
    if fuzzy_seeds:
        for k, seeds in fuzzy_seeds.items():
            if k not in k_values:
                continue
            for canon in seeds:
                seed = canon
                for seed_str, strand in [(seed, "+"), (revcomp(seed), "-")]:
                    for start, end, copies in iter_kmer_runs_fuzzy(seq_u, k, seed_str, max_mismatch, min_run_bp):
                        ent = motif_entropy(canon)
                        conf = repeat_confidence(copies, k, "fuzzy", ent)
                        blocks.append(
                            RepeatBlock(
                                read=read_id,
                                start=start,
                                end=end,
                                k=k,
                                canonical=canon,
                                mode="fuzzy",
                                strand=strand,
                                copies=copies,
                                run_bp=copies * k,
                                confidence=conf,
                                entropy=ent,
                            )
                        )

    return _merge_blocks(blocks)


def block_to_bed_name(canonical: str, k: int) -> str:
    return f"TELREP:{canonical}:k{k}"


# -------------------------
# Chunking + multiprocessing helpers
# -------------------------
def chunk_records(records: Iterable[Tuple[str, str]], chunk_size: int) -> Iterator[List[Tuple[str, str]]]:
    chunk: List[Tuple[str, str]] = []
    for r in records:
        chunk.append(r)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def _pass1_worker(args):
    """Pass 1 worker: perfect discovery counts only (no block outputs)."""
    chunk, k_values, min_run_bp = args
    out = {k: Counter() for k in k_values}
    for read_id, seq in chunk:
        blocks = detect_blocks_for_read(read_id, seq, k_values, min_run_bp, fuzzy_seeds=None)
        for b in blocks:
            out[b.k][b.canonical] += 1
    return out


def _pass2_worker(args):
    """Pass 2 worker: produce blocks (perfect + fuzzy seeds)."""
    chunk, k_values, min_run_bp, fuzzy_seeds, max_mismatch = args
    out_blocks: List[RepeatBlock] = []
    for read_id, seq in chunk:
        out_blocks.extend(
            detect_blocks_for_read(
                read_id,
                seq,
                k_values,
                min_run_bp,
                fuzzy_seeds=fuzzy_seeds,
                max_mismatch=max_mismatch,
            )
        )
    return out_blocks


def run_teloscan_pass1(
    records: Iterable[Tuple[str, str]],
    k_values: List[int],
    min_run_bp: int,
    threads: int,
    chunk_size: int,
) -> Dict[int, Counter]:
    """
    Streaming pass 1: discover perfect blocks, count canonical motifs per k.
    """
    motif_counts_by_k: Dict[int, Counter] = {k: Counter() for k in k_values}
    work_iter = ((chunk, k_values, min_run_bp) for chunk in chunk_records(records, chunk_size))

    if threads > 1:
        with mp.Pool(processes=threads) as pool:
            for local in pool.imap_unordered(_pass1_worker, work_iter, chunksize=1):
                for k in k_values:
                    motif_counts_by_k[k].update(local[k])
    else:
        for args in work_iter:
            local = _pass1_worker(args)
            for k in k_values:
                motif_counts_by_k[k].update(local[k])

    return motif_counts_by_k


def run_teloscan_pass2(
    records: Iterable[Tuple[str, str]],
    k_values: List[int],
    min_run_bp: int,
    fuzzy_seeds: Optional[Dict[int, List[str]]],
    max_mismatch: int,
    threads: int,
    chunk_size: int,
) -> Iterator[List[RepeatBlock]]:
    """
    Streaming pass 2: yield blocks per chunk (perfect + fuzzy seeds).
    """
    work_iter = (
        (chunk, k_values, min_run_bp, fuzzy_seeds, max_mismatch)
        for chunk in chunk_records(records, chunk_size)
    )

    if threads > 1:
        with mp.Pool(processes=threads) as pool:
            for blocks in pool.imap_unordered(_pass2_worker, work_iter, chunksize=1):
                yield blocks
    else:
        for args in work_iter:
            yield _pass2_worker(args)


def summarize_blocks_incremental(
    blocks: List[RepeatBlock],
    total_bp_by_motif: Counter,
    runs_by_motif: Counter,
    by_k: Dict[int, Counter],
    by_strand: Optional[Dict[str, Counter]] = None,
) -> None:
    """Update summary counters in-place for a batch of blocks."""
    for b in blocks:
        total_bp_by_motif[b.canonical] += b.run_bp
        runs_by_motif[b.canonical] += 1
        if b.k not in by_k:
            by_k[b.k] = Counter()
        by_k[b.k][b.canonical] += 1
        if by_strand is not None:
            key = f"{b.canonical}:{b.strand}"
            by_strand[key] = by_strand.get(key, Counter())
            by_strand[key]["runs"] += 1
            by_strand[key]["bp"] += b.run_bp
