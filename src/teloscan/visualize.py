"""Visualization utilities for TeloScan results."""

from __future__ import annotations

import io
import base64
from collections import Counter, defaultdict
from typing import Dict, List, Optional

from teloscan.engine import RepeatBlock


def _check_matplotlib():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except ImportError:
        raise ImportError(
            "matplotlib is required for visualization. Install with: pip install matplotlib"
        )


def plot_length_distribution(
    blocks: List[RepeatBlock],
    output_path: Optional[str] = None,
    title: str = "Telomere Block Length Distribution",
) -> Optional[str]:
    """
    Plot histogram of telomere block lengths.
    Returns base64-encoded PNG if output_path is None.
    """
    plt = _check_matplotlib()

    lengths = [b.run_bp for b in blocks]
    if not lengths:
        return None

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(lengths, bins=min(50, len(set(lengths))), color="#2196F3", edgecolor="white", alpha=0.8)
    ax.set_xlabel("Block length (bp)")
    ax.set_ylabel("Count")
    ax.set_title(title)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return None
    else:
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        return base64.b64encode(buf.read()).decode("ascii")


def plot_motif_abundance(
    total_bp_by_motif: Counter,
    output_path: Optional[str] = None,
    top_n: int = 20,
    title: str = "Top Telomeric Motifs by Total BP",
) -> Optional[str]:
    """Bar chart of top motifs by total base pairs covered."""
    plt = _check_matplotlib()

    if not total_bp_by_motif:
        return None

    items = total_bp_by_motif.most_common(top_n)
    labels = [m for m, _ in items]
    values = [v for _, v in items]

    fig, ax = plt.subplots(figsize=(12, 6))
    bars = ax.barh(range(len(labels)), values, color="#4CAF50", edgecolor="white")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8, fontfamily="monospace")
    ax.set_xlabel("Total base pairs")
    ax.set_title(title)
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return None
    else:
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        return base64.b64encode(buf.read()).decode("ascii")


def plot_per_read_heatmap(
    blocks: List[RepeatBlock],
    output_path: Optional[str] = None,
    max_reads: int = 50,
    title: str = "Telomeric Block Locations per Read",
) -> Optional[str]:
    """
    Heatmap showing block positions within reads.
    Each row is a read, x-axis is position, colored bars represent blocks.
    """
    plt = _check_matplotlib()
    from matplotlib.patches import Rectangle

    if not blocks:
        return None

    # Group blocks by read, pick top reads by total bp
    by_read = defaultdict(list)
    for b in blocks:
        by_read[b.read].append(b)

    read_totals = {r: sum(b.run_bp for b in bs) for r, bs in by_read.items()}
    top_reads = sorted(read_totals, key=read_totals.get, reverse=True)[:max_reads]

    if not top_reads:
        return None

    # Determine max position for x-axis
    max_pos = max(b.end for b in blocks if b.read in set(top_reads))

    fig, ax = plt.subplots(figsize=(14, max(4, len(top_reads) * 0.3)))
    colors = {"perfect": "#2196F3", "fuzzy": "#FF9800", "mixed": "#9C27B0"}

    for idx, read_id in enumerate(top_reads):
        for b in by_read[read_id]:
            color = colors.get(b.mode, "#607D8B")
            rect = Rectangle((b.start, idx - 0.35), b.run_bp, 0.7, color=color, alpha=0.7)
            ax.add_patch(rect)

    ax.set_xlim(0, max_pos)
    ax.set_ylim(-0.5, len(top_reads) - 0.5)
    ax.set_yticks(range(len(top_reads)))
    ax.set_yticklabels([r[:20] for r in top_reads], fontsize=7)
    ax.set_xlabel("Position (bp)")
    ax.set_title(title)
    ax.invert_yaxis()

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=c, label=m) for m, c in colors.items()]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=8)

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return None
    else:
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        return base64.b64encode(buf.read()).decode("ascii")


def plot_strand_distribution(
    blocks: List[RepeatBlock],
    output_path: Optional[str] = None,
    title: str = "Strand Distribution of Telomeric Blocks",
) -> Optional[str]:
    """Pie chart of strand distribution."""
    plt = _check_matplotlib()

    if not blocks:
        return None

    strand_counts = Counter(b.strand for b in blocks)

    fig, ax = plt.subplots(figsize=(6, 6))
    labels = list(strand_counts.keys())
    sizes = list(strand_counts.values())
    colors = ["#2196F3", "#FF5722"][:len(labels)]
    ax.pie(sizes, labels=labels, colors=colors, autopct="%1.1f%%", startangle=90)
    ax.set_title(title)

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return None
    else:
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        return base64.b64encode(buf.read()).decode("ascii")
