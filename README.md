# telomerediff

**telomerediff** is an open-source Python tool for **de novo discovery of telomeric repeat blocks**
directly from **long-read sequencing data (Oxford Nanopore / PacBio)** when the telomeric motif is unknown.

Unlike reference-based approaches, telomerediff operates in an **alignment-free** manner, enabling
robust telomere discovery in non-model organisms and incomplete genome assemblies.

---

## Key Features

- **Alignment-free** telomeric repeat discovery — no reference genome required
- **Multi-format input** — FASTQ, FASTA, BAM, CRAM (plain or gzip-compressed), stdin
- **De novo k-mer detection** with configurable k-mer sizes (single, range, or mixed)
- **Error-tolerant refinement** (fuzzy matching) for noisy ONT reads
- **Strand-aware** detection with proper forward/reverse complement handling
- **Confidence scoring** combining copy count, sequence entropy, and detection mode
- **Quality filtering** — filter reads by average base quality
- **Multi-threaded** processing with configurable thread count and chunk size
- **Rich output formats** — TSV, BED, GFF3, VCF, self-contained HTML report, PNG plots
- **Library API** — use telomerediff programmatically in Python scripts

---

## Installation

```bash
git clone https://github.com/glbala87/telomerediff.git
cd telomerediff
pip install -e .
```

### Optional dependencies

```bash
# Progress bar
pip install -e ".[progress]"

# Visualization (plots + HTML report)
pip install -e ".[plots]"

# BAM/CRAM support
pip install -e ".[bam]"

# Everything
pip install -e ".[all]"

# Development (testing)
pip install -e ".[dev]"
```

---

## Quick Start

```bash
# Basic run
telomerediff -i reads.fastq.gz

# With error-tolerant refinement (recommended for ONT)
telomerediff -i reads.fastq.gz --refine

# Full analysis with all outputs
telomerediff -i reads.fastq.gz \
  --refine \
  --out-html report.html \
  --out-gff3 telomeres.gff3 \
  --out-vcf telomeres.vcf \
  --out-plots plots/
```

### Try with sample data

```bash
telomerediff -i tests/data/sample.fastq -k 6-7 --refine --out-html report.html
```

---

## Usage

```
telomerediff -i <input> [options]
```

### Input options

| Flag | Description |
|------|-------------|
| `-i`, `--input` | Input file: FASTQ, FASTA, BAM, CRAM (.gz supported), or `-` for stdin |
| `--reference` | Reference FASTA for CRAM decoding |
| `--format` | Force input format: `fasta`, `fastq`, `bam`, `cram` (default: auto-detect) |
| `--min-quality` | Minimum average base quality to keep a read |

### Algorithm options

| Flag | Default | Description |
|------|---------|-------------|
| `-k`, `--k` | `4-15` | k-mer sizes: `6`, `4-15`, or `5,6,7` |
| `--min-run-bp` | `150` | Minimum repeat block length in bp |
| `--refine` | off | Enable error-tolerant fuzzy matching pass |
| `--max-mismatch` | `1` | Max mismatches per k-mer copy (fuzzy mode) |
| `--top-motifs-per-k` | `10` | Top motifs per k used as fuzzy seeds |

### Performance options

| Flag | Default | Description |
|------|---------|-------------|
| `-t`, `--threads` | CPU/2 | Number of worker threads |
| `--chunk-size` | `2000` | Reads per processing chunk |
| `--no-progress` | off | Disable progress bar |
| `--verbose` | off | Enable debug logging |

### Output options

| Flag | Default | Description |
|------|---------|-------------|
| `--out-per-read` | `telomerediff.per_read.tsv` | Per-read blocks TSV |
| `--out-bed` | `telomerediff.blocks.bed` | BED format coordinates |
| `--out-summary` | `telomerediff.summary.tsv` | Motif-level summary TSV |
| `--out-gff3` | *(off)* | GFF3 annotation file |
| `--out-vcf` | *(off)* | VCF annotation file |
| `--out-html` | *(off)* | Self-contained HTML report with embedded plots |
| `--out-plots` | *(off)* | Directory for PNG plot files |

---

## Output Files

### per_read.tsv

One row per detected telomeric block:

| Column | Description |
|--------|-------------|
| `read` | Read identifier |
| `start` | Block start position (0-based) |
| `end` | Block end position |
| `k` | k-mer size |
| `canonical` | Canonical (rotation + strand normalized) repeat motif |
| `mode` | Detection mode: `perfect`, `fuzzy`, or `mixed` |
| `strand` | `+` or `-` |
| `copies` | Number of repeat copies |
| `run_bp` | Total block length in base pairs |
| `confidence` | Confidence score (0–1) |
| `entropy` | Shannon entropy of motif base composition |

### summary.tsv

Aggregated statistics per motif with strand breakdown:

| Column | Description |
|--------|-------------|
| `canonical` | Canonical repeat motif |
| `k` | k-mer size |
| `observed_runs` | Number of blocks found |
| `total_bp` | Total base pairs across all blocks |
| `strand_plus_runs` | Blocks on + strand |
| `strand_minus_runs` | Blocks on - strand |

### blocks.bed

Standard BED6 format for genome browsers (IGV, UCSC).

### GFF3 / VCF

Standard annotation formats for integration with genomics pipelines. Include confidence scores, copy counts, and motif information in attributes/INFO fields.

### HTML Report

Self-contained HTML file (open in any browser) with:
- Summary statistics cards
- Block length distribution histogram
- Top motif abundance bar chart
- Per-read block location heatmap
- Strand distribution pie chart
- Sortable motif and block tables

---

## How It Works

telomerediff uses a **two-pass algorithm**:

**Pass 1 — Motif Discovery:**
Scans all reads for in-frame consecutive identical k-mers. Counts canonical motif occurrences across all reads to identify the most abundant telomeric repeat candidates.

**Pass 2 — Block Detection:**
Re-scans reads using the top motifs from Pass 1. In perfect mode, finds exact tandem repeats. With `--refine`, additionally performs fuzzy matching (Hamming distance) to recover blocks with sequencing errors. Outputs are written incrementally for memory efficiency.

**Canonicalization:**
Repeat motifs are normalized by selecting the lexicographically smallest string among all rotations and the reverse complement. This merges equivalent representations (e.g., `TTAGGG`, `GGGTTA`, `CCCTAA` all map to the same canonical form).

---

## Python Library API

```python
from telomerediff import read_fasta_fastq, parse_k, detect_blocks_for_read

for read_id, seq in read_fasta_fastq("reads.fastq.gz"):
    blocks = detect_blocks_for_read(
        read_id, seq,
        k_values=parse_k("6"),
        min_run_bp=150,
    )
    for b in blocks:
        print(f"{b.read}  {b.canonical}  {b.strand}  {b.run_bp}bp  conf={b.confidence}")
```

### Exported API

| Name | Type | Description |
|------|------|-------------|
| `read_fasta_fastq()` | function | Universal sequence reader (FASTQ/FASTA/BAM/CRAM) |
| `parse_k()` | function | Parse k-mer size specification |
| `detect_blocks_for_read()` | function | Detect repeat blocks in a single read |
| `canonical_repeat_unit()` | function | Canonicalize a repeat motif |
| `RepeatBlock` | dataclass | Result object with all block attributes |

---

## Examples

### Stdin piping

```bash
cat reads.fastq | telomerediff -i -
cat sequences.fasta | telomerediff -i - --format fasta
```

### Specific k-mer sizes

```bash
# Human telomere (TTAGGG, k=6)
telomerediff -i reads.fastq.gz -k 6 --refine

# Arabidopsis (TTTAGGG, k=7)
telomerediff -i reads.fastq.gz -k 7 --refine

# Scan a range
telomerediff -i reads.fastq.gz -k 4-15 --refine
```

### BAM input

```bash
telomerediff -i aligned.bam -k 6 --refine
telomerediff -i aligned.cram --reference ref.fa -k 6
```

### Quality filtering

```bash
telomerediff -i reads.fastq.gz --min-quality 10 --refine
```

---

## Project Structure

```
telomerediff/
├── src/telomerediff/
│   ├── __init__.py       # Package exports and version
│   ├── cli.py            # Command-line interface
│   ├── engine.py         # Core detection algorithms
│   ├── io.py             # Multi-format sequence reader
│   ├── repeats.py        # k-mer specification parser
│   ├── report.py         # HTML report generator
│   └── visualize.py      # Matplotlib plotting utilities
├── tests/
│   ├── data/
│   │   └── sample.fastq  # Bundled sample data (5 reads)
│   └── test_basic.py     # 57 tests covering all modules
├── .github/workflows/
│   └── ci.yml            # CI: Python 3.9–3.12 matrix
├── pyproject.toml
├── CITATION.cff
└── LICENSE
```

---

## Testing

```bash
pip install -e ".[dev]"
pytest -v
```

---

## Scientific Rationale

Telomeric regions are composed of tandem repeats and are often poorly represented in reference
genomes due to their repetitive nature. Long-read sequencing (ONT / PacBio) frequently captures
telomeric arrays directly within reads, making alignment-free discovery both feasible and
advantageous — especially for non-model organisms where the telomeric motif may be unknown.

---

## Citation

If you use telomerediff in your research, please cite:

```
Gattu Linga, B.S. (2025). telomerediff: De novo discovery of telomeric repeat blocks
from long-read sequencing. https://github.com/glbala87/telomerediff
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.
