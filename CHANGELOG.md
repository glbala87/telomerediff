# Changelog

All notable changes to TeloScan will be documented in this file.

## [0.2.1] - 2025

### Added
- Strand-aware detection with proper forward/reverse complement handling
- Confidence scoring combining copy count, sequence entropy, and detection mode
- Multi-threaded processing with configurable thread count and chunk size
- GFF3 and VCF output formats
- Self-contained HTML report with embedded plots (matplotlib)
- PNG plot generation (block length distribution, motif abundance, per-read heatmap, strand distribution)
- Quality filtering by average base quality
- BAM/CRAM input support via pysam
- Gzip-compressed FASTA/FASTQ support
- Stdin input with automatic buffering for two-pass analysis
- Format auto-detection (extension + magic byte peek)
- `--format` flag to override auto-detection
- Python library API (`detect_blocks_for_read`, `canonical_repeat_unit`, `read_fasta_fastq`, `parse_k`, `RepeatBlock`)
- Comprehensive test suite (57 tests) covering all modules
- Security: XSS prevention in HTML reports, field sanitization in GFF3/VCF/BED output
- CI: GitHub Actions with Python 3.9–3.12 matrix
- Input validation for all CLI arguments with early error reporting
- k-mer size validation warnings for unreasonable values
- Documentation of stdin memory implications and other limitations

### Changed
- Two-pass algorithm: Pass 1 discovers motifs, Pass 2 detects blocks (memory-efficient streaming)
- k-mer specification now supports single values, ranges, comma-separated, and mixed formats
- Incremental output writing during Pass 2 for memory efficiency

## [0.1.0] - 2025

### Added
- Initial release
- Alignment-free de novo telomeric repeat discovery from long-read sequencing
- Basic k-mer scanning with configurable k-mer sizes
- Error-tolerant fuzzy matching (Hamming distance)
- Canonical motif normalization (all rotations + reverse complement)
- TSV and BED output formats
- FASTA/FASTQ input support
