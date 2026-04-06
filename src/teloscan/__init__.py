__all__ = [
    "__version__",
    "RepeatBlock",
    "canonical_repeat_unit",
    "detect_blocks_for_read",
    "read_fasta_fastq",
    "parse_k",
]
__version__ = "0.2.1"

from teloscan.engine import RepeatBlock, canonical_repeat_unit, detect_blocks_for_read
from teloscan.io import read_fasta_fastq
from teloscan.repeats import parse_k
