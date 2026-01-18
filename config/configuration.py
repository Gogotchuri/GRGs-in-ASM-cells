"""
Configuration for fetching and analyzing data

Configures and defines values for (roughly) two-part approach, intending to verify the result of the study and modernize the pipeline.

Part1: hg19 (Original reference genome, used in the study) + New, more modern processing pipeline with STAR and DESeq2
This is done mostly for sanity check of my setup and pipeline, the same data even with faster and improved tooling should result in roughly the same findings.
Part2: hg38 + modern pipeline mentioned above
This is the current reference genome and pipeline, the results should be more accurate. Hopefully we will discover something interesting
"""

from pathlib import Path
from dataclasses import dataclass

ROOT = Path(__file__).parent.parent

@dataclass
class ReferenceConfig:
	"""Reference genome configuration."""
	name: str
	version: str
	genome_url: str
	gtf_url: str
	output_dir: Path

	@property
	def genome_file(self) -> Path:
		return self.output_dir / "genome.fa.gz"

	@property
	def gtf_file(self) -> Path:
		return self.output_dir / "annotation.gtf.gz"

HG19_CONFIG = ReferenceConfig(
	name="hg19",
	version="GRCh37",
	genome_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz",
	gtf_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
	output_dir=ROOT / "data" / "reference" / "hg19",
)

HG38_CONFIG = ReferenceConfig(
	name="hg38",
	version="GRCh38",
	# Gencode v44 (current stable)
	genome_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz",
	gtf_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz",
	output_dir=ROOT / "data" / "reference" / "hg38",
)
