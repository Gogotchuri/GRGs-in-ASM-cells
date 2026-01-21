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
import subprocess

ROOT = Path(__file__).parent.parent
CONFIG_DIR = ROOT / "config"
DATA_DIR = ROOT / "data"

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

	@property
	def star_index(self) -> Path:
		return self.output_dir / "star_index"


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

@dataclass
class AnalysisConfig:
	"""Analysis run configuration. Separates files for each analysis part. Under designamted genome directories."""
	part: int
	reference: ReferenceConfig
	output_suffix: str
	description: str

	@property
	def bam_dir(self) -> Path:
		return ROOT / "data" / "processed" / f"bam_{self.reference.name}"

	@property
	def counts_dir(self) -> Path:
		return ROOT / "data" / "processed" / f"counts_{self.reference.name}"

	@property
	def results_dir(self) -> Path:
		return ROOT / "results" / self.reference.name

	@property
	def general_results_tables_dir(self) -> Path:
		return ROOT / "results" / "general" / "tables"

	@property
	def results_tables_dir(self) -> Path:
		return self.results_dir / "tables"

	@property
	def results_figures_dir(self) -> Path:
		return self.results_dir / "figures"

	@property
	def deseq2_results_file(self) -> Path:
		return self.results_tables_dir / "deseq2_results.csv"

	@property
	def top_deseq2_results_file(self) -> Path:
		return self.results_tables_dir / "top_results_from_deseq2.csv"

	@property
	def validated_genes_file(self) -> Path:
		return self.results_tables_dir / "validated_genes.csv"

	@property
	def top_validated_genes_file(self) -> Path:
		return self.results_tables_dir / "top_validated_genes.csv"

# Analysis configurations
PART1_CONFIG = AnalysisConfig(
	part=1,
	reference=HG19_CONFIG,
	output_suffix="hg19",
	description="Part 1: hg19",
)

PART2_CONFIG = AnalysisConfig(
	part=2,
	reference=HG38_CONFIG,
	output_suffix="hg38",
	description="Part 2: hg38",
)


def get_config(part: int) -> AnalysisConfig:
	"""Get configuration for the specified analysis part."""
	if part == 1:
		return PART1_CONFIG
	elif part == 2:
		return PART2_CONFIG
	else:
		raise ValueError(f"Invalid part: {part}. Must be 1 or 2.")


# Pipeline parameters
STAR_PARAMS = {
	"threads": 16,
	"sjdb_overhang": 62,  # read_length - 1
	"out_sam_type": "BAM SortedByCoordinate", # For featureCount
}

FEATURECOUNTS_PARAMS = {
	"threads": 16,
	"feature_type": "exon",
	"attribute": "gene_id",
}
