import subprocess
from pathlib import Path
import pandas as pd
from configuration import STAR_PARAMS, AnalysisConfig, CONFIG_DIR, DATA_DIR
from utilities import get_config_based_on_args

DATA_RAW = DATA_DIR / "raw"

def align_sample(srr_id: str, sample_id: str, config: AnalysisConfig) -> Path:
	"""
	Align a single sample with STAR and produce BAM.
	Building the STAR index is assumed to be done beforehand.
	"""
	outdir = config.bam_dir
	outdir.mkdir(parents=True, exist_ok=True)

	bam_out = outdir / f"{sample_id}.Aligned.sortedByCoord.out.bam"
	if bam_out.exists():
		print(f"  {sample_id}: BAM exists, skipping")
		return bam_out

	print(f"  Aligning {sample_id} to {config.reference.name}...")
	fq1 = DATA_RAW / f"{srr_id}_1.fastq.gz"
	fq2 = DATA_RAW / f"{srr_id}_2.fastq.gz"

	cmd = [
		"STAR",
		"--genomeDir", str(config.reference.star_index),
		"--readFilesIn", str(fq1), str(fq2),
		"--readFilesCommand", "zcat", # for reading gzipped FASTQ files
		"--outFileNamePrefix", str(outdir / f"{sample_id}."),
		"--outSAMtype", "BAM", "SortedByCoordinate", # makes BAM work with featureCounts
		"--runThreadN", str(STAR_PARAMS["threads"]),
		"--quantMode", "GeneCounts",
		"--outSAMattributes", "Standard",
		"--outSAMunmapped", "Within",

		"--outFilterType", "BySJout",
		"--outFilterMultimapNmax", "20",
		"--outFilterMismatchNoverReadLmax", "0.04",
	]
	subprocess.run(cmd, check=True)

	return bam_out

def parse_star_log(log_file: Path) -> dict:
	"""
	Parse STAR Log.final.out file.
	Returns columns:
	 ['sample_id', 'input_reads', 'uniquely_mapped', 'uniquely_mapped_pct', 'splice_junctions',
	  'multi_mapped_pct', 'total_mapped_pct', 'junction_pct']
	"""
	sample_id = log_file.name.replace(".Log.final.out", "")
	with open(log_file) as f:
		content = f.read()

	# Line pattern: <Metric Name> | <Value>
	def extract_value(line_n):
		# Read nth line
		line = content.split('\n')[line_n-1]
		# Extract value
		value_str = line.split('|')[1].strip().replace('%', '')
		if '.' in value_str:
			return float(value_str)
		else:
			return int(value_str)

	stats = {
		'sample_id': sample_id,
		'input_reads': extract_value(6),
		'uniquely_mapped': extract_value(9),
		'uniquely_mapped_pct': extract_value(10),
		'splice_junctions': extract_value(12),
		'multi_mapped_pct': extract_value(25),
	}

	stats['total_mapped_pct'] = stats['uniquely_mapped_pct'] + stats['multi_mapped_pct']
	stats['junction_pct'] = (stats['splice_junctions'] / stats['uniquely_mapped']) * 100

	return stats


def collect_alignment_stats(config) -> pd.DataFrame:
	"""Collect alignment statistics from STAR logs."""
	stats = []
	for log_file in config.bam_dir.glob("*.Log.final.out"):
		stats.append(parse_star_log(log_file))

	return pd.DataFrame(stats)


def main():
	config = get_config_based_on_args("STAR alignment")

	# Align samples
	print("\nAligning samples...")
	samples = pd.read_csv(CONFIG_DIR / "samples_table.tsv", sep="\t")

	for _, row in samples.iterrows():
		align_sample(row["srr_id"], row["sample_id"], config)

	# Collect stats
	stats_file = config.bam_dir / "alignment_stats.tsv"
	print("\nCollecting alignment stats into", stats_file)
	stats = collect_alignment_stats(config)
	stats.to_csv(stats_file, sep="\t", index=False)

	print("\n=== Summary ===")
	print(f"Reference: {config.reference.name}")
	print(f"Avg mapping rate: {stats['total_mapped_pct'].mean():.1f}%")
	print(f"Stats saved: {stats_file}")
	print(f"BAMs stored ad: {config.bam_dir}")


if __name__ == "__main__":
	main()
