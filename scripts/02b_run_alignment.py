import subprocess
import argparse
from pathlib import Path
import pandas as pd
from configuration import get_config, STAR_PARAMS, ROOT, AnalysisConfig

DATA_RAW = ROOT / "data" / "raw"
CONFIG_DIR = ROOT / "config"

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


def collect_alignment_stats(config) -> pd.DataFrame:
	"""Collect alignment statistics from STAR logs."""
	stats = []
	for log_file in config.bam_dir.glob("*.Log.final.out"):
		sample_id = log_file.name.replace(".Log.final.out", "")
		with open(log_file) as f:
			content = f.read()

		# Parse key metrics
		for line in content.split("\n"):
			if "Number of input reads" in line:
				input_reads = int(line.split("|")[1].strip())
			elif "Uniquely mapped reads %" in line:
				unique_pct = float(line.split("|")[1].strip().replace("%", ""))
			elif "% of reads mapped to multiple loci" in line:
				multi_pct = float(line.split("|")[1].strip().replace("%", ""))

		stats.append({
			"sample_id": sample_id,
			"input_reads": input_reads,
			"unique_mapped_pct": unique_pct,
			"multi_mapped_pct": multi_pct,
			"total_mapped_pct": unique_pct + multi_pct
		})

	return pd.DataFrame(stats)


def main():
	parser = argparse.ArgumentParser(description="STAR alignment")
	parser.add_argument("--part", type=int, required=True, choices=[1, 2],
						help="Analysis part, alignment with genome STAR index: 1=hg19, 2=hg38")
	args = parser.parse_args()

	config = get_config(args.part)

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
