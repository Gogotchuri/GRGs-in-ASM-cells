#!/usr/bin/env python3
"""
Build STAR index
Usage:
    python 02_b_build_star_index.py --part 1 // For hg19
    python 02_b_build_star_index.py --part 2 // For hg38
"""

from pathlib import Path
import argparse
from configuration import get_config, STAR_PARAMS, decompress_if_needed
import subprocess

def build_star_index(config) -> None:
	"""Build STAR genome index."""
	star_index = config.reference.star_index

	if (star_index / "SA").exists():
		print(f"STAR index for {config.reference.name} exists, skipping")
		return

	print(f"Building STAR index for {config.reference.name}...")

	star_index.mkdir(parents=True, exist_ok=True)

	# Decompress genome and GTF
	genome_fa = decompress_if_needed(config.reference.genome_file)
	gtf = decompress_if_needed(config.reference.gtf_file)

	cmd = [
		"STAR", "--runMode", "genomeGenerate",
		"--genomeDir", str(star_index),
		"--genomeFastaFiles", str(genome_fa),
		"--sjdbGTFfile", str(gtf),
		"--sjdbOverhang", str(STAR_PARAMS["sjdb_overhang"]),
		"--runThreadN", str(STAR_PARAMS["threads"])
	]
	subprocess.run(cmd, check=True)
	print(f"STAR index built: {star_index}")

def main():
	parser = argparse.ArgumentParser(description="STAR index builder")
	parser.add_argument("--part", type=int, required=True, choices=[1, 2],
						help="For analysis part: 1=hg19, 2=hg38")
	args = parser.parse_args()

	config = get_config(args.part)
	print(f"{config.description}\n")

	# Build index
	build_star_index(config)

if __name__ == "__main__":
	main()
