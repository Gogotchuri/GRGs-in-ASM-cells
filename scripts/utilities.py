
import pandas as pd
from pathlib import Path
import subprocess
import argparse
from configuration import AnalysisConfig, get_config, CONFIG_DIR


def load_original_genes() -> pd.DataFrame:
	return pd.read_csv(CONFIG_DIR / "modified_original_genes.tsv", sep="\t")

def load_top_deseq2_results(config) -> pd.DataFrame:
	return pd.read_csv(config.top_deseq2_results_file)

def load_deseq2_results(config) -> pd.DataFrame:
	return pd.read_csv(config.deseq2_results_file)

def get_config_based_on_args(description) -> AnalysisConfig:
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument("--part", type=int, required=True, choices=[1, 2],
						help="Analysis part: 1=hg19, 2=hg38")
	args = parser.parse_args()

	config = get_config(args.part)
	print(f"=== {config.description} ===\n")
	return config

def decompress_if_needed(gz_file: Path) -> Path:
	"""Decompress .gz file if uncompressed version doesn't exist yet."""
	uncompressed = gz_file.with_suffix("")
	if not uncompressed.exists() and gz_file.exists():
		print(f"Decompressing {gz_file.name}...")
		subprocess.run(["gunzip", "-k", str(gz_file)], check=True)
	return uncompressed
