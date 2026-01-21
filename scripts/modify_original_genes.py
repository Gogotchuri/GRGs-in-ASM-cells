"""
This script is used to modify the original genes, adding column with log2 change instead of natural logarithm
Since the deseq2 uses log2 this makes the comparison between our results and the ones reported from original paper easier
"""


import pandas as pd
import numpy as np
from configuration import CONFIG_DIR

def add_log2_change(df: pd.DataFrame) -> pd.DataFrame:
	""" Convert Ln[Fold Change] to log2"""
	return df.assign(log2FoldChange=df["Ln[Fold Change]"] / np.log(2))

def main():
	df = pd.read_csv(CONFIG_DIR / "original_genes.tsv", sep="\t")
	df = add_log2_change(df)
	df.to_csv(CONFIG_DIR / "modified_original_genes.tsv", sep="\t", index=False)

if __name__ == "__main__":
	main()
