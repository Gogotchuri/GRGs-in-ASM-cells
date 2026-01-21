import argparse
import pandas as pd
from configuration import ROOT, get_config

CONFIG_DIR = ROOT / "config"

def load_original_genes() -> pd.DataFrame:
	return pd.read_csv(CONFIG_DIR / "modified_original_genes.tsv", sep="\t")

def load_top_deseq2_results(config) -> pd.DataFrame:
	return pd.read_csv(config.top_deseq2_results_file)

def load_deseq2_results(config) -> pd.DataFrame:
	return pd.read_csv(config.deseq2_results_file)

def validate_original_genes(original_genes: pd.DataFrame, new_genes: pd.DataFrame) -> pd.DataFrame:
	"""
	Validate that new genes are in the original gene list.
	Creates a new dataframe with validated genes, genes that have a match in new genes from original list,
	Rcords gene_name, padj_original, padj_new, log2FoldChange_original, log2FoldChange_new, direction_match (if the FC has the same sign)
	At the end gathers stats: validated_perc, correct_direction_perc
	"""
	# Inner join genes on gene_name
	merged = original_genes.merge(
		new_genes[["gene_name", "padj", "log2FoldChange"]],
		on="gene_name",
		how="inner",
		suffixes=("_original", "_new")
	)

	# Check if fold change direction matches (same sign)
	merged["direction_match"] = (
			(merged["log2FoldChange_original"] * merged["log2FoldChange_new"]) > 0
	)

	# Keep only the required columns
	merged = merged[[
		"gene_name",
		"padj_original",
		"padj_new",
		"log2FoldChange_original",
		"log2FoldChange_new",
		"direction_match"
	]]

	# Calculate statistics
	total_original = len(original_genes)
	validated_count = len(merged)
	correct_direction_count = merged["direction_match"].sum()

	validated_perc = (validated_count / total_original) * 100
	correct_direction_perc = (correct_direction_count / validated_count) * 100 if validated_count > 0 else 0

	print(f"Original genes: {total_original}")
	print(f"Validated genes (found in new results): {validated_count} ({validated_perc:.1f}%)")
	print(f"Correct direction (same FC sign): {correct_direction_count} ({correct_direction_perc:.1f}%)")

	return merged.sort_values("padj_new")

def main():
	parser = argparse.ArgumentParser(description="DESeq2 analysis")
	parser.add_argument("--part", type=int, required=True, choices=[1, 2],
						help="Analysis part: 1=hg19, 2=hg38")
	args = parser.parse_args()

	config = get_config(args.part)
	print(f"=== {config.description} ===\n")
	config.results_tables_dir.mkdir(parents=True, exist_ok=True) # Just in case

	original_genes = load_original_genes()
	new_top_genes = load_top_deseq2_results(config)

	print(f"\n=== Validation Statistics (Top Genes from DESeq2) ===")
	validated_df = validate_original_genes(original_genes, new_top_genes)
	validated_df.to_csv(config.top_validated_genes_file, index=False)
	print(f"Validated genes saved to {config.top_validated_genes_file}")

	print(f"\n=== Validation Statistics (From all significant genes) ===")
	new_all_genes = load_deseq2_results(config)
	validated_df = validate_original_genes(original_genes, new_all_genes)
	validated_df.to_csv(config.validated_genes_file, index=False)
	print(f"Validated genes saved to {config.validated_genes_file}")

if __name__ == "__main__":
	main()
