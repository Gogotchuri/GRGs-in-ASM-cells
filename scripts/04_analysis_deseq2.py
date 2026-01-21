import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from configuration import CONFIG_DIR
from utilities import get_config_based_on_args

MIN_GENE_COUNTS = 10

def load_counts(config) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
	"""Load count matrix and sample metadata."""
	counts_file = config.counts_dir / "counts_matrix.tsv"
	counts = pd.read_csv(counts_file, sep="\t", index_col=0)

	# Store gene names separately
	gene_names = counts["gene_name"].copy() if "gene_name" in counts.columns else None
	if "gene_name" in counts.columns:
		counts = counts.drop(columns=["gene_name"])

	samples = pd.read_csv(CONFIG_DIR / "samples_table.tsv", sep="\t", index_col="sample_id")

	# Ensure counts columns match sample order
	counts = counts[samples.index.tolist()]

	return counts, samples, gene_names

def filter_genes(counts: pd.DataFrame) -> pd.DataFrame:
	"""Filter genes with low counts."""
	counts_filtered = counts[counts.sum(axis=1) >= MIN_GENE_COUNTS]
	return counts_filtered

def run_deseq2(counts: pd.DataFrame, samples: pd.DataFrame) -> tuple[pd.DataFrame, DeseqDataSet]:
	"""Run PyDESeq2 differential expression analysis."""
	print("Running PyDESeq2...")

	# PyDESeq2 expects samples as rows, genes as columns, so we transpose
	# We also filter genes with too low counts to speed up analysis
	counts_T = filter_genes(counts).T

	# Create DESeq2 dataset
	dds = DeseqDataSet(
		counts=counts_T,
		metadata=samples,
		# I have opted to not include cell line in the design
		design="~ condition",
		refit_cooks=True,
		n_cpus=8
	)

	# default deseq2 pipeline https://pydeseq2.readthedocs.io/en/stable/api/docstrings/pydeseq2.dds.DeseqDataSet.html#pydeseq2.dds.DeseqDataSet.deseq2
	dds.deseq2()

	# Statistical testing
	stat_res = DeseqStats(dds, contrast=["condition", "untreated", "dex"], n_cpus=8, alpha=0.05)
	stat_res.summary()

	return stat_res.results_df, dds

def annotate_results(results: pd.DataFrame, gene_names: pd.Series) -> pd.DataFrame:
	"""Add gene names and sort by adjusted p-value."""
	results = results.copy()

	if gene_names is not None:
		results["gene_name"] = results.index.map(lambda x: gene_names.get(x, x))
	else:
		results["gene_name"] = results.index

	return results.sort_values("padj")

def summarize_results(results: pd.DataFrame) -> dict:
	"""Generate summary statistics."""
	total = len(results)
	sig_05 = (results["padj"] < 0.05).sum()
	sig_strict = ((results["padj"] < 0.05) & (abs(results["log2FoldChange"]) > 1)).sum()
	up = ((results["padj"] < 0.05) & (results["log2FoldChange"] > 1)).sum()
	down = ((results["padj"] < 0.05) & (results["log2FoldChange"] < -1)).sum()

	print("\n=== Analysis Summary ===")
	print(f"Total genes tested:        {total:,}")
	print(f"Significant (padj<0.05):   {sig_05:,}  [Original: 316]")
	print(f"Top Genes (padj<0.05, |LFC| > 1): {sig_strict:,}")
	print(f"   Upregulated:           {up:,}")
	print(f"   Downregulated:         {down:,}")


def main():
	config = get_config_based_on_args("DESeq2 analysis")

	# Create output directory
	config.results_tables_dir.mkdir(parents=True, exist_ok=True)

	# Load data
	print("Loading count data...")
	counts, samples, gene_names = load_counts(config)
	print(f"Loaded {counts.shape[0]} genes, {counts.shape[1]} samples")

	# Run DESeq2
	results, dds = run_deseq2(counts, samples)
	results = annotate_results(results, gene_names)

	# Summary
	summarize_results(results)

	results.to_csv(config.deseq2_results_file)
	results[(results["padj"] < 0.05) & (abs(results["log2FoldChange"]) > 1)].to_csv(config.top_deseq2_results_file)

if __name__ == "__main__":
	main()
