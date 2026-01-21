import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import numpy as np
import pandas as pd
from scipy import stats

from utilities import get_config, load_deseq2_results, load_top_deseq2_results, load_original_genes, CONFIG_DIR

def plot_gene_comparison(original_genes, h19_genes, h38_genes, output_file):
	"""Bar plot comparing fold changes across analyses."""
	genes = original_genes['gene_name'][:15]
	# If we don't have a gene name in modern sets, we can drop it
	genes = genes[genes.isin(h19_genes['gene_name']) & genes.isin(h38_genes['gene_name'])]
	gene_data = {
		'original': original_genes.set_index('gene_name').loc[genes]['log2FoldChange'].values,
		'hg19_modern': h19_genes.set_index('gene_name').loc[genes]['log2FoldChange'].values,
		'hg38_modern': h38_genes.set_index('gene_name').loc[genes]['log2FoldChange'].values,
	}

	x = np.arange(len(genes))
	width = 0.25

	fig, ax = plt.subplots(figsize=(12, len(genes)))

	ax.bar(x - width, gene_data['original'], width, label='Original', color='gray')
	ax.bar(x, gene_data['hg19_modern'], width, label='hg19 Modern', color='steelblue')
	ax.bar(x + width, gene_data['hg38_modern'], width, label='hg38 Modern', color='coral')

	ax.set_ylabel('log2 Fold Change')
	ax.set_xticks(x)
	ax.set_xticklabels(genes, rotation=45)
	ax.legend()
	ax.set_title('Comparison of Key Glucocorticoid-Responsive Genes')

	plt.tight_layout()
	plt.savefig(output_file, dpi=300, bbox_inches='tight')

def plot_deg_overlap(original_degs, hg19_degs, hg38_degs, filename):
	"""Three-way Venn diagram of DEG overlap."""

	plt.figure(figsize=(10, 8))

	venn3([
		set(original_degs['gene_name']),
		set(hg19_degs['gene_name']),
		set(hg38_degs['gene_name']),
	], set_labels=('Original\n(TopHat/Cufflinks+Cuffdiff/hg19)',
				   'Modern\n(STAR/featureCount+DESeq2/hg19)',
				   'Modern\n(STAR/featureCount+DESeq2/hg38)'))

	plt.title('Overlap of Differentially Expressed Genes')
	plt.savefig(filename, dpi=300, bbox_inches='tight')

def plot_lfc_correlation(validated_df, title, filepath):
	"""Scatter plot comparing fold changes."""
	fig, ax = plt.subplots(figsize=(8, 8))

	ax.scatter(validated_df['log2FoldChange_original'], validated_df['log2FoldChange_new'],
			   alpha=0.5, s=10, c='steelblue')

	# Add diagonal line
	lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
			max(ax.get_xlim()[1], ax.get_ylim()[1])]
	ax.plot(lims, lims, 'r--', alpha=0.75, label='y=x')

	# Correlation
	r, p = stats.pearsonr(validated_df['log2FoldChange_original'], validated_df['log2FoldChange_new'])
	print(f"Pearson correlation: {r:.3f}")
	ax.text(0.05, 0.95, f'r = {r:.3f}', transform=ax.transAxes, fontsize=12)

	ax.set_xlabel('log2 Fold Change (Original)')
	ax.set_ylabel('log2 Fold Change (New)')
	ax.set_title(title)

	plt.savefig(filepath, dpi=300, bbox_inches='tight')


def plot_volcano_comparison(results_list: list[pd.DataFrame], titles, output_file):
	"""Side-by-side volcano plots."""

	fig, axes = plt.subplots(1, 3, figsize=(15, 5))

	for ax, results, title in zip(axes, results_list, titles):
		# Significance coloring
		colors = np.where(
			(results['padj'] < 0.05) & (results['log2FoldChange'] > 1), 'red',
			np.where(
				(results['padj'] < 0.05) & (results['log2FoldChange'] < -1), 'blue',
				'gray'
			)
		)
		# if padj is 0 set it to 1e-60 to keep the plots comparable
		results['padj'] = results['padj'].replace(0, 1e-60)

		ax.scatter(results['log2FoldChange'],
				   -np.log10(results['padj']),
				   c=colors, alpha=0.5, s=5)

		ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.5)
		ax.axvline(-1, color='black', linestyle='--', linewidth=0.5)
		ax.axvline(1, color='black', linestyle='--', linewidth=0.5)

		ax.set_xlabel('log2 Fold Change')
		ax.set_ylabel('-log10(padj)')
		ax.set_title(title)

	plt.tight_layout()
	plt.savefig(output_file, dpi=300, bbox_inches='tight')
	plt.close()
	print(f"Saved: {output_file}")
	return fig

def main():
	# Volcano plots
	original_results = load_original_genes()
	hg19_config = get_config(1)
	hg38_config = get_config(2)
	hg19_results = load_top_deseq2_results(hg19_config)
	hg38_results = load_top_deseq2_results(hg38_config)

	plot_volcano_comparison(
		results_list=[original_results, hg19_results[:316], hg38_results[:316]],
		titles=['Original (TopHat/hg19)', 'Modern (STAR/hg19)', 'Modern (STAR/hg38)'],
		output_file=hg19_config.general_results_figures_dir / "volcano_comparison.png")

	# Pearson correlation
	validated_genes_hg19 = pd.read_csv(hg19_config.validated_genes_file)
	validated_genes_hg38 = pd.read_csv(hg38_config.validated_genes_file)

	plot_lfc_correlation(validated_genes_hg38, "log2 Fold Change (Original) vs. New (hg38)", hg38_config.general_results_figures_dir / "lfc_correlation_hg38.png")
	plot_lfc_correlation(validated_genes_hg19, "log2 Fold Change (Original) vs. New (hg19)", hg19_config.general_results_figures_dir / "lfc_correlation_hg19.png")

	# Venn diagram
	plot_deg_overlap(original_results, hg19_results, hg38_results, hg19_config.general_results_figures_dir / "deg_overlap.png")

	# Gene comparison
	plot_gene_comparison(original_results, load_deseq2_results(hg19_config), load_deseq2_results(hg38_config), hg19_config.general_results_figures_dir / "gene_comparison.png")

if __name__ == "__main__":
	main()