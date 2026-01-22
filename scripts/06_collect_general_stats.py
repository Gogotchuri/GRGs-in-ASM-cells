#!/usr/bin/env python3
"""
Extract alignment statistics from STAR logs and featureCounts summary.
The idea is to get the same stats as reported in the original paper Results section.
"""

from pathlib import Path
import pandas as pd
from configuration import AnalysisConfig
from utilities import get_config_based_on_args


def parse_featurecounts_summary(summary_file: Path) -> pd.DataFrame:
	"""Parse featureCounts summary file."""
	df = pd.read_csv(summary_file, sep="\t")
	# Clean column names
	new_cols = [df.columns[0]]
	for c in df.columns[1:]:
		new_cols.append(c.split("/")[-1].replace('.Aligned.sortedByCoord.out.bam', ''))
	df.columns = new_cols
	df = df.rename(columns={df.columns[0]: 'sample_id'})
	summary = df.reset_index(drop=True).set_index('sample_id').T.rename(columns={'Assigned':'assigned_reads'})

	summary['assigned_pct'] = 100 * summary['assigned_reads'] / summary.sum(axis=1) # Divide assigned with total reads for sample
	return summary[['assigned_reads', 'assigned_pct']]

def collect_stats(config: AnalysisConfig) -> pd.DataFrame:
	"""
	Collect all statistics slightly differently than in 02b alignment step, we want additional stats here
	And we need averages over everything
	"""
	# stats from STAR alignment built in 02b
	# 	 ['sample_id', 'input_reads', 'uniquely_mapped', 'uniquely_mapped_pct', 'splice_junctions',
	# 	  'multi_mapped_pct', 'total_mapped_pct', 'junction_pct']
	stats_file = config.bam_dir / "alignment_stats.tsv"
	star_stats = pd.read_csv(stats_file, sep="\t").set_index('sample_id')
	# Add featureCounts stats if available
	fc_df = parse_featurecounts_summary(config.counts_dir / "gene_counts.txt.summary")
	return star_stats.join(fc_df)

def main():
	config = get_config_based_on_args("General stat collection")
	df = collect_stats(config)
	st = generate_summary_table(df)
	st.to_csv(config.results_tables_dir / "general_stats_summary.tsv", sep="\t", index=False)

def generate_summary_table(df: pd.DataFrame) -> pd.DataFrame:
	"""Generate summary table matching paper format."""
	summary = {
		'avg_reads': df['input_reads'].mean()*2,
		'min_reads': df['input_reads'].min()*2,
		'max_reads': df['input_reads'].max()*2,

		'avg_align_pct': df['total_mapped_pct'].mean(),
		'min_align_pct': df['total_mapped_pct'].min(),
		'max_align_pct': df['total_mapped_pct'].max(),

		'avg_junction_pct': df['junction_pct'].mean(),
		'min_junction_pct': df['junction_pct'].min(),
		'max_junction_pct': df['junction_pct'].max(),

		'avg_assigned_pct': df['assigned_pct'].mean(),
		'min_assigned_pct': df['assigned_pct'].min(),
		'max_assigned_pct': df['assigned_pct'].max(),
	}

	# Round each value to 2 decimals
	for k in summary.keys():
		summary[k] = round(summary[k], 2)

	return pd.DataFrame([summary], columns=summary.keys())

if __name__ == "__main__":
	main()