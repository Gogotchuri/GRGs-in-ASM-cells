import pandas as pd
from configuration import get_config

def load_top_genes(config) -> pd.DataFrame:
   return pd.read_csv(config.top_deseq2_results_file)

def load_validated_genes(config) -> pd.DataFrame:
   return pd.read_csv(config.top_validated_genes_file)

def load_not_replicated_genes(config) -> pd.DataFrame:
   return pd.read_csv(config.results_tables_dir / "not_replicated_genes.csv")

def left_outer_join(left_df, right_df):
    return left_df[~left_df["gene_name"].isin(right_df["gene_name"])]

def adjust_validation_columns(df: pd.DataFrame) -> pd.DataFrame:
    return df[["gene_name", "padj_original", "padj_new", "log2FoldChange_original", "log2FoldChange_new"]]

def adjust_results_columns(df: pd.DataFrame) -> pd.DataFrame:
    return df[["gene_name", "padj", "log2FoldChange"]]

def main():
    hg19_config = get_config(1)
    hg38_config = get_config(2)
    # Validated
    hg19_validated_genes = load_validated_genes(hg19_config)
    hg38_validated_genes = load_validated_genes(hg38_config)
    # Top
    hg19_top_genes = load_top_genes(hg19_config)
    hg38_top_genes = load_top_genes(hg38_config)

    ## Evaluate the symmetric difference between the results
    # Validated genes in hg19 that has been discarded after analysis with the new genome
    validated_in_hg19_discarded_in_hg38 = adjust_validation_columns(left_outer_join(hg19_validated_genes, hg38_validated_genes))
    # Additional genes that has been validated due to new genome (We haven't validated them with hg19)
    new_validations_with_hg38 = adjust_validation_columns(left_outer_join(hg38_validated_genes, hg19_validated_genes))
    # Top genes discarded after analysis with the new genome
    top_genes_from_h19_discarded_with_hg38 = adjust_results_columns(left_outer_join(hg19_top_genes, hg38_top_genes))
    # Top genes discovered after analysis with the new genome
    new_genes_with_hg38 = adjust_results_columns(left_outer_join(hg38_top_genes, hg19_top_genes))
    # Completely new significant genes discovered after analysis with the new genome
    completely_new_genes_with_hg38 = adjust_results_columns(left_outer_join(new_genes_with_hg38, hg38_validated_genes))

    print("Validated genes in hg19 that has been discarded after analysis with hg38:\n", validated_in_hg19_discarded_in_hg38)
    print("\nAdditional genes validated in hg38 that has been discarded initially with hg19:\n", new_validations_with_hg38)
    print("\nTop genes discarded after analysis with hg38:\n", top_genes_from_h19_discarded_with_hg38[:10])
    print("\nTop new genes discovered after analysis with hg38:\n", new_genes_with_hg38[:10])
    print("\nCompletely new significant genes discovered after analysis with hg38:\n", completely_new_genes_with_hg38[:10])

    # Create results directory if it doesn't exist
    results_dir = hg19_config.general_results_tables_dir
    results_dir.mkdir(parents=True, exist_ok=True)

    # Write csv files
    validated_in_hg19_discarded_in_hg38.to_csv(results_dir / "validated_in_hg19_discarded_in_hg38.csv", index=False)
    new_validations_with_hg38.to_csv(results_dir / "new_validations_with_hg38.csv", index=False)
    top_genes_from_h19_discarded_with_hg38.to_csv(results_dir / "top_genes_from_h19_discarded_with_hg38.csv", index=False)
    new_genes_with_hg38.to_csv(results_dir / "new_genes_with_hg38_from_h19.csv", index=False)
    completely_new_genes_with_hg38.to_csv(results_dir / "completely_new_genes_with_hg38.csv", index=False)

    print("Results written to:", results_dir)









if __name__ == "__main__":
    main()