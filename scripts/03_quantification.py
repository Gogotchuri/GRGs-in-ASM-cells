#!/usr/bin/env python3
#!/usr/bin/env python3
"""
Run featureCounts for gene-level quantification.

Usage:
    python 03_quantification.py --part 1  # hg19 reference
    python 03_quantification.py --part 2  # hg38 reference
"""

import subprocess
import pandas as pd
from pathlib import Path
from configuration import FEATURECOUNTS_PARAMS, CONFIG_DIR
from utilities import decompress_if_needed, get_config_based_on_args

def run_featurecounts(config) -> Path:
    """Run featureCounts on all BAM files."""
    outdir = config.counts_dir
    outdir.mkdir(parents=True, exist_ok=True)

    outfile = outdir / "gene_counts.txt"
    if outfile.exists():
        print("Count matrix exists, skipping featureCounts")
        return outfile

    samples = pd.read_csv(CONFIG_DIR / "samples.tsv", sep="\t")
    bam_files = [str(config.bam_dir / f"{s}.Aligned.sortedByCoord.out.bam")
                 for s in samples["sample_id"]]

    gtf = decompress_if_needed(config.reference.gtf_file)

    print(f"Running featureCounts with {config.reference.name} annotation...")
    cmd = [
              "featureCounts",
              "-p", "--countReadPairs",
              "-T", str(FEATURECOUNTS_PARAMS["threads"]),
              "-a", str(gtf),
              "-o", str(outfile),
              "-t", FEATURECOUNTS_PARAMS["feature_type"],
              "-g", FEATURECOUNTS_PARAMS["attribute"],
              "--extraAttributes", "gene_name",
          ] + bam_files

    subprocess.run(cmd, check=True)
    return outfile

def clean_count_matrix(config) -> Path:
    """Clean featureCounts output"""
    counts_raw = config.counts_dir / "gene_counts.txt"
    counts_clean = config.counts_dir / "counts_matrix.tsv"

    if counts_clean.exists():
        print("Clean count matrix exists")
        return counts_clean

    print("Cleaning count matrix...")
    df = pd.read_csv(counts_raw, sep="\t", comment="#")

    # Rename columns (remove path prefix from BAM names)
    new_cols = []
    for c in df.columns:
        if ".Aligned.sortedByCoord.out.bam" in c:
            new_cols.append(c.split("/")[-1].replace(".Aligned.sortedByCoord.out.bam", ""))
        else:
            new_cols.append(c)
    df.columns = new_cols

    # Keep gene_id, gene_name, and count columns
    samples = pd.read_csv(CONFIG_DIR / "samples_table.tsv", sep="\t")
    keep_cols = ["Geneid", "gene_name"] + list(samples["sample_id"])
    df = df[keep_cols]
    df = df.rename(columns={"Geneid": "gene_id"})

    df.to_csv(counts_clean, sep="\t", index=False)
    print(f"Saved: {counts_clean}")
    return counts_clean

def main():
    config = get_config_based_on_args("Gene quantification with featureCounts")

    # Run featureCounts
    run_featurecounts(config)

    # Clean output
    clean_count_matrix(config)

    print(f"\n=== Quantification complete ===")
    print(f"Counts in: {config.counts_dir}")


if __name__ == "__main__":
    main()
