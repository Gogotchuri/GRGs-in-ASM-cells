#!/usr/bin/env python3
"""
Download FASTQ files from SRA and reference genomes.
"""

import subprocess
import pandas as pd

from configuration import HG19_CONFIG, HG38_CONFIG, ROOT

CONFIG = ROOT / "config"


def download_sra(srr_id: str) -> None:
    """Download and convert SRA to FASTQ using fasterq-dump (sra-tools)."""
    outdir = ROOT / "data" / "raw"
    outdir.mkdir(parents=True, exist_ok=True)

    fq1 = outdir / f"{srr_id}_1.fastq.gz"
    if fq1.exists():
        print(f"  {srr_id}: exists, skipping")
        return

    print(f"  {srr_id}: downloading...")
    subprocess.run([
        "fasterq-dump", srr_id,
        "--outdir", str(outdir),
        "--threads", "4",
        "--split-files",
        "--progress"
    ], check=True)

    for fq in outdir.glob(f"{srr_id}*.fastq"):
        print(f"  Compressing {fq.name}...")
        subprocess.run(["pigz", "-p", "4", str(fq)], check=True)


def download_reference(config) -> None:
    """Download reference genome and annotation."""
    config.output_dir.mkdir(parents=True, exist_ok=True)

    for url, outfile in [(config.genome_url, config.genome_file),
                         (config.gtf_url, config.gtf_file)]:
        # Checking for safety for retrial safety
        if outfile.exists():
            print(f"  {outfile.name}: exists, skipping")
            continue
        print(f"  Downloading {outfile.name}...")
        subprocess.run(["wget", "-q", "--show-progress", "-O", str(outfile), url], check=True)


def main():
    print("=== Downloading FASTQ files (based on samples_table.tsv) ===")
    samples = pd.read_csv(CONFIG / "samples_table.tsv", sep="\t")
    for _, row in samples.iterrows():
        download_sra(row["srr_id"])

    print("\n=== Downloading hg19 reference ===")
    download_reference(HG19_CONFIG)

    print("\n=== Downloading hg38 reference ===")
    download_reference(HG38_CONFIG)

    print("\nDownload Complete!")


if __name__ == "__main__":
    main()
