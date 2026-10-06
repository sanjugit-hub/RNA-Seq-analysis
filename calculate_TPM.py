# ============================================================
# Purpose: Calculate gene-level TPM from FeatureCounts raw
#          counts using union exon lengths from the GTF file.
#          Produces a gene × sample TPM matrix.
# Author: SBM
# ============================================================

import pandas as pd
import re
import os
import sys


# ------------------------------------------------------------
# INPUT / OUTPUT FILES
# ------------------------------------------------------------

GTF = "Reference_genome/Homo_sapiens.GRCh38.116.gtf"

FEATURECOUNTS = "counts/gene_counts/testis_MIsoTEx_featurecounts_gene_counts.txt"

OUTPUT = "counts/gene_counts/testis_MIsoTEx/gene_TPM_testis_MIsoTEx.txt"


# ------------------------------------------------------------
# STEP 1: Calculate union exon length for each gene
# ------------------------------------------------------------

print("===== STEP 1: Calculating gene exon lengths =====")

exons = {}

with open(GTF, "r") as f:

    for line in f:

        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) < 9:
            continue

        feature = fields[2]

        # Only use exon features
        if feature != "exon":
            continue

        chrom = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        attributes = fields[8]

        # Extract gene_id
        match = re.search(r'gene_id "([^"]+)"', attributes)

        if not match:
            continue

        gene_id = match.group(1)

        if gene_id not in exons:
            exons[gene_id] = []

        exons[gene_id].append(
            (chrom, start, end)
        )


# ------------------------------------------------------------
# Merge overlapping exons and calculate union length
# ------------------------------------------------------------

gene_lengths = {}

for gene_id, regions in exons.items():

    # Sort by chromosome and start
    regions.sort(key=lambda x: (x[0], x[1]))

    merged = []

    for chrom, start, end in regions:

        if not merged:

            merged.append(
                [chrom, start, end]
            )

        else:

            last_chrom = merged[-1][0]
            last_start = merged[-1][1]
            last_end = merged[-1][2]

            # Merge overlapping/adjacent exons
            if chrom == last_chrom and start <= last_end + 1:

                merged[-1][2] = max(
                    last_end,
                    end
                )

            else:

                merged.append(
                    [chrom, start, end]
                )

    # Calculate total union exon length
    length = sum(
        end - start + 1
        for chrom, start, end in merged
    )

    gene_lengths[gene_id] = length


print(
    f"Genes with exon lengths: {len(gene_lengths)}"
)


# ------------------------------------------------------------
# STEP 2: Read FeatureCounts output
# ------------------------------------------------------------

print("\n===== STEP 2: Reading FeatureCounts counts =====")

counts = pd.read_csv(
    FEATURECOUNTS,
    sep="\t",
    comment="#"
)

# Remove leading/trailing spaces
counts.columns = counts.columns.str.strip()

print(
    f"Genes in FeatureCounts: {len(counts)}"
)


# ------------------------------------------------------------
# Identify annotation and sample columns
# ------------------------------------------------------------

annotation_columns = [
    "Geneid",
    "Chr",
    "Start",
    "End",
    "Strand",
    "Length",
    "gene_name"
]

sample_columns = [
    col
    for col in counts.columns
    if col not in annotation_columns
]


print(
    f"Number of samples: {len(sample_columns)}"
)

print("Samples:")

for sample in sample_columns:
    print("  ", sample)


# ------------------------------------------------------------
# Convert sample count columns to numeric
# ------------------------------------------------------------

for sample in sample_columns:

    counts[sample] = pd.to_numeric(
        counts[sample],
        errors="coerce"
    ).fillna(0)


# ------------------------------------------------------------
# STEP 3: Match FeatureCounts genes with GTF
# ------------------------------------------------------------

print("\n===== STEP 3: Matching genes to GTF =====")

counts["gene_length"] = counts["Geneid"].map(
    gene_lengths
)

matched = counts["gene_length"].notna().sum()

unmatched = counts["gene_length"].isna().sum()

print(
    f"Genes matched with GTF exon lengths: {matched}"
)

print(
    f"Genes without GTF exon length: {unmatched}"
)


# ------------------------------------------------------------
# Keep genes with valid exon lengths
# ------------------------------------------------------------

counts = counts[
    counts["gene_length"].notna()
].copy()

counts = counts[
    counts["gene_length"] > 0
].copy()


# ------------------------------------------------------------
# STEP 4: Calculate TPM
# ------------------------------------------------------------

print("\n===== STEP 4: Calculating TPM =====")

# Convert gene length to kilobases
gene_length_kb = (
    counts["gene_length"] / 1000
)


# Keep Geneid and gene_name in output
tpm = counts[
    ["Geneid", "gene_name"]
].copy()


for sample in sample_columns:

    print(
        f"Calculating TPM: {sample}"
    )

    # Raw gene counts
    raw_counts = counts[sample]

    # --------------------------------------------------------
    # RPK = raw count / gene length in kb
    # --------------------------------------------------------

    rpk = (
        raw_counts /
        gene_length_kb
    )

    # --------------------------------------------------------
    # TPM = RPK / sum(RPK) × 1,000,000
    # --------------------------------------------------------

    total_rpk = rpk.sum()

    if total_rpk > 0:

        tpm[sample] = (
            rpk /
            total_rpk
        ) * 1_000_000

    else:

        tpm[sample] = 0


# ------------------------------------------------------------
# STEP 5: Save TPM matrix
# ------------------------------------------------------------

print("\n===== STEP 5: Writing TPM output =====")

output_dir = os.path.dirname(OUTPUT)

if output_dir:
    os.makedirs(
        output_dir,
        exist_ok=True
    )


tpm.to_csv(
    OUTPUT,
    sep="\t",
    index=False
)


# ------------------------------------------------------------
# STEP 6: Basic validation
# ------------------------------------------------------------

print("\n===== STEP 6: TPM validation =====")

for sample in sample_columns:

    total_tpm = tpm[sample].sum()

    print(
        f"{sample}: total TPM = {total_tpm:.2f}"
    )


print("\n==============================================")
print("TPM calculation completed successfully.")
print(f"Output file: {OUTPUT}")
print(f"Genes: {len(tpm)}")
print(f"Samples: {len(sample_columns)}")
print("==============================================")

