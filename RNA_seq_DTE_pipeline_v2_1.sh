#!/bin/bash
set -euo pipefail

# ==========================
# Input arguments
# ==========================
JOB_ID="$1"
STUDY_ID="$2"

# ==========================
# Project paths
# ==========================
PROJECT_ROOT="${PROJECT_ROOT:-/home/group1/cgi-bin/rna_seq}"

if [ -z "${ROOT_FOLDER:-}" ]; then
    echo "ERROR: ROOT_FOLDER not set"
    exit 1
fi

BASE="$PROJECT_ROOT/$ROOT_FOLDER"
RAW_FASTQ_ROOT="$BASE/raw_fastq"

MASTER_XLS="$PROJECT_ROOT/SRA/rnaseq_master_1.xlsx"
HISAT_BASE="$PROJECT_ROOT/hisat_deseq_stringtei"
METADATA_FILE="$HISAT_BASE/metadata/${JOB_ID}_metadata.csv"
# ==========================
# Tool paths / settings
# ==========================
THREADS=8

REF_GTF="$HISAT_BASE/Reference_genome/Homo_sapiens.GRCh38.113.gtf"
HISAT2_INDEX="$HISAT_BASE/Reference_genome/hg38_index"
THREADS=8
TRIMMOMATIC_JAR="$HISAT_BASE/trimmomatic-0.39.jar"
ADAPTERS="$HISAT_BASE/TruSeq3-PE.fa"
TRIMMED_DIR="$HISAT_BASE/trimmed_reads"
mkdir -p "$TRIMMED_DIR"
BAM_DIR="$HISAT_BASE/bam"

# ==========================
# Output folders
# ==========================
mkdir -p "$HISAT_BASE/logs" "$HISAT_BASE/bam"
mkdir -p "$HISAT_BASE/counts/gene_counts" "$HISAT_BASE/counts/transcript_counts"
mkdir -p "$HISAT_BASE/gtf_individual" "$HISAT_BASE/gtf_merged"
mkdir -p "$HISAT_BASE/ballgown" "$HISAT_BASE/deseq2"
mkdir -p "$HISAT_BASE/metadata"
# ==========================
# Log files
# ==========================
HISAT_LOG="$HISAT_BASE/logs/hisat_log.txt"
STRINGTIE_LOG="$HISAT_BASE/logs/stringtie_log.txt"
FEATURECOUNT_LOG="$HISAT_BASE/logs/featurecount_log.txt"
DESEQ_LOG="$HISAT_BASE/logs/deseq_log.txt"
BALLGOWN_LOG="$HISAT_BASE/logs/ballgown_log.txt"
TRIM_LOG="$HISAT_BASE/logs/trim_log.txt"

touch "$TRIM_LOG"

touch "$HISAT_LOG" "$STRINGTIE_LOG" "$FEATURECOUNT_LOG" "$DESEQ_LOG" "$BALLGOWN_LOG"

# ==========================
# Check input files
# ==========================
[ -f "$MASTER_XLS" ] || { echo "ERROR: Master Excel not found: $MASTER_XLS"; exit 1; }
[ -f "$METADATA_FILE" ] || { echo "ERROR: Metadata file not found: $METADATA_FILE"; exit 1; }
[ -f "$REF_GTF" ] || { echo "ERROR: GTF not found: $REF_GTF"; exit 1; }
[ -f "${HISAT2_INDEX}.1.ht2" ] || { echo "ERROR: HISAT2 index not found: $HISAT2_INDEX"; exit 1; }

# ==========================
# Get sample IDs from master Excel
# ==========================
SAMPLE_LIST=$(python3 <<EOF
import pandas as pd
import sys

master = "$MASTER_XLS"
study_id = "$STUDY_ID"

df = pd.read_excel(master)
df["Study ID"] = df["Study ID"].ffill()

df_sub = df[df["Study ID"].astype(str).str.strip() == study_id]

if df_sub.empty:
    sys.exit("ERROR: No samples found for STUDY_ID: " + study_id)

samples = df_sub["Sample ID"].dropna().astype(str).str.strip().tolist()

for s in samples:
    print(s)
EOF
)

if [ -z "$SAMPLE_LIST" ]; then
    echo "ERROR: No samples found from master file"
    exit 1
fi

echo "Samples selected from master file:"
echo "$SAMPLE_LIST"

echo "===== STEP 0: TRIMMOMATIC + HISAT2 ====="

while read -r SAMPLE; do

    RAW_R1="$RAW_FASTQ_ROOT/$SAMPLE/${SAMPLE}_1.fastq.gz"
    RAW_R2="$RAW_FASTQ_ROOT/$SAMPLE/${SAMPLE}_2.fastq.gz"

    TRIM_R1="$TRIMMED_DIR/${SAMPLE}_R1.trimmed.fastq.gz"
    TRIM_R2="$TRIMMED_DIR/${SAMPLE}_R2.trimmed.fastq.gz"

    UNP_R1="$TRIMMED_DIR/${SAMPLE}_R1.unpaired.fastq.gz"
    UNP_R2="$TRIMMED_DIR/${SAMPLE}_R2.unpaired.fastq.gz"

    BAM_OUT="$BAM_DIR/${SAMPLE}.bam"

    if ! grep -Fxq "$SAMPLE" "$TRIM_LOG"; then
        java -jar "$TRIMMOMATIC_JAR" PE \
            -threads "$THREADS" \
            "$RAW_R1" "$RAW_R2" \
            "$TRIM_R1" "$UNP_R1" \
            "$TRIM_R2" "$UNP_R2" \
            ILLUMINACLIP:"$ADAPTERS":2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

        echo "$SAMPLE" >> "$TRIM_LOG"
    fi

    if ! grep -Fxq "$SAMPLE" "$HISAT_LOG"; then
        hisat2 -p "$THREADS" \
            -x "$HISAT2_INDEX" \
            -1 "$TRIM_R1" \
            -2 "$TRIM_R2" | \
            samtools sort -@ "$THREADS" -o "$BAM_OUT"

        samtools index "$BAM_OUT"
        echo "$SAMPLE" >> "$HISAT_LOG"
    fi

done <<< "$SAMPLE_LIST"

# ==========================
# Build BAM list only for current study samples
# ==========================
BAM_FILES=()

while read -r SAMPLE; do
    [[ -z "$SAMPLE" || "$SAMPLE" == \#* ]] && continue
    BAM_FILES+=("$BAM_DIR/${SAMPLE}.bam")
done <<< "$SAMPLE_LIST"

echo "===== ALIGNMENT COMPLETED ====="

# ==========================
# STEP 1: FeatureCounts
# ==========================
echo "===== STEP 1: FeatureCounts gene quantification ====="

GENE_DIR="$HISAT_BASE/counts/gene_counts"
TRANSCRIPT_DIR="$HISAT_BASE/counts/transcript_counts"

FEATURECOUNTS_OUT="$GENE_DIR/${STUDY_ID}_featurecounts_gene_counts.txt"

if grep -Fxq "$STUDY_ID" "$FEATURECOUNT_LOG"; then
    echo "Already processed: FeatureCounts for $STUDY_ID"
else
    featureCounts \
        -T "$THREADS" \
        -p \
        -t exon \
        -g gene_id \
        --extraAttributes gene_name \
        -a "$REF_GTF" \
        -o "$FEATURECOUNTS_OUT" \
        "${BAM_FILES[@]}"

    echo "$STUDY_ID" >> "$FEATURECOUNT_LOG"
fi

# ==========================
# STEP 2: Convert featureCounts output to DESeq2 matrix
# ==========================
echo "===== STEP 2: Preparing DESeq2 gene matrix ====="

GENE_MATRIX="$GENE_DIR/${STUDY_ID}_gene_count_matrix_featurecounts.csv"
TEMP_COUNTS="$HISAT_BASE/counts/temp_counts_${JOB_ID}.tsv"

if [ -f "$GENE_MATRIX" ]; then
    echo "Gene matrix already exists: $GENE_MATRIX"
else
    grep -v '^#' "$FEATURECOUNTS_OUT" | cut -f1,7- > "$TEMP_COUNTS"

    python3 <<EOF
import pandas as pd

df = pd.read_csv("$TEMP_COUNTS", sep="\t")
df = df.set_index(df.columns[0])
df.columns = [c.split("/")[-1].replace(".bam","") for c in df.columns]
df.to_csv("$GENE_MATRIX")

print("Matrix shape:", df.shape)
EOF

    rm "$TEMP_COUNTS"
fi

echo "FeatureCounts matrix: $GENE_MATRIX"

# ==========================
# STEP 3: StringTie individual GTFs
# ==========================
echo "===== STEP 3: StringTie per BAM ====="

IND_GTF_DIR="$HISAT_BASE/gtf_individual"
MERGED_GTF_DIR="$HISAT_BASE/gtf_merged"
BALLGOWN_DIR="$HISAT_BASE/ballgown"

for BAM in "${BAM_FILES[@]}"; do
    SAMPLE=$(basename "$BAM" .bam)

    if grep -Fxq "$SAMPLE" "$STRINGTIE_LOG"; then
        echo "Already processed: StringTie for $SAMPLE"
    else
        stringtie "$BAM" \
            -G "$REF_GTF" \
            -o "$IND_GTF_DIR/${SAMPLE}.gtf" \
            -p "$THREADS" \
            -l "$SAMPLE"

        echo "$SAMPLE" >> "$STRINGTIE_LOG"
    fi
done

# ==========================
# STEP 4: Merge StringTie GTFs
# ==========================
echo "===== STEP 4: Merge StringTie GTFs ====="

GTF_LIST="$MERGED_GTF_DIR/${STUDY_ID}_gtf_list.txt"
MERGED_GTF="$MERGED_GTF_DIR/${STUDY_ID}_stringtie_merged.gtf"

if [ -f "$MERGED_GTF" ]; then
    echo "Merged GTF already exists: $MERGED_GTF"
else
    > "$GTF_LIST"

    for BAM in "${BAM_FILES[@]}"; do
        SAMPLE=$(basename "$BAM" .bam)
        GTF_FILE="$IND_GTF_DIR/${SAMPLE}.gtf"

        if [ -f "$GTF_FILE" ]; then
            echo "$GTF_FILE" >> "$GTF_LIST"
        else
            echo "ERROR: Missing GTF file: $GTF_FILE"
            exit 1
        fi
    done

    stringtie --merge \
        -G "$REF_GTF" \
        -o "$MERGED_GTF" \
        "$GTF_LIST"
fi

# ==========================
# STEP 5: Ballgown prep
# ==========================
echo "===== STEP 5: Ballgown prep ====="

for BAM in "${BAM_FILES[@]}"; do
    SAMPLE=$(basename "$BAM" .bam)

    if grep -Fxq "$SAMPLE" "$BALLGOWN_LOG"; then
        echo "Already processed: Ballgown prep for $SAMPLE"
    else
        mkdir -p "$BALLGOWN_DIR/$SAMPLE"

        stringtie "$BAM" \
            -G "$MERGED_GTF" \
            -o "$BALLGOWN_DIR/$SAMPLE/${SAMPLE}.gtf" \
            -p "$THREADS" \
            -B \
            -e

        echo "$SAMPLE" >> "$BALLGOWN_LOG"
    fi
done

# ==========================
# STEP 6: Generate sample list for prepDE.py3
# ==========================
echo "===== STEP 6: Generate sample list for prepDE.py3 ====="

SAMPLE_LIST_PREPDE="$HISAT_BASE/counts/${STUDY_ID}_prepDE_sample_list.txt"
> "$SAMPLE_LIST_PREPDE"

for BAM in "${BAM_FILES[@]}"; do
    SAMPLE=$(basename "$BAM" .bam)
    GTF="$BALLGOWN_DIR/$SAMPLE/${SAMPLE}.gtf"

    if [ -f "$GTF" ]; then
        echo "$SAMPLE $(realpath "$GTF")" >> "$SAMPLE_LIST_PREPDE"
    else
        echo "ERROR: Missing Ballgown GTF: $GTF"
        exit 1
    fi
done
cat "$SAMPLE_LIST_PREPDE"
# ==========================
# STEP 7: prepDE.py3
# ==========================
echo "===== STEP 7: Run prepDE.py3 ====="

GENE_PREPDE="$GENE_DIR/${STUDY_ID}_gene_count_matrix_prepDE.csv"
TX_MATRIX="$TRANSCRIPT_DIR/${STUDY_ID}_transcript_count_matrix.csv"

if [ -f "$GENE_PREPDE" ] && [ -f "$TX_MATRIX" ]; then
    echo "prepDE outputs already exist"
else
    python3 $HISAT_BASE/prepDE.py3 \
        -i "$SAMPLE_LIST_PREPDE" \
        -g "$GENE_PREPDE" \
        -t "$TX_MATRIX"
fi

# ==========================
# STEP 8: DESeq2
# ==========================
echo "===== STEP 8: Run DESeq2 ====="

DESEQ2_DIR="$HISAT_BASE/deseq2/${STUDY_ID}"
mkdir -p "$DESEQ2_DIR"

if grep -Fxq "$STUDY_ID" "$DESEQ_LOG"; then
    echo "Already processed: DESeq2 for $STUDY_ID"
else
    Rscript $HISAT_BASE/scripts/deseq2_analysis_v2_2.R \
        "$GENE_MATRIX" \
        "$TX_MATRIX" \
        "$METADATA_FILE" \
        "$DESEQ2_DIR"

    echo "$STUDY_ID" >> "$DESEQ_LOG"
fi

# ==========================
# STEP 9: Ballgown analysis
# ==========================
echo "===== STEP 9: Ballgown analysis ====="

BALLGOWN_RESULTS_DIR="$BALLGOWN_DIR/results/${STUDY_ID}"
mkdir -p "$BALLGOWN_RESULTS_DIR"

BALLGOWN_ANALYSIS_ID="${STUDY_ID}_analysis"

if grep -Fxq "$BALLGOWN_ANALYSIS_ID" "$BALLGOWN_LOG"; then
    echo "Already processed: Ballgown analysis for $STUDY_ID"
else
    Rscript $HISAT_BASE/scripts/run_ballgown.R \
        "$BALLGOWN_DIR" \
        "$METADATA_FILE" \
        "$BALLGOWN_RESULTS_DIR"

    echo "$BALLGOWN_ANALYSIS_ID" >> "$BALLGOWN_LOG"
fi

echo "===== PIPELINE COMPLETED SUCCESSFULLY ====="
