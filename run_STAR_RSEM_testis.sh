#!/bin/bash

# ------------------------------------------------------------
# Purpose: Build a STAR-based RSEM reference and perform
#          gene-level RSEM quantification for all testis
#          paired-end RNA-seq samples using existing trimmed
#          FASTQ files.
# Author: SBM
# ------------------------------------------------------------

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------

BASE="/home/group1/cgi-bin/rna_seq/hisat_deseq_stringtei"

TRIMMED_DIR="$BASE/trimmed_reads"

# New STAR-based RSEM reference
RSEM_REF_DIR="$BASE/RSEM/star_reference"
RSEM_REF="$RSEM_REF_DIR/GRCh38_116"

# RSEM quantification results
RSEM_OUT="$BASE/RSEM/star_results"

# ------------------------------------------------------------
# INPUT REFERENCE FILES
# ------------------------------------------------------------
# IMPORTANT:
# Use the SAME GRCh38 release 116 FASTA and GTF that were
# used for your existing GRCh38_116 RSEM reference.

GENOME_FA="/home/group1/cgi-bin/rna_seq/hisat_deseq_stringtei/Reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
ANNOTATION_GTF="/home/group1/cgi-bin/rna_seq/hisat_deseq_stringtei/Reference_genome/Homo_sapiens.GRCh38.116.gtf"

# ------------------------------------------------------------
# CREATE DIRECTORIES
# ------------------------------------------------------------

mkdir -p "$RSEM_REF_DIR"
mkdir -p "$RSEM_OUT"

# ------------------------------------------------------------
# CHECK PROGRAMS
# ------------------------------------------------------------

if ! command -v STAR >/dev/null 2>&1; then
    echo "ERROR: STAR is not installed."
    exit 1
fi

if ! command -v rsem-prepare-reference >/dev/null 2>&1; then
    echo "ERROR: RSEM is not installed."
    exit 1
fi

if ! command -v rsem-calculate-expression >/dev/null 2>&1; then
    echo "ERROR: RSEM is not installed."
    exit 1
fi

# ------------------------------------------------------------
# CHECK INPUT FILES
# ------------------------------------------------------------

if [ ! -f "$GENOME_FA" ]; then
    echo "ERROR: Genome FASTA not found:"
    echo "$GENOME_FA"
    exit 1
fi

if [ ! -f "$ANNOTATION_GTF" ]; then
    echo "ERROR: GTF annotation not found:"
    echo "$ANNOTATION_GTF"
    exit 1
fi

if [ ! -d "$TRIMMED_DIR" ]; then
    echo "ERROR: Trimmed reads directory not found:"
    echo "$TRIMMED_DIR"
    exit 1
fi

# ------------------------------------------------------------
# DISPLAY VERSIONS
# ------------------------------------------------------------

echo "============================================================"
echo "STAR version:"
STAR --version

echo ""
echo "RSEM version:"
rsem-calculate-expression --version
echo "============================================================"

# ============================================================
# STEP 1: BUILD STAR-BASED RSEM REFERENCE
# ============================================================

echo ""
echo "============================================================"
echo "STEP 1: Building STAR-based RSEM reference"
echo "============================================================"

echo ""
echo "Genome:"
echo "$GENOME_FA"

echo ""
echo "GTF:"
echo "$ANNOTATION_GTF"

echo ""
echo "RSEM reference:"
echo "$RSEM_REF"

# ------------------------------------------------------------
# Check if reference already exists
# ------------------------------------------------------------

if [ -f "${RSEM_REF}.transcripts.fa" ]; then

    echo ""
    echo "✔ STAR-based RSEM reference already exists."
    echo "Skipping reference construction."

else

    echo ""
    echo "⚡ Building RSEM reference using STAR..."

    rsem-prepare-reference \
        --star \
        --gtf "$ANNOTATION_GTF" \
        "$GENOME_FA" \
        "$RSEM_REF"

    if [ $? -ne 0 ]; then
        echo ""
        echo "ERROR: RSEM STAR reference construction failed."
        exit 1
    fi

    echo ""
    echo "✅ STAR-based RSEM reference created."

fi

# ------------------------------------------------------------
# VERIFY RSEM REFERENCE
# ------------------------------------------------------------

if [ ! -f "${RSEM_REF}.transcripts.fa" ]; then

    echo ""
    echo "ERROR: RSEM STAR reference is incomplete."
    exit 1

fi

# ============================================================
# STEP 2: FIND TESTIS SAMPLES
# ============================================================

echo ""
echo "============================================================"
echo "STEP 2: Finding testis samples"
echo "============================================================"

cd "$TRIMMED_DIR" || exit 1

shopt -s nullglob

R1_FILES=(*testis*_R1.trimmed.fastq.gz)

if [ ${#R1_FILES[@]} -eq 0 ]; then

    echo ""
    echo "ERROR: No testis R1 files found."

    echo "Expected pattern:"
    echo "*testis*_R1.trimmed.fastq.gz"

    exit 1

fi

echo ""
echo "Found ${#R1_FILES[@]} testis R1 files."

# ============================================================
# STEP 3: RSEM + STAR GENE-LEVEL QUANTIFICATION
# ============================================================

echo ""
echo "============================================================"
echo "STEP 3: RSEM + STAR quantification"
echo "============================================================"

for R1 in "${R1_FILES[@]}"; do

    # --------------------------------------------------------
    # Get sample name
    # --------------------------------------------------------

    SAMPLE=$(basename "$R1" _R1.trimmed.fastq.gz)

    R2="${SAMPLE}_R2.trimmed.fastq.gz"

    echo ""
    echo "============================================================"
    echo "Processing: $SAMPLE"
    echo "============================================================"

    echo "R1:"
    echo "$R1"

    echo "R2:"
    echo "$R2"

    # --------------------------------------------------------
    # Check R2
    # --------------------------------------------------------

    if [ ! -f "$R2" ]; then

        echo ""
        echo "WARNING: R2 not found for $SAMPLE"
        echo "Expected:"
        echo "$R2"
        echo "Skipping..."

        continue

    fi

    # --------------------------------------------------------
    # Skip completed samples
    # --------------------------------------------------------

    if [ -f "$RSEM_OUT/${SAMPLE}.genes.results" ]; then

        echo ""
        echo "✔ Gene-level RSEM result already exists."
        echo "Skipping $SAMPLE."

        continue

    fi

    # --------------------------------------------------------
    # Run RSEM using STAR
    # --------------------------------------------------------

    echo ""
    echo "⚡ Running RSEM using STAR..."

    rsem-calculate-expression \
        --star \
        --star-gzipped-read-file \
        --paired-end \
        -p 8 \
        "$R1" \
        "$R2" \
        "$RSEM_REF" \
        "$RSEM_OUT/$SAMPLE" \
        2>&1 | tee "$RSEM_OUT/${SAMPLE}.log"

    RSEM_STATUS=${PIPESTATUS[0]}

    # --------------------------------------------------------
    # Check result
    # --------------------------------------------------------

    if [ $RSEM_STATUS -eq 0 ]; then

        echo ""
        echo "✅ SUCCESS: $SAMPLE"

        if [ -f "$RSEM_OUT/${SAMPLE}.genes.results" ]; then

            echo ""
            echo "Gene-level result:"
            echo "$RSEM_OUT/${SAMPLE}.genes.results"

        else

            echo ""
            echo "WARNING: genes.results not found."

        fi

    else

        echo ""
        echo "❌ ERROR: RSEM failed for $SAMPLE"
        echo "Check:"
        echo "$RSEM_OUT/${SAMPLE}.log"

    fi

done

# ============================================================
# FINAL SUMMARY
# ============================================================

echo ""
echo "============================================================"
echo "RSEM + STAR PIPELINE COMPLETED"
echo "============================================================"

echo ""
echo "STAR-based RSEM reference:"
echo "$RSEM_REF"

echo ""
echo "Gene-level results:"
echo "$RSEM_OUT"

echo ""
echo "Generated gene-level files:"
ls -lh "$RSEM_OUT"/*.genes.results 2>/dev/null

echo ""
echo "DONE"
echo "============================================================"
