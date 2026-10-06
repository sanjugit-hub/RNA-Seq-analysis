#!/bin/bash

# ------------------------------------------------------------
# Purpose: Build Bowtie2 index from the prepared RSEM transcript
# reference and run RSEM gene-level TPM quantification for all
# testis paired-end RNA-seq samples.
# Author: SBM
# ------------------------------------------------------------

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------

TRIMMED_DIR="/home/group1/cgi-bin/rna_seq/hisat_deseq_stringtei/trimmed_reads"

RSEM_REF_DIR="/home/group1/cgi-bin/rna_seq/hisat_deseq_stringtei/RSEM/reference"

RSEM_REF="$RSEM_REF_DIR/GRCh38_116"

RSEM_OUT="/home/group1/cgi-bin/rna_seq/hisat_deseq_stringtei/RSEM/results"


# ------------------------------------------------------------
# CREATE OUTPUT DIRECTORIES
# ------------------------------------------------------------

mkdir -p "$RSEM_OUT"


# ------------------------------------------------------------
# CHECK REQUIRED PROGRAMS
# ------------------------------------------------------------

if ! command -v rsem-calculate-expression >/dev/null 2>&1; then
    echo "ERROR: RSEM is not available."
    exit 1
fi

if ! command -v bowtie2 >/dev/null 2>&1; then
    echo "ERROR: Bowtie2 is not available."
    exit 1
fi

if ! command -v bowtie2-build >/dev/null 2>&1; then
    echo "ERROR: bowtie2-build is not available."
    exit 1
fi


# ------------------------------------------------------------
# CHECK RSEM REFERENCE
# ------------------------------------------------------------

TRANSCRIPT_FASTA="${RSEM_REF}.transcripts.fa"

if [ ! -f "$TRANSCRIPT_FASTA" ]; then
    echo "ERROR: RSEM transcript FASTA not found:"
    echo "$TRANSCRIPT_FASTA"
    exit 1
fi

echo "RSEM transcript FASTA:"
echo "$TRANSCRIPT_FASTA"


# ------------------------------------------------------------
# BUILD BOWTIE2 INDEX
# ------------------------------------------------------------

echo "============================================================"
echo "Building Bowtie2 index"
echo "============================================================"

# Check whether Bowtie2 index already exists
if [ -f "${RSEM_REF}.1.bt2" ] || [ -f "${RSEM_REF}.1.bt2l" ]; then

    echo "Bowtie2 index already exists."
    echo "Skipping index construction."

else

    echo "Bowtie2 index not found."
    echo "Building index from:"
    echo "$TRANSCRIPT_FASTA"
    echo ""

    bowtie2-build \
        "$TRANSCRIPT_FASTA" \
        "$RSEM_REF"

    if [ $? -ne 0 ]; then
        echo "ERROR: Bowtie2 index construction failed."
        exit 1
    fi

    echo ""
    echo "Bowtie2 index construction completed."

fi


# ------------------------------------------------------------
# VERIFY BOWTIE2 INDEX
# ------------------------------------------------------------

echo ""
echo "Checking Bowtie2 index files..."

if [ -f "${RSEM_REF}.1.bt2" ]; then
    echo "Bowtie2 index detected:"
    ls -lh "${RSEM_REF}"*.bt2
elif [ -f "${RSEM_REF}.1.bt2l" ]; then
    echo "Large Bowtie2 index detected:"
    ls -lh "${RSEM_REF}"*.bt2l
else
    echo "ERROR: Bowtie2 index files were not created."
    exit 1
fi


# ------------------------------------------------------------
# MOVE TO TRIMMED READS DIRECTORY
# ------------------------------------------------------------

cd "$TRIMMED_DIR" || exit 1


# ------------------------------------------------------------
# RUN RSEM FOR ALL TESTIS SAMPLES
# ------------------------------------------------------------

echo ""
echo "============================================================"
echo "Starting RSEM quantification"
echo "============================================================"

for R1 in *testis*_R1.trimmed.fastq.gz; do

    [ -e "$R1" ] || continue

    SAMPLE=$(basename "$R1" _R1.trimmed.fastq.gz)

    R2="${SAMPLE}_R2.trimmed.fastq.gz"

    echo ""
    echo "============================================================"
    echo "Processing sample: $SAMPLE"
    echo "============================================================"
    echo "R1: $R1"
    echo "R2: $R2"

    # --------------------------------------------------------
    # Check R2
    # --------------------------------------------------------

    if [ ! -f "$R2" ]; then

        echo "WARNING: R2 not found for $SAMPLE"
        echo "Expected:"
        echo "$R2"
        echo "Skipping sample."

        continue
    fi


    # --------------------------------------------------------
    # Create sample output directories
    # --------------------------------------------------------

    mkdir -p "$RSEM_OUT/${SAMPLE}.temp"
    mkdir -p "$RSEM_OUT/${SAMPLE}.stat"


    # --------------------------------------------------------
    # Run RSEM
    # --------------------------------------------------------

    rsem-calculate-expression \
        --bowtie2 \
        --paired-end \
        -p 8 \
        "$R1" \
        "$R2" \
        "$RSEM_REF" \
        "$RSEM_OUT/$SAMPLE" \
        2>&1 | tee "$RSEM_OUT/${SAMPLE}.log"


    # --------------------------------------------------------
    # Check result
    # --------------------------------------------------------

    if [ ${PIPESTATUS[0]} -eq 0 ]; then

        echo ""
        echo "SUCCESS: $SAMPLE"

        if [ -f "$RSEM_OUT/${SAMPLE}.genes.results" ]; then
            echo "Gene-level result:"
            echo "$RSEM_OUT/${SAMPLE}.genes.results"
        else
            echo "WARNING: RSEM completed but genes.results was not found."
        fi

    else

        echo ""
        echo "ERROR: RSEM failed for $SAMPLE"
        echo "Check log:"
        echo "$RSEM_OUT/${SAMPLE}.log"

    fi

done


# ------------------------------------------------------------
# FINAL SUMMARY
# ------------------------------------------------------------

echo ""
echo "============================================================"
echo "RSEM quantification completed"
echo "============================================================"

echo ""
echo "Results directory:"
echo "$RSEM_OUT"

echo ""
echo "Gene-level RSEM result files:"
ls -lh "$RSEM_OUT"/*.genes.results 2>/dev/null

echo ""
echo "Log files:"
ls -lh "$RSEM_OUT"/*.log 2>/dev/null

echo ""
echo "DONE"
echo "============================================================"
