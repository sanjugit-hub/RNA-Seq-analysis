#!/usr/bin/env bash

###############################################################################
# RNA-seq Transcript Quantification Pipeline
#
# Author      : Sanjana B M
# Version     : 2.1  (fixed: pipeline no longer aborts after the first sample)
#
# Workflow
#
# samples.csv
#      │
#      ▼
# Group SRRs by GSM
#      │
#      ▼
# Process ONE GSM completely
#      │
#      ├── fasterq-dump
#      ├── Trim Galore
#      ├── Delete raw FASTQ
#      ├── Merge trimmed FASTQs
#      ├── Run Kallisto
#      ├── Delete downloaded SRRs
#      └── Next GSM
#
###############################################################################

set -Eeuo pipefail

###############################################################################
# Activate Conda Environment
###############################################################################

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

###############################################################################
# USER CONFIGURATION
###############################################################################

PROJECT_ROOT="/home/ibab/RNAseq_pipeline"

CSV="$PROJECT_ROOT/samples.csv"

DOWNLOAD_DIR="$PROJECT_ROOT/downloads"

RAW_DIR="$PROJECT_ROOT/raw_fastq"

RESULT_DIR="$PROJECT_ROOT/kallisto/results"

LOG_DIR="$PROJECT_ROOT/logs"

TMP_DIR="$PROJECT_ROOT/tmp"

INDEX="$PROJECT_ROOT/kallisto/index/Hs_cdna.index"

THREADS=8

BOOTSTRAPS=50

###############################################################################
# DIRECTORY SETUP
###############################################################################

mkdir -p "$DOWNLOAD_DIR"
mkdir -p "$RAW_DIR"
mkdir -p "$RESULT_DIR"

LOG_DIR="$PROJECT_ROOT/logs"

mkdir -p "$LOG_DIR"

SUCCESS_LOG="$LOG_DIR/success.log"

FAILED_LOG="$LOG_DIR/failed.log"

SKIPPED_LOG="$LOG_DIR/skipped.log"

PIPELINE_LOG="$LOG_DIR/pipeline.log"

touch "$SUCCESS_LOG"
touch "$FAILED_LOG"
touch "$SKIPPED_LOG"
touch "$PIPELINE_LOG"

###############################################################################
# CHECK INPUT FILES
###############################################################################

if [[ ! -f "$CSV" ]]; then
    echo "ERROR : samples.csv not found."
    exit 1
fi

if [[ ! -f "$INDEX" ]]; then
    echo "ERROR : Kallisto index not found."
    exit 1
fi

###############################################################################
# CHECK REQUIRED SOFTWARE
###############################################################################

PROGRAMS=(
    fasterq-dump
    trim_galore
    kallisto
    fastqc
)

for PROGRAM in "${PROGRAMS[@]}"
do

    if ! command -v "$PROGRAM" >/dev/null 2>&1
    then
        echo "ERROR : $PROGRAM not found."
        exit 1
    fi

done

###############################################################################
# LOGGING FUNCTIONS
###############################################################################

timestamp() {

    date +"%Y-%m-%d %H:%M:%S"

}

log() {

    local MESSAGE="$1"

    echo "[$(timestamp)] $MESSAGE"

    echo "[$(timestamp)] $MESSAGE" >> "$PIPELINE_LOG"

}

success() {

    echo "$1" >> "$SUCCESS_LOG"

}

failed() {

    echo "$1" >> "$FAILED_LOG"

}

skipped() {

    echo "$1" >> "$SKIPPED_LOG"

}

###############################################################################
# ERROR HANDLER
#
# NOTE: This trap is intentionally NOT allowed to kill the whole pipeline
# anymore. process_sample() is always invoked inside its own subshell (see
# MAIN section below), so when this trap calls `exit 1`, it only terminates
# that subshell / that one sample - the outer while-loop keeps running and
# moves on to the next GSM.
###############################################################################

error_exit() {

    log "Error encountered while processing this sample. Aborting THIS sample only."

    exit 1

}

trap error_exit ERR

###############################################################################
# DISPLAY SETTINGS
###############################################################################

echo

echo "=========================================================="

echo "RNA-seq Transcript Quantification Pipeline"

echo "=========================================================="

echo "Project Root : $PROJECT_ROOT"

echo "Metadata     : $CSV"

echo "Downloads    : $DOWNLOAD_DIR"

echo "Raw FASTQ    : $RAW_DIR"

echo "Results      : $RESULT_DIR"

echo "Index        : $INDEX"

echo "Threads      : $THREADS"

echo "Bootstraps   : $BOOTSTRAPS"

echo

log "Pipeline Started"

###############################################################################
# GLOBAL VARIABLES
###############################################################################

CURRENT_GSM=""

CURRENT_TISSUE=""

declare -a SRR_LIST=()

###############################################################################
# FUNCTION : process_sample
###############################################################################

process_sample() {

    local SAMPLE="$1"
    local TISSUE="$2"

    shift 2

    local SRRS=("$@")

    local SAMPLE_DIR="$RAW_DIR/$SAMPLE"

    local OUTDIR="$RESULT_DIR/${SAMPLE}_kal_results"

    mkdir -p "$SAMPLE_DIR"

    mkdir -p "$OUTDIR"

    log "------------------------------------------------------------"
    log "Processing Sample : $SAMPLE"
    log "Tissue            : $TISSUE"
    log "SRRs              : ${SRRS[*]}"

    ###########################################################################
    # Resume support
    ###########################################################################

    if [[ -f "$OUTDIR/abundance.tsv" ]]; then

        log "Sample already completed. Skipping."

        skipped "$SAMPLE"

        return

    fi

    ###########################################################################
    # Process each SRR
    ###########################################################################

    for SRR in "${SRRS[@]}"
    do

        log "Processing SRR : $SRR"

        SRR_DIR="$SAMPLE_DIR/$SRR"

        mkdir -p "$SRR_DIR"

        DOWNLOAD_FILE="$DOWNLOAD_DIR/$SRR"

        #######################################################################
        # Check downloaded SRR
        #######################################################################

        if [[ ! -e "$DOWNLOAD_FILE" ]]; then

            log "Downloaded SRR not found : $SRR"

            failed "$SRR"

            continue

        fi

        #######################################################################
        # Skip if already trimmed
        #######################################################################

        if [[ -f "$SRR_DIR/${SRR}_1_val_1.fq" && \
              -f "$SRR_DIR/${SRR}_2_val_2.fq" ]]; then

            log "Trimmed files already exist."

            continue

        fi

        #######################################################################
        # fasterq-dump
        #######################################################################

        log "Running fasterq-dump..."

        if ! fasterq-dump \
            "$DOWNLOAD_FILE" \
            --split-files \
            --threads "$THREADS" \
            --outdir "$SRR_DIR"; then

            log "fasterq-dump failed for $SRR."

            failed "$SRR"

            continue

        fi

        #######################################################################
        # Verify FASTQ
        #######################################################################

        if [[ ! -f "$SRR_DIR/${SRR}_1.fastq" ]]; then

            log "FASTQ R1 missing."

            failed "$SRR"

            continue

        fi

        if [[ ! -f "$SRR_DIR/${SRR}_2.fastq" ]]; then

            log "FASTQ R2 missing."

            failed "$SRR"

            continue

        fi

        log "FASTQ generated successfully."

        #######################################################################
        # Run Trim Galore
        #######################################################################

        log "Running Trim Galore..."

        if ! trim_galore \
            --paired \
            --fastqc \
            -j "$THREADS" \
            -o "$SRR_DIR" \
            "$SRR_DIR/${SRR}_1.fastq" \
            "$SRR_DIR/${SRR}_2.fastq"; then

            log "Trim Galore failed for $SRR."

            failed "$SRR"

            continue

        fi

        #######################################################################
        # Verify Trimmed FASTQs
        #######################################################################

        if [[ ! -f "$SRR_DIR/${SRR}_1_val_1.fq" ]]; then

            log "Trimmed R1 missing."

            failed "$SRR"

            continue

        fi

        if [[ ! -f "$SRR_DIR/${SRR}_2_val_2.fq" ]]; then

            log "Trimmed R2 missing."

            failed "$SRR"

            continue

        fi

        log "Trim Galore completed."

        #######################################################################
        # Remove raw FASTQs
        #######################################################################

        rm -f "$SRR_DIR/${SRR}_1.fastq"
        rm -f "$SRR_DIR/${SRR}_2.fastq"

        log "Temporary FASTQs removed."

    done

    ###########################################################################
    # Merge trimmed FASTQs
    ###########################################################################

    MERGED_R1="$SAMPLE_DIR/${SAMPLE}_R1.fq"
    MERGED_R2="$SAMPLE_DIR/${SAMPLE}_R2.fq"

    rm -f "$MERGED_R1"
    rm -f "$MERGED_R2"

    declare -a R1_FILES=()
    declare -a R2_FILES=()

    for SRR in "${SRRS[@]}"
    do

        TRIM1="$SAMPLE_DIR/$SRR/${SRR}_1_val_1.fq"
        TRIM2="$SAMPLE_DIR/$SRR/${SRR}_2_val_2.fq"

        if [[ -f "$TRIM1" && -f "$TRIM2" ]]; then

            R1_FILES+=("$TRIM1")
            R2_FILES+=("$TRIM2")

        fi

    done

    ###########################################################################
    # Check merged input
    ###########################################################################

    if [[ ${#R1_FILES[@]} -eq 0 ]]; then

        log "No trimmed FASTQs available."

        failed "$SAMPLE"

        return

    fi

    ###########################################################################
    # Single SRR
    ###########################################################################

    if [[ ${#R1_FILES[@]} -eq 1 ]]; then

        cp "${R1_FILES[0]}" "$MERGED_R1"
        cp "${R2_FILES[0]}" "$MERGED_R2"

        log "Single SRR detected."

    ###########################################################################
    # Multiple SRRs
    ###########################################################################

    else

        log "Merging ${#R1_FILES[@]} SRRs..."

        cat "${R1_FILES[@]}" > "$MERGED_R1"

        cat "${R2_FILES[@]}" > "$MERGED_R2"

    fi

    ###########################################################################
    # Verify merged FASTQs
    ###########################################################################

    if [[ ! -f "$MERGED_R1" || ! -f "$MERGED_R2" ]]; then

        log "Merged FASTQs missing."

        failed "$SAMPLE"

        return

    fi

    log "Merged FASTQs created."

    ###########################################################################
    # Run Kallisto
    ###########################################################################

    log "Running Kallisto..."

    if ! kallisto quant \
        -i "$INDEX" \
        -o "$OUTDIR" \
        -t "$THREADS" \
        -b "$BOOTSTRAPS" \
        "$MERGED_R1" \
        "$MERGED_R2"; then

        log "Kallisto command failed for $SAMPLE."

        failed "$SAMPLE"

        return

    fi

    ###########################################################################
    # Verify Kallisto
    ###########################################################################

    if [[ ! -f "$OUTDIR/abundance.tsv" ]]; then

        log "Kallisto failed."

        failed "$SAMPLE"

        return

    fi

    log "Kallisto completed successfully."

    success "$SAMPLE"

    ###########################################################################
    # Cleanup downloaded SRRs
    ###########################################################################

    log "Cleaning temporary files..."

    for SRR in "${SRRS[@]}"
    do

        rm -rf "$DOWNLOAD_DIR/$SRR"

    done

    ###########################################################################
    # Remove merged FASTQs
    ###########################################################################

    rm -f "$MERGED_R1"
    rm -f "$MERGED_R2"

    ###########################################################################
    # Remove SRR directories
    ###########################################################################

    for SRR in "${SRRS[@]}"
    do

        SRR_DIR="$SAMPLE_DIR/$SRR"

        ###############################################################
        # Remove trimmed FASTQs
        ###############################################################

        rm -f "$SRR_DIR/${SRR}_1_val_1.fq"
        rm -f "$SRR_DIR/${SRR}_2_val_2.fq"

        ###############################################################
        # Remove FastQC reports
        ###############################################################

        rm -f "$SRR_DIR"/*_fastqc.html
        rm -f "$SRR_DIR"/*_fastqc.zip

        ###############################################################
        # Remove Trim Galore reports
        ###############################################################

        rm -f "$SRR_DIR"/*_trimming_report.txt

        ###############################################################
        # Remove any leftover FASTQs
        ###############################################################

        rm -f "$SRR_DIR/${SRR}_1.fastq"
        rm -f "$SRR_DIR/${SRR}_2.fastq"

        ###############################################################
        # Remove SRR directory if empty
        ###############################################################

        rmdir "$SRR_DIR" 2>/dev/null || true

    done

    ###########################################################################
    # Remove GSM directory if empty
    ###########################################################################

    rmdir "$SAMPLE_DIR" 2>/dev/null || true

    log "Cleanup completed."

    log "Sample completed : $SAMPLE"

    log "------------------------------------------------------------"

}

###############################################################################
# WRAPPER : run_sample
#
# Runs process_sample() inside its OWN subshell and never lets a failure
# escape and kill the outer loop. This is the key fix: previously, if any
# command inside process_sample failed, the global ERR trap called `exit 1`
# in the SAME shell that was running the whole while-loop over samples.csv,
# which terminated the entire pipeline after just one sample. Now, a failure
# is contained to that one sample, logged, and the loop continues.
###############################################################################

run_sample() {

    local SAMPLE="$1"

    if ( process_sample "$@" ); then
        log "run_sample: $SAMPLE finished (see logs above for success/skip/fail detail)."
    else
        log "run_sample: $SAMPLE hit an unrecoverable error. Marking as FAILED and moving to next sample."
        failed "$SAMPLE"
    fi

}

###############################################################################
# MAIN
###############################################################################

log "Reading samples.csv..."

CURRENT_GSM=""
CURRENT_TISSUE=""

declare -a SRR_LIST=()

# NOTE: using process substitution ( < <(...) ) instead of a pipe
# ( ... | while ... ) so the loop runs in the *current* shell rather than a
# subshell. This avoids surprises with variable scoping and makes it crystal
# clear that only run_sample()'s internal subshell is isolated - not the loop
# driving it.

while IFS=',' read -r GSM SRR TISSUE
do

    ###############################################################
    # Remove leading/trailing spaces
    ###############################################################

    GSM=$(echo "$GSM" | xargs)

    SRR=$(echo "$SRR" | xargs)

    TISSUE=$(echo "$TISSUE" | xargs)

    ###############################################################
    # If GSM missing, use SRR as Sample ID
    ###############################################################

    if [[ -z "$GSM" ]]; then
        GSM="$SRR"
    fi

    ###############################################################
    # First sample
    ###############################################################

    if [[ -z "$CURRENT_GSM" ]]; then

        CURRENT_GSM="$GSM"

        CURRENT_TISSUE="$TISSUE"

    fi

    ###############################################################
    # Same GSM
    ###############################################################

    if [[ "$GSM" == "$CURRENT_GSM" ]]; then

        SRR_LIST+=("$SRR")

    ###############################################################
    # New GSM encountered
    ###############################################################

    else

        run_sample \
            "$CURRENT_GSM" \
            "$CURRENT_TISSUE" \
            "${SRR_LIST[@]}"

        ###########################################################
        # Reset
        ###########################################################

        CURRENT_GSM="$GSM"

        CURRENT_TISSUE="$TISSUE"

        SRR_LIST=("$SRR")

    fi

done < <(tail -n +2 "$CSV")

###############################################################################
# VERY IMPORTANT
# Process the LAST GSM
###############################################################################

if [[ -n "$CURRENT_GSM" ]]; then

    run_sample \
        "$CURRENT_GSM" \
        "$CURRENT_TISSUE" \
        "${SRR_LIST[@]}"

fi

###############################################################################
# SUMMARY
###############################################################################

SUCCESS_COUNT=$(sort -u "$SUCCESS_LOG" | wc -l)

FAILED_COUNT=$(sort -u "$FAILED_LOG" | wc -l)

SKIPPED_COUNT=$(sort -u "$SKIPPED_LOG" | wc -l)

TOTAL=$((SUCCESS_COUNT + FAILED_COUNT + SKIPPED_COUNT))

echo
echo "=========================================================="

echo "RNA-seq Pipeline Finished"

echo "=========================================================="

echo "Total Samples : $TOTAL"

echo "Successful    : $SUCCESS_COUNT"

echo "Skipped       : $SKIPPED_COUNT"

echo "Failed        : $FAILED_COUNT"

echo

echo "Success Log : $SUCCESS_LOG"

echo "Failed Log  : $FAILED_LOG"

echo "Skipped Log : $SKIPPED_LOG"

echo

log "Pipeline Finished"

exit 0
