# ==========================
# Load libraries
# ==========================
library(ballgown)
#library(RColorBrewer)


# ==========================
# Command line arguments
# ==========================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript run_ballgown.R <ballgown_dir> <sample_metadata.csv> <output_dir>")
}

ballgown_dir <- args[1]
sample_info_file <- args[2]
output_dir <- args[3]

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==========================
# Load metadata
# ==========================
sample_info <- read.csv(
  sample_info_file,
  header = TRUE,
  stringsAsFactors = FALSE
)


# ==========================
# Load Ballgown object
# ==========================
sample_paths <- file.path(ballgown_dir, sample_info$sample)
missing_paths <- sample_paths[!dir.exists(sample_paths)]

if (length(missing_paths) > 0) {
  stop("Missing Ballgown sample folders:\n", paste(missing_paths, collapse = "\n"))
}

bg <- ballgown(
  samples = sample_paths,
  meas = "all"
)

pData(bg) <- sample_info

# ==========================
# Differential transcript expression
# ==========================
results_transcripts <- stattest(
  bg,
  feature = "transcript",
  covariate = "condition",
  getFC = TRUE,
  meas = "FPKM"
)

results_transcripts <- data.frame(
  transcript_id = transcriptNames(bg),
  gene_id       = geneIDs(bg),
  gene_name     = geneNames(bg),
  results_transcripts,
  stringsAsFactors = FALSE
)



results_transcripts <- results_transcripts[order(results_transcripts$pval), ]

write.table(
  results_transcripts,
  file = file.path(output_dir, "ballgown_transcript_diffexp_results.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
# ==========================
# Save Ballgown object
# ==========================
saveRDS(
  bg,
  file = file.path(output_dir, "ballgown_object.rds")
)
