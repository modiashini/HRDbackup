#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(BSgenome.Hsapiens.UCSC.hg19)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- ifelse(length(args) >= 1, args[[1]], "config/config.yml")
cfg <- yaml::read_yaml(config_path)
source("R/repair_scar_detection.R")

out_dir <- cfg$paths$output_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ref_genome <- BSgenome.Hsapiens.UCSC.hg19

# ---------- Alt-EJ from all MAF/RDS files ----------
maf_files <- list.files(cfg$paths$maf_dir, pattern = "\\.rds(\\.gz)?$", full.names = TRUE)
if (length(maf_files) == 0) warning("No MAF RDS files found in ", cfg$paths$maf_dir)

altej_events <- rbindlist(lapply(maf_files, function(path) {
  message("Alt-EJ: ", basename(path))
  maf <- load_maf_rds(path)
  cohort <- sub("\\.rds(\\.gz)?$", "", basename(path))
  res <- classify_altej(
    maf, ref_genome,
    min_deletion_bp = cfg$thresholds$altej_min_deletion_bp,
    min_microhomology_bp = cfg$thresholds$altej_min_microhomology_bp,
    max_microhomology_bp = cfg$thresholds$altej_max_microhomology_bp
  )
  res[, Cohort := cohort]
  res
}), fill = TRUE)

if (nrow(altej_events) > 0) {
  altej_summary <- summarize_altej_by_sample(altej_events)
  fwrite(altej_events, file.path(out_dir, "PCAWG_AltEJ_like_deletions.tsv"), sep = "\t")
  fwrite(altej_summary, file.path(out_dir, "PCAWG_AltEJ_summary_by_sample.tsv"), sep = "\t")
}

# ---------- SSA from SV deletions ----------
if (file.exists(cfg$paths$sv_file)) {
  message("SSA: ", cfg$paths$sv_file)
  sv <- fread(cfg$paths$sv_file)
  ssa_events <- classify_ssa_from_sv(
    sv, ref_genome,
    min_deletion_bp = cfg$thresholds$ssa_min_deletion_bp,
    min_repeat_bp = cfg$thresholds$ssa_min_repeat_bp,
    min_percent_identity = cfg$thresholds$ssa_min_percent_identity,
    flank_bp = cfg$thresholds$ssa_flank_bp
  )
  ssa_summary <- summarize_ssa_by_sample(ssa_events)
  fwrite(ssa_events, file.path(out_dir, "PCAWG_SSA_like_deletions.tsv"), sep = "\t")
  fwrite(ssa_summary, file.path(out_dir, "PCAWG_SSA_summary_by_sample.tsv"), sep = "\t")
} else {
  warning("SV file not found: ", cfg$paths$sv_file)
}
