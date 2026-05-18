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

get_cfg <- function(x, path, default = NULL) {
  cur <- x
  for (p in path) {
    if (is.null(cur[[p]])) return(default)
    cur <- cur[[p]]
  }
  cur
}

message("Writing outputs to: ", out_dir)

# ----------------------------
# Alt-EJ detection from MAF/RDS files
# ----------------------------

maf_dir <- cfg$paths$maf_dir

maf_files <- list.files(
  maf_dir,
  pattern = "\\.rds(\\.gz)?$",
  full.names = TRUE
)

if (length(maf_files) == 0) {
  warning("No MAF RDS files found in ", maf_dir)
} else {
  message("Found ", length(maf_files), " MAF/RDS files.")

  altej_annotated_all <- rbindlist(lapply(maf_files, function(path) {
    cohort <- sub("\\.rds(\\.gz)?$", "", basename(path))
    message("Alt-EJ: ", cohort)

    maf <- load_maf_rds(path)

    res <- classify_altej(
      maf,
      ref_genome,
      min_deletion_bp = get_cfg(cfg, c("thresholds", "altej_min_deletion_bp"), 6),
      min_microhomology_bp = get_cfg(cfg, c("thresholds", "altej_min_microhomology_bp"), 3),
      max_microhomology_bp = get_cfg(cfg, c("thresholds", "altej_max_microhomology_bp"), 25),
      flank_bp = get_cfg(cfg, c("thresholds", "altej_flank_bp"), 250)
    )

    res[, Cohort := cohort]
    res
  }), fill = TRUE)

  altej_candidates <- altej_annotated_all[AltEJ_Candidate == TRUE]
  altej_summary <- summarize_altej_by_sample(altej_annotated_all)

  fwrite(
    altej_annotated_all,
    file.path(out_dir, "PCAWG_AltEJ_annotated_indels.tsv"),
    sep = "\t"
  )

  fwrite(
    altej_candidates,
    file.path(out_dir, "PCAWG_AltEJ_like_deletions.tsv"),
    sep = "\t"
  )

  fwrite(
    altej_summary,
    file.path(out_dir, "PCAWG_AltEJ_summary_by_sample.tsv"),
    sep = "\t"
  )

  message("Alt-EJ candidate events: ", nrow(altej_candidates))
}

# ----------------------------
# SSA detection from PCAWG SV deletions
# ----------------------------

sv_file <- cfg$paths$sv_file

if (!file.exists(sv_file)) {
  warning("SV file not found: ", sv_file)
} else {
  message("SSA: ", sv_file)

  sv <- fread(sv_file)

  ssa_annotated <- classify_ssa_from_sv(
    sv,
    ref_genome,
    min_deletion_bp = get_cfg(cfg, c("thresholds", "ssa_min_deletion_bp"), 50),
    min_repeat_bp = get_cfg(cfg, c("thresholds", "ssa_min_repeat_bp"), 30),
    min_percent_identity = get_cfg(cfg, c("thresholds", "ssa_min_percent_identity"), 80),
    flank_bp = get_cfg(cfg, c("thresholds", "ssa_flank_bp"), 500),
    effect_col = get_cfg(cfg, c("columns", "sv_effect_col"), "effect"),
    progress_every = get_cfg(cfg, c("runtime", "progress_every"), 1000)
  )

  ssa_candidates <- ssa_annotated[SSA_Candidate == TRUE]
  ssa_summary <- summarize_ssa_by_sample(ssa_annotated)

  fwrite(
    ssa_annotated,
    file.path(out_dir, "PCAWG_SSA_annotated_large_deletions.tsv"),
    sep = "\t"
  )

  fwrite(
    ssa_candidates,
    file.path(out_dir, "PCAWG_SSA_like_deletions.tsv"),
    sep = "\t"
  )

  fwrite(
    ssa_summary,
    file.path(out_dir, "PCAWG_SSA_summary_by_sample.tsv"),
    sep = "\t"
  )

  message("SSA candidate events: ", nrow(ssa_candidates))
}

# ----------------------------
# Optional combined summary for downstream figures
# ----------------------------

altej_summary_file <- file.path(out_dir, "PCAWG_AltEJ_summary_by_sample.tsv")
ssa_summary_file <- file.path(out_dir, "PCAWG_SSA_summary_by_sample.tsv")

if (file.exists(altej_summary_file) && file.exists(ssa_summary_file)) {
  alt <- fread(altej_summary_file)
  ssa <- fread(ssa_summary_file)

  combined <- merge(ssa, alt, by = "SampleID", all = TRUE)

  for (col in c("Num_Large_DEL", "Num_SSA", "Num_DEL", "Num_AltEJ")) {
    if (col %in% names(combined)) {
      combined[is.na(get(col)), (col) := 0]
    }
  }

  if (all(c("Num_Large_DEL", "Num_DEL") %in% names(combined))) {
    combined[, Total_Deletions := Num_Large_DEL + Num_DEL]
  }

  fwrite(
    combined,
    file.path(out_dir, "PCAWG_composite_SSA_AltEJ_summary_by_sample.tsv"),
    sep = "\t"
  )
}
