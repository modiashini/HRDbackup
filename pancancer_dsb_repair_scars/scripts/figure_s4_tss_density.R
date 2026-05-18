#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(ggplot2)
})

cfg_path <- commandArgs(trailingOnly = TRUE)
if (length(cfg_path) == 0) cfg_path <- "config/config.yml"
cfg <- yaml::read_yaml(cfg_path)

get_cfg <- function(x, path, default = NULL) {
  cur <- x
  for (p in path) {
    if (is.null(cur[[p]])) return(default)
    cur <- cur[[p]]
  }
  cur
}

output_dir <- get_cfg(cfg, c("paths", "output_dir"), "output")
figure_dir <- get_cfg(cfg, c("paths", "supplement_dir"),
                      file.path(get_cfg(cfg, c("paths", "figure_dir"), "figures"), "supplement"))
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

ssa_tsv <- get_cfg(cfg, c("paths", "ssa_events"), file.path(output_dir, "PCAWG_SSA_like_deletions.tsv"))
alt_tsv <- get_cfg(cfg, c("paths", "altej_events"), file.path(output_dir, "PCAWG_AltEJ_like_deletions.tsv"))
hrd_csv <- get_cfg(cfg, c("paths", "hrd_file"))
donor_file <- get_cfg(cfg, c("paths", "donor_file"))

window_kb <- get_cfg(cfg, c("supplement", "figS4_window_kb"), 20)
bin_kb <- get_cfg(cfg, c("supplement", "figS4_bin_kb"), 0.5)
ymax <- get_cfg(cfg, c("supplement", "figS4_ymax"), 3)
keep_cancers <- unlist(get_cfg(cfg, c("supplement", "figS4_cohorts"), list("Lymph-BNHL", "Lymph-CLL")))

norm_chr <- function(x) {
  x <- as.character(x)
  x[x %in% c("MT", "M")] <- "chrM"
  x[!grepl("^chr", x)] <- paste0("chr", x[!grepl("^chr", x)])
  x
}

standardize_events <- function(dt) {
  setnames(dt, "Tumor_Sample_Barcode", "SampleID", skip_absent = TRUE)
  setnames(dt, "Sample", "SampleID", skip_absent = TRUE)
  setnames(dt, "Start_Position", "Start", skip_absent = TRUE)
  setnames(dt, "End_Position", "End", skip_absent = TRUE)

  needed <- c("SampleID", "Chromosome", "Start", "End")
  missing <- setdiff(needed, names(dt))
  if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

  dt[, Chromosome := norm_chr(Chromosome)]
  dt[, Start := as.integer(Start)]
  dt[, End := as.integer(End)]
  dt <- dt[is.finite(Start) & is.finite(End)]
  dt
}

add_hrd_bins <- function(dt, hrd_csv) {
  hrd <- fread(hrd_csv)
  setnames(hrd, "HRD-sum", "HRD_sum", skip_absent = TRUE)
  hrd[, HRD_sum := as.numeric(HRD_sum)]

  dt <- merge(dt, hrd[, .(SampleID, HRD_sum)], by = "SampleID", all.x = TRUE)
  dt <- dt[!is.na(HRD_sum)]
  dt[, HRD_bin := fifelse(HRD_sum <= 21, "HRD 0-21",
                          fifelse(HRD_sum <= 42, "HRD 21-42", "HRD > 42"))]
  dt[, HRD_bin := factor(HRD_bin, levels = c("HRD 0-21", "HRD 21-42", "HRD > 42"))]
  dt
}

add_cohort <- function(dt, donor_file) {
  if (is.null(donor_file) || !file.exists(donor_file)) return(dt)
  donor <- fread(donor_file)
  if (!all(c("icgc_specimen_id", "histology_abbreviation") %in% names(donor))) return(dt)

  donor <- donor[, .(
    SampleID = as.character(icgc_specimen_id),
    Cohort = as.character(histology_abbreviation)
  )]
  unique(donor)
  merge(dt, donor, by = "SampleID", all.x = TRUE)
}

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
tx_gr <- transcripts(txdb)
tss_gr <- resize(tx_gr, width = 1, fix = "start")
tss_gr <- keepStandardChromosomes(tss_gr, pruning.mode = "coarse")

signed_dist_to_nearest_tss_bp <- function(bp, tss_gr) {
  if (length(bp) == 0) return(numeric(0))
  nn <- distanceToNearest(bp, tss_gr, ignore.strand = TRUE)
  if (length(nn) == 0) return(numeric(0))
  q <- queryHits(nn)
  s <- subjectHits(nn)
  start(bp)[q] - start(tss_gr)[s]
}

signed_dist_event_closest_end <- function(dt, tss_gr) {
  dt <- as.data.table(copy(dt))
  dt[, Chromosome := norm_chr(Chromosome)]

  gr_l <- GRanges(seqnames = dt$Chromosome, ranges = IRanges(start = dt$Start, width = 1))
  gr_r <- GRanges(seqnames = dt$Chromosome, ranges = IRanges(start = dt$End, width = 1))

  gr_l <- keepStandardChromosomes(gr_l, pruning.mode = "coarse")
  gr_r <- keepStandardChromosomes(gr_r, pruning.mode = "coarse")

  d_l <- signed_dist_to_nearest_tss_bp(gr_l, tss_gr)
  d_r <- signed_dist_to_nearest_tss_bp(gr_r, tss_gr)

  n <- min(length(d_l), length(d_r))
  if (n == 0) return(numeric(0))

  d_l <- d_l[seq_len(n)]
  d_r <- d_r[seq_len(n)]

  ifelse(abs(d_l) <= abs(d_r), d_l, d_r)
}

density_df_from_dist <- function(d_bp, group, window_kb, bin_kb) {
  d_kb <- d_bp / 1000
  d_kb <- d_kb[is.finite(d_kb) & d_kb >= -window_kb & d_kb <= window_kb]

  breaks <- seq(-window_kb, window_kb, by = bin_kb)
  mids <- head(breaks, -1) + diff(breaks) / 2
  counts <- hist(d_kb, breaks = breaks, plot = FALSE)$counts

  data.table(
    Dist_kb = mids,
    Density = counts / sum(counts, na.rm = TRUE),
    Group = group
  )
}

plot_tss_density <- function(dt, title, out_pdf, restrict_cohorts = NULL) {
  if (!is.null(restrict_cohorts) && "Cohort" %in% names(dt)) {
    dt <- dt[Cohort %in% restrict_cohorts]
  }

  bins <- levels(dt$HRD_bin)

  curves <- rbindlist(lapply(bins, function(bin) {
    sub <- dt[HRD_bin == bin]
    d <- signed_dist_event_closest_end(sub, tss_gr)
    x <- density_df_from_dist(d, bin, window_kb, bin_kb)
    x[, HRD_bin := bin]
    x
  }), fill = TRUE)

  curves[, HRD_bin := factor(HRD_bin, levels = bins)]

  p <- ggplot(curves, aes(x = Dist_kb, y = Density, color = HRD_bin)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, color = "gray50") +
    geom_line(linewidth = 1.05, na.rm = TRUE) +
    coord_cartesian(xlim = c(-window_kb, window_kb), ylim = c(0, ymax)) +
    scale_color_manual(
      values = c("HRD 0-21" = "#c6dbef", "HRD 21-42" = "#6baed6", "HRD > 42" = "#08306b"),
      breaks = bins
    ) +
    labs(
      title = title,
      x = "Distance to nearest TSS (kb)",
      y = "Event density",
      color = NULL
    ) +
    theme_classic(base_size = 14)

  ggsave(out_pdf, p, width = 7.5, height = 5)
}

ssa <- fread(ssa_tsv)
ssa <- standardize_events(ssa)
ssa <- add_cohort(ssa, donor_file)
ssa <- add_hrd_bins(ssa, hrd_csv)

alt <- fread(alt_tsv)
alt <- standardize_events(alt)
alt <- add_cohort(alt, donor_file)
alt <- add_hrd_bins(alt, hrd_csv)

plot_tss_density(
  ssa,
  "SSA-like event density vs nearest TSS, lymphoid cohorts",
  file.path(figure_dir, "figS4A_ssa_tss_density_lymphoid.pdf"),
  restrict_cohorts = keep_cancers
)

plot_tss_density(
  alt,
  "Alt-EJ-like event density vs nearest TSS",
  file.path(figure_dir, "figS4B_altej_tss_density_all_cohorts.pdf"),
  restrict_cohorts = NULL
)
