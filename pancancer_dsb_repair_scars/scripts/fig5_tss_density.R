#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(GenomicRanges)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- ifelse(length(args) >= 1, args[[1]], "config/config.yml")
cfg <- yaml::read_yaml(config_path)
out_dir <- cfg$paths$figure_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
proc <- cfg$paths$output_dir

read_hrd <- function(path) {
  h <- fread(path)
  setnames(h, "HRD-sum", "HRD_sum", skip_absent = TRUE)
  if (!("SampleID" %in% names(h))) setnames(h, "Sample", "SampleID", skip_absent = TRUE)
  h[, HRD_sum := as.numeric(HRD_sum)]
  h[, HRD_bin := cut(HRD_sum, breaks = c(-Inf, 21, 42, Inf), labels = c("HRD 0-21", "HRD 21-42", "HRD >42"))]
  h[, .(SampleID, HRD_sum, HRD_bin)]
}

read_tss_bed <- function(path) {
  t <- fread(path, header = FALSE)
  if (ncol(t) < 3) stop("TSS BED must contain at least chr, start, end.")
  GRanges(seqnames = t[[1]], ranges = IRanges(start = as.integer(t[[2]]) + 1L, end = as.integer(t[[3]])))
}

signed_dist_closest_end <- function(events, tss_gr) {
  left <- GRanges(events$Chromosome, IRanges(events$Start, events$Start))
  right <- GRanges(events$Chromosome, IRanges(events$End, events$End))
  dl <- distanceToNearest(left, tss_gr)
  dr <- distanceToNearest(right, tss_gr)
  out <- rep(NA_real_, nrow(events))
  if (length(dl) > 0) out[queryHits(dl)] <- mcols(dl)$distance
  if (length(dr) > 0) {
    i <- queryHits(dr)
    out[i] <- pmin(out[i], mcols(dr)$distance, na.rm = TRUE)
  }
  out / 1000
}

make_curve <- function(dt, tss_gr, window_kb = 20, bin_kb = 0.25) {
  bins <- seq(-window_kb, window_kb, by = bin_kb)
  rbindlist(lapply(levels(dt$HRD_bin), function(b) {
    sub <- dt[HRD_bin == b]
    if (nrow(sub) == 0) return(NULL)
    dist <- signed_dist_closest_end(sub, tss_gr)
    dist <- dist[!is.na(dist) & abs(dist) <= window_kb]
    h <- hist(dist, breaks = bins, plot = FALSE)
    data.table(HRD_bin = b, Dist_kb = h$mids, Density = h$counts / sum(h$counts))
  }), fill = TRUE)
}

plot_curve <- function(curves, title, outfile) {
  p <- ggplot(curves, aes(Dist_kb, Density, color = HRD_bin)) +
    geom_line(linewidth = 1.0, na.rm = TRUE) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
    theme_classic(base_size = 13) +
    labs(title = title, x = "Distance to nearest TSS (kb)", y = "Normalized event density", color = NULL)
  ggsave(file.path(out_dir, outfile), p, width = 7, height = 4.5)
}

hrd <- read_hrd(cfg$paths$hrd_file)
tss_gr <- read_tss_bed(cfg$paths$tss_bed)

ssa <- fread(file.path(proc, "PCAWG_SSA_like_deletions.tsv"))
ssa <- ssa[SSA_Candidate == TRUE]
setnames(ssa, "Sample", "SampleID", skip_absent = TRUE)
ssa <- merge(ssa, hrd, by = "SampleID")
plot_curve(make_curve(ssa, tss_gr), "SSA-like deletions relative to TSS", "fig5A_ssa_tss_density.pdf")

alt <- fread(file.path(proc, "PCAWG_AltEJ_like_deletions.tsv"))
alt <- alt[AltEJ_Candidate == TRUE]
setnames(alt, c("Tumor_Sample_Barcode", "Start_Position", "End_Position"), c("SampleID", "Start", "End"), skip_absent = TRUE)
alt <- merge(alt, hrd, by = "SampleID")
plot_curve(make_curve(alt, tss_gr), "Alt-EJ-like deletions relative to TSS", "fig5C_altej_tss_density.pdf")
