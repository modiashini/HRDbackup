# Core Alt-EJ and SSA scar detection functions.
# Source this file from scripts/00_run_detection.R or from an interactive R session.

suppressPackageStartupMessages({
  library(data.table)
  library(BSgenome)
  library(Biostrings)
})

standardize_variant_table <- function(x) {
  dt <- as.data.table(x)
  setnames(dt, make.names(names(dt), unique = TRUE))
  rename_map <- c(
    Start_position = "Start_Position", start = "Start_Position", Start = "Start_Position",
    End_position = "End_Position", end = "End_Position", End = "End_Position",
    chr = "Chromosome", chromosome = "Chromosome",
    sample = "Tumor_Sample_Barcode", Sample = "Tumor_Sample_Barcode", SampleID = "Tumor_Sample_Barcode"
  )
  for (old in names(rename_map)) {
    if (old %in% names(dt) && !(rename_map[[old]] %in% names(dt))) setnames(dt, old, rename_map[[old]])
  }
  required <- c("Chromosome", "Start_Position", "End_Position", "Tumor_Sample_Barcode")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))
  dt[, Chromosome := as.character(Chromosome)]
  dt[!grepl("^chr", Chromosome), Chromosome := paste0("chr", Chromosome)]
  dt[, Start_Position := as.integer(Start_Position)]
  dt[, End_Position := as.integer(End_Position)]
  dt[]
}

load_maf_rds <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path) else path
  standardize_variant_table(readRDS(con))
}

microhomology_length <- function(ref_allele, upstream, downstream) {
  ref_allele <- toupper(as.character(ref_allele))
  upstream <- toupper(as.character(upstream))
  downstream <- toupper(as.character(downstream))
  n <- nchar(ref_allele)
  if (is.na(n) || n < 2) return(0L)

  down_len <- 0L
  candidate <- substr(ref_allele, 1, n - 1)
  while (nchar(candidate) > 0) {
    if (startsWith(downstream, candidate)) { down_len <- nchar(candidate); break }
    candidate <- substr(candidate, 1, nchar(candidate) - 1)
  }

  up_len <- 0L
  candidate <- substring(ref_allele, 2)
  while (nchar(candidate) > 0) {
    if (endsWith(upstream, candidate)) { up_len <- nchar(candidate); break }
    candidate <- substring(candidate, 2)
  }
  as.integer(max(up_len, down_len))
}

annotate_indel_context <- function(maf_dt, ref_genome, flank_bp = 250) {
  dt <- standardize_variant_table(maf_dt)
  if (!("Variant_Type" %in% names(dt))) stop("MAF table must contain Variant_Type.")
  dt <- dt[toupper(Variant_Type) %in% c("DEL", "INS")]
  if (nrow(dt) == 0) stop("No DEL/INS variants found.")

  dt[, Variant_Type := toupper(Variant_Type)]
  if (!("Reference_Allele" %in% names(dt))) stop("MAF table must contain Reference_Allele.")
  dt[, Deletion_Size := fifelse(Variant_Type == "DEL", nchar(as.character(Reference_Allele)), NA_integer_)]

  dt[, upstream := as.character(BSgenome::getSeq(
    ref_genome, names = Chromosome,
    start = pmax(Start_Position - flank_bp, 1L),
    end = pmax(Start_Position - 1L, 1L)
  ))]
  dt[, downstream := as.character(BSgenome::getSeq(
    ref_genome, names = Chromosome,
    start = fifelse(Variant_Type == "INS", Start_Position, End_Position + 1L),
    end = fifelse(Variant_Type == "INS", Start_Position + flank_bp - 1L, End_Position + flank_bp)
  ))]
  dt[Variant_Type == "DEL", Microhomology_Length := mapply(microhomology_length, Reference_Allele, upstream, downstream)]
  dt[Variant_Type != "DEL", Microhomology_Length := NA_integer_]
  dt[]
}

classify_altej <- function(maf_dt, ref_genome,
                           min_deletion_bp = 6,
                           min_microhomology_bp = 2,
                           max_microhomology_bp = 25,
                           flank_bp = 250) {
  dt <- annotate_indel_context(maf_dt, ref_genome, flank_bp = flank_bp)
  dt[, AltEJ_Candidate := Variant_Type == "DEL" &
       Deletion_Size >= min_deletion_bp &
       Microhomology_Length >= min_microhomology_bp &
       Microhomology_Length <= max_microhomology_bp]
  dt[]
}

summarize_altej_by_sample <- function(altej_dt) {
  as.data.table(altej_dt)[, .(
    Num_DEL = sum(Variant_Type == "DEL", na.rm = TRUE),
    Num_AltEJ = sum(AltEJ_Candidate, na.rm = TRUE),
    AltEJ_Fraction = fifelse(sum(Variant_Type == "DEL", na.rm = TRUE) > 0,
                              sum(AltEJ_Candidate, na.rm = TRUE) / sum(Variant_Type == "DEL", na.rm = TRUE),
                              NA_real_)
  ), by = .(SampleID = Tumor_Sample_Barcode)]
}

standardize_sv_table <- function(x) {
  dt <- as.data.table(x)
  setnames(dt, make.names(names(dt), unique = TRUE))
  rename_map <- c(chr = "Chromosome", chromosome = "Chromosome", start = "Start", end = "End", sample = "Sample", SampleID = "Sample")
  for (old in names(rename_map)) {
    if (old %in% names(dt) && !(rename_map[[old]] %in% names(dt))) setnames(dt, old, rename_map[[old]])
  }
  required <- c("Chromosome", "Start", "End", "Sample")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) stop("Missing required SV columns: ", paste(missing, collapse = ", "))
  dt[, Chromosome := as.character(Chromosome)]
  dt[!grepl("^chr", Chromosome), Chromosome := paste0("chr", Chromosome)]
  dt[, Start := as.integer(Start)]
  dt[, End := as.integer(End)]
  dt[]
}

find_homeologous_repeat <- function(upstream, downstream, min_repeat_bp = 30, min_pid = 80) {
  upstream <- DNAString(toupper(as.character(upstream)))
  downstream <- DNAString(toupper(as.character(downstream)))
  aln <- pairwiseAlignment(upstream, downstream, type = "local", gapOpening = -10, gapExtension = -0.5)
  aln_len <- width(alignedPattern(aln))
  pid_val <- Biostrings::pid(aln, type = "PID1")
  if (!is.na(aln_len) && !is.na(pid_val) && aln_len >= min_repeat_bp && pid_val >= min_pid) {
    list(length = as.integer(aln_len), seq = as.character(alignedPattern(aln)), percent_identity = as.numeric(pid_val))
  } else {
    list(length = 0L, seq = "", percent_identity = 0)
  }
}

classify_ssa_from_sv <- function(sv_dt, ref_genome,
                                 min_deletion_bp = 50,
                                 min_repeat_bp = 30,
                                 min_percent_identity = 80,
                                 flank_bp = 500,
                                 effect_col = "effect") {
  dt <- standardize_sv_table(sv_dt)
  if (effect_col %in% names(dt)) dt <- dt[toupper(get(effect_col)) == "DEL"]
  dt <- dt[!is.na(Start) & !is.na(End)]
  dt[, Deletion_Size := abs(End - Start) + 1L]
  dt <- dt[Deletion_Size >= min_deletion_bp]
  if (nrow(dt) == 0) stop("No SV deletions remain after filtering.")

  out <- vector("list", nrow(dt))
  for (i in seq_len(nrow(dt))) {
    row <- dt[i]
    upstream <- as.character(BSgenome::getSeq(ref_genome, row$Chromosome, max(1L, row$Start - flank_bp), row$Start - 1L))
    downstream <- as.character(BSgenome::getSeq(ref_genome, row$Chromosome, row$End + 1L, row$End + flank_bp))
    rep <- find_homeologous_repeat(upstream, downstream, min_repeat_bp, min_percent_identity)
    out[[i]] <- data.table(
      Sample = row$Sample,
      Chromosome = row$Chromosome,
      Start = row$Start,
      End = row$End,
      Deletion_Size = row$Deletion_Size,
      Repeat_Length = rep$length,
      Repeat_Seq = rep$seq,
      Repeat_PID = rep$percent_identity
    )
    if (i %% 1000 == 0) message("Processed ", i, " / ", nrow(dt), " SV deletions")
  }
  res <- rbindlist(out, fill = TRUE)
  res[, SSA_Candidate := Repeat_Length >= min_repeat_bp & Repeat_PID >= min_percent_identity]
  res[]
}

summarize_ssa_by_sample <- function(ssa_dt) {
  as.data.table(ssa_dt)[, .(
    Num_Large_DEL = .N,
    Num_SSA = sum(SSA_Candidate, na.rm = TRUE),
    SSA_Fraction = sum(SSA_Candidate, na.rm = TRUE) / .N
  ), by = .(SampleID = Sample)]
}
