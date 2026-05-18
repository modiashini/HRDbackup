# Core functions for detecting Alt-EJ-like and SSA-like deletion scars.

suppressPackageStartupMessages({
  library(data.table)
  library(BSgenome)
  library(Biostrings)
})

standardize_chromosome <- function(x, keep_chr_prefix = TRUE) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("23", "X")] <- "X"
  x[x %in% c("24", "Y")] <- "Y"
  x[x %in% c("MT", "M", "chrMT")] <- "M"

  if (keep_chr_prefix) {
    x[!grepl("^chr", x)] <- paste0("chr", x[!grepl("^chr", x)])
  } else {
    x <- sub("^chr", "", x)
  }

  x
}

standard_chromosomes <- paste0("chr", c(1:22, "X", "Y"))

standardize_variant_table <- function(x, sample_col = NULL) {
  dt <- as.data.table(x)
  setnames(dt, make.names(names(dt), unique = TRUE))

  rename_map <- c(
    Start_position = "Start_Position",
    start = "Start_Position",
    Start = "Start_Position",
    End_position = "End_Position",
    end = "End_Position",
    End = "End_Position",
    chr = "Chromosome",
    chromosome = "Chromosome",
    sample = "Tumor_Sample_Barcode",
    Sample = "Tumor_Sample_Barcode",
    SampleID = "Tumor_Sample_Barcode"
  )

  for (old in names(rename_map)) {
    new <- rename_map[[old]]
    if (old %in% names(dt) && !(new %in% names(dt))) {
      setnames(dt, old, new)
    }
  }

  if (!is.null(sample_col) && sample_col %in% names(dt)) {
    setnames(dt, sample_col, "Tumor_Sample_Barcode")
  }

  required <- c(
    "Chromosome", "Start_Position", "End_Position",
    "Tumor_Sample_Barcode", "Variant_Type", "Reference_Allele"
  )

  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop("Missing required variant columns: ", paste(missing, collapse = ", "))
  }

  dt[, Chromosome := standardize_chromosome(Chromosome)]
  dt[, Start_Position := as.integer(Start_Position)]
  dt[, End_Position := as.integer(End_Position)]
  dt[, Variant_Type := toupper(as.character(Variant_Type))]
  dt[, Tumor_Sample_Barcode := as.character(Tumor_Sample_Barcode)]

  dt <- dt[Chromosome %in% standard_chromosomes]
  dt <- dt[!is.na(Start_Position) & !is.na(End_Position)]
  dt[]
}

load_maf_rds <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path) else path
  x <- readRDS(con)
  standardize_variant_table(x)
}

microhomology_length <- function(deleted_seq, upstream, downstream) {
  deleted_seq <- toupper(as.character(deleted_seq))
  upstream <- toupper(as.character(upstream))
  downstream <- toupper(as.character(downstream))

  if (is.na(deleted_seq) || is.na(upstream) || is.na(downstream)) return(0L)

  n <- nchar(deleted_seq)
  if (is.na(n) || n < 2) return(0L)

  down_len <- 0L
  candidate <- substr(deleted_seq, 1, n - 1)

  while (nchar(candidate) > 0) {
    if (startsWith(downstream, candidate)) {
      down_len <- nchar(candidate)
      break
    }
    candidate <- substr(candidate, 1, nchar(candidate) - 1)
  }

  up_len <- 0L
  candidate <- substring(deleted_seq, 2)

  while (nchar(candidate) > 0) {
    if (endsWith(upstream, candidate)) {
      up_len <- nchar(candidate)
      break
    }
    candidate <- substring(candidate, 2)
  }

  as.integer(max(up_len, down_len))
}

annotate_indel_context <- function(maf_dt, ref_genome, flank_bp = 250) {
  dt <- standardize_variant_table(maf_dt)
  dt <- dt[Variant_Type %in% c("DEL", "INS")]

  if (nrow(dt) == 0) stop("No DEL/INS variants found.")

  dt[, Deleted_Seq := fifelse(
    Variant_Type == "DEL",
    as.character(Reference_Allele),
    NA_character_
  )]

  dt[, Deletion_Size := fifelse(
    Variant_Type == "DEL",
    pmax(
      nchar(as.character(Reference_Allele)),
      abs(End_Position - Start_Position + 1L),
      na.rm = TRUE
    ),
    NA_integer_
  )]

  dt[, upstream := as.character(BSgenome::getSeq(
    ref_genome,
    names = Chromosome,
    start = pmax(Start_Position - flank_bp, 1L),
    end = pmax(Start_Position - 1L, 1L)
  ))]

  dt[, downstream := as.character(BSgenome::getSeq(
    ref_genome,
    names = Chromosome,
    start = fifelse(Variant_Type == "INS", Start_Position, End_Position + 1L),
    end = fifelse(Variant_Type == "INS", Start_Position + flank_bp - 1L, End_Position + flank_bp)
  ))]

  dt[Variant_Type == "DEL",
     Microhomology_Length := mapply(microhomology_length, Deleted_Seq, upstream, downstream)]

  dt[Variant_Type != "DEL", Microhomology_Length := NA_integer_]

  dt[, Del_Length := Deletion_Size]
  dt[, MH_Length := Microhomology_Length]

  dt[]
}

classify_altej <- function(maf_dt,
                           ref_genome,
                           min_deletion_bp = 6,
                           min_microhomology_bp = 3,
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
  dt <- as.data.table(altej_dt)

  dt[, .(
    Num_DEL = sum(Variant_Type == "DEL", na.rm = TRUE),
    Num_AltEJ = sum(AltEJ_Candidate, na.rm = TRUE),
    AltEJ_Fraction = fifelse(
      sum(Variant_Type == "DEL", na.rm = TRUE) > 0,
      sum(AltEJ_Candidate, na.rm = TRUE) / sum(Variant_Type == "DEL", na.rm = TRUE),
      NA_real_
    ),
    Mean_MH_AltEJ = fifelse(
      sum(AltEJ_Candidate, na.rm = TRUE) > 0,
      mean(Microhomology_Length[AltEJ_Candidate == TRUE], na.rm = TRUE),
      NA_real_
    ),
    Mean_DelLen_AltEJ = fifelse(
      sum(AltEJ_Candidate, na.rm = TRUE) > 0,
      mean(Deletion_Size[AltEJ_Candidate == TRUE], na.rm = TRUE),
      NA_real_
    )
  ), by = .(SampleID = Tumor_Sample_Barcode)]
}

standardize_sv_table <- function(x, sample_col = NULL, effect_col = NULL) {
  dt <- as.data.table(x)
  setnames(dt, make.names(names(dt), unique = TRUE))

  rename_map <- c(
    chr = "Chromosome",
    chromosome = "Chromosome",
    start = "Start",
    end = "End",
    sample = "Sample",
    SampleID = "Sample"
  )

  for (old in names(rename_map)) {
    new <- rename_map[[old]]
    if (old %in% names(dt) && !(new %in% names(dt))) {
      setnames(dt, old, new)
    }
  }

  if (!is.null(sample_col) && sample_col %in% names(dt)) {
    setnames(dt, sample_col, "Sample")
  }

  if (!is.null(effect_col) && effect_col %in% names(dt) && effect_col != "effect") {
    setnames(dt, effect_col, "effect")
  }

  required <- c("Chromosome", "Start", "End", "Sample")
  missing <- setdiff(required, names(dt))

  if (length(missing) > 0) {
    stop("Missing required SV columns: ", paste(missing, collapse = ", "))
  }

  dt[, Chromosome := standardize_chromosome(Chromosome)]
  dt[, Start := as.integer(Start)]
  dt[, End := as.integer(End)]
  dt[, Sample := as.character(Sample)]

  dt <- dt[Chromosome %in% standard_chromosomes]
  dt <- dt[!is.na(Start) & !is.na(End)]
  dt[]
}

find_homeologous_repeat <- function(upstream,
                                    downstream,
                                    min_repeat_bp = 30,
                                    min_percent_identity = 80,
                                    gap_opening = -10,
                                    gap_extension = -0.5) {
  upstream <- DNAString(toupper(as.character(upstream)))
  downstream <- DNAString(toupper(as.character(downstream)))

  aln <- pairwiseAlignment(
    upstream,
    downstream,
    type = "local",
    gapOpening = gap_opening,
    gapExtension = gap_extension
  )

  aln_len <- width(alignedPattern(aln))
  pid_val <- Biostrings::pid(aln, type = "PID1")

  if (!is.na(aln_len) &&
      !is.na(pid_val) &&
      aln_len >= min_repeat_bp &&
      pid_val >= min_percent_identity) {
    list(
      length = as.integer(aln_len),
      seq = as.character(alignedPattern(aln)),
      percent_identity = as.numeric(pid_val)
    )
  } else {
    list(length = 0L, seq = "", percent_identity = 0)
  }
}

classify_ssa_from_sv <- function(sv_dt,
                                 ref_genome,
                                 min_deletion_bp = 50,
                                 min_repeat_bp = 30,
                                 min_percent_identity = 80,
                                 flank_bp = 500,
                                 effect_col = "effect",
                                 progress_every = 1000) {
  dt <- standardize_sv_table(sv_dt, effect_col = effect_col)

  if (effect_col %in% names(dt)) {
    dt <- dt[toupper(get(effect_col)) == "DEL"]
  }

  dt[, Deletion_Size := abs(End - Start) + 1L]
  dt <- dt[Deletion_Size >= min_deletion_bp]

  if (nrow(dt) == 0) stop("No SV deletions remain after filtering.")

  out <- vector("list", nrow(dt))

  for (i in seq_len(nrow(dt))) {
    row <- dt[i]

    upstream <- as.character(BSgenome::getSeq(
      ref_genome,
      names = row$Chromosome,
      start = max(1L, row$Start - flank_bp),
      end = row$Start - 1L
    ))

    downstream <- as.character(BSgenome::getSeq(
      ref_genome,
      names = row$Chromosome,
      start = row$End + 1L,
      end = row$End + flank_bp
    ))

    rep <- find_homeologous_repeat(
      upstream,
      downstream,
      min_repeat_bp = min_repeat_bp,
      min_percent_identity = min_percent_identity
    )

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

    if (!is.null(progress_every) && i %% progress_every == 0) {
      message("Processed ", i, " / ", nrow(dt), " SV deletions")
    }
  }

  res <- rbindlist(out, fill = TRUE)

  res[, SSA_Candidate := Repeat_Length >= min_repeat_bp &
       Repeat_PID >= min_percent_identity]

  res[]
}

summarize_ssa_by_sample <- function(ssa_dt) {
  dt <- as.data.table(ssa_dt)

  dt[, .(
    Num_Large_DEL = .N,
    Num_SSA = sum(SSA_Candidate, na.rm = TRUE),
    SSA_Fraction = sum(SSA_Candidate, na.rm = TRUE) / .N,
    Mean_Repeat_Length_SSA = fifelse(
      sum(SSA_Candidate, na.rm = TRUE) > 0,
      mean(Repeat_Length[SSA_Candidate == TRUE], na.rm = TRUE),
      NA_real_
    ),
    Mean_Deletion_Size_SSA = fifelse(
      sum(SSA_Candidate, na.rm = TRUE) > 0,
      mean(Deletion_Size[SSA_Candidate == TRUE], na.rm = TRUE),
      NA_real_
    )
  ), by = .(SampleID = Sample)]
}
