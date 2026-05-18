# Pan-cancer genomic scars of Alt-EJ and SSA

Cleaned, reusable code for detecting alternative end joining (Alt-EJ/MMEJ) and single-strand annealing (SSA) deletion scars from PCAWG/ICGC-style mutation files, then reproducing manuscript figures.

## Repository layout

```text
R/
  repair_scar_detection.R      # core reusable Alt-EJ + SSA detection functions
scripts/
  00_run_detection.R           # runs the main detection pipeline
  common_figures.py            # shared plotting/data-loading helpers
  fig2_hrd_scaling.py          # Figure 2 panels: HRD scaling + feature bins
  fig3_tumor_type_hrd.py       # Figure 3 panels: tumor type + HRD bins
  fig4_genomewide_density.py   # Figure 4 panels: genome-wide density tracks
  fig5_tss_density.R           # Figure 5 panels: TSS-distance density curves
config/config.yml              # all input/output paths and thresholds
```

## What the main detection code does

`R/repair_scar_detection.R` is the main reusable code. It defines everything needed to:

1. standardize ICGC/PCAWG MAF-like indel tables;
2. calculate deletion length and uncapped microhomology length;
3. classify Alt-EJ-like deletions using configurable thresholds;
4. scan structural-variant deletions for flanking homeologous repeats;
5. classify SSA-like deletions using configurable repeat length and percent identity thresholds;
6. write both event-level and sample-level summary tables.

The defaults match the manuscript logic: Alt-EJ-like deletions are deletion events >5 bp with 2–25 bp microhomology, while SSA-like deletions are larger deletions flanked by repeats ≥30 bp with ≥80% identity.

## Quick start

Install R dependencies:

```r
install.packages(c("data.table", "dplyr", "yaml", "BiocManager"))
BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg19", "Biostrings", "GenomicRanges"))
```

Install Python dependencies:

```bash
pip install pandas numpy matplotlib seaborn scipy statsmodels pyyaml
```

Edit `config/config.yml`, then run:

```bash
Rscript scripts/00_run_detection.R config/config.yml
python scripts/fig2_hrd_scaling.py --config config/config.yml
python scripts/fig3_tumor_type_hrd.py --config config/config.yml
python scripts/fig4_genomewide_density.py --config config/config.yml
Rscript scripts/fig5_tss_density.R config/config.yml
```

## Expected processed outputs

The detection pipeline writes:

```text
data/processed/PCAWG_AltEJ_like_deletions.tsv
data/processed/PCAWG_AltEJ_summary_by_sample.tsv
data/processed/PCAWG_SSA_like_deletions.tsv
data/processed/PCAWG_SSA_summary_by_sample.tsv
```

Figure scripts read these processed tables and write PDFs/PNGs to `figures/`.

## Notes

This cleaned version removes duplicated code blocks, dummy placeholder functions, notebook/Colab-specific paths, and overwritten variables. Each script loads its own dependencies and reads paths from one config file, so figures can be rerun independently.
