# Pan-cancer genomic scars of Alt-EJ and SSA

Code accompanying *“Pan-Cancer Genomic Scars of Alternative End Joining and Single-Strand Annealing”* by Ashini Modi, Alessandro Zito, and Giovanni Parmigiani.

This repository identifies genomic deletion scars associated with two mutagenic double-strand break (DSB) repair pathways:

- **Alternative end joining (Alt-EJ/MMEJ)** — characterized by deletions flanked by short microhomology tracts;
- **Single-strand annealing (SSA)** — characterized by larger deletions flanked by long homologous repeats.

The code reproduces the main analyses and figures from the manuscript, including:

- pan-cancer quantification of Alt-EJ- and SSA-like deletions,
- relationships between repair-pathway usage and homologous recombination deficiency (HRD),
- tumor-type-specific repair patterns,
- genome-wide hotspot analyses,
- and transcription-associated enrichment near transcription start sites (TSSs).

---

# Repository structure

```text
R/
  repair_scar_detection.R

scripts/
  00_run_detection.R
  common_figures.py
  fig2_hrd_scaling.py
  fig3_tumor_type_hrd.py
  fig4_genomewide_density.py
  fig5_tss_density.R

config/
  config.yml
