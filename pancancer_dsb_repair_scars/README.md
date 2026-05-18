Pan-cancer genomic scars of Alt-EJ and SSA

Code accompanying “Pan-Cancer Genomic Scars of Alternative End Joining and Single-Strand Annealing” by Ashini Modi, Alessandro Zito, and Giovanni Parmigiani.

This repository identifies genomic deletion scars associated with two mutagenic double-strand break (DSB) repair pathways:

Alternative end joining (Alt-EJ/MMEJ) — characterized by deletions flanked by short microhomology tracts
Single-strand annealing (SSA) — characterized by larger deletions flanked by long homologous repeats

The code reproduces the main analyses and figures from the manuscript, including:

pan-cancer quantification of Alt-EJ- and SSA-like deletions
relationships between repair-pathway usage and homologous recombination deficiency (HRD)
tumor-type-specific repair patterns
genome-wide hotspot analyses
transcription-associated enrichment near transcription start sites (TSSs)
Repository structure
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
Overview

The main detection pipeline classifies deletion events from PCAWG/ICGC mutation data into repair-associated genomic scar categories.

Alt-EJ-like deletions

Alt-EJ-like deletions are identified from indel calls using:

deletion length
junction microhomology length

Following the manuscript definition, deletions are classified as Alt-EJ-like if they:

are larger than 5 bp
contain 3–25 bp of microhomology at the breakpoint junction
SSA-like deletions

SSA-like deletions are identified from structural variant deletion calls by searching for homologous repeat tracts flanking deletion breakpoints.

Events are classified as SSA-like if they contain:

repeat tracts ≥30 bp
≥80% sequence identity
Main analyses reproduced
Figure 2 — HRD scaling of SSA and Alt-EJ

The repository reproduces analyses showing that:

Alt-EJ-like deletions increase strongly with HRD score
SSA-like deletions also increase with HRD, though more modestly
specific deletion-size and homology regimes correlate most strongly with HR deficiency
Figure 3 — Tumor-type-specific repair pathway usage

The code reproduces comparisons of SSA and Alt-EJ burden across tumor types and HRD bins, including:

elevated SSA usage in prostate and liver tumors despite low HRD
enrichment of Alt-EJ scars in melanoma, kidney, and liver cancers
Figure 4 — Genome-wide deletion landscapes

Genome-wide density analyses identify:

recurrent Alt-EJ hotspots in lymphoid malignancies at the IGH locus
recurrent SSA enrichments in repeat-rich chromosomal regions
Figure 5 — Transcription-associated enrichment

The repository reproduces TSS-distance analyses demonstrating:

strong enrichment of SSA-like deletions near transcription start sites in HR-proficient lymphoid tumors
moderate TSS enrichment for Alt-EJ-like deletions
Running the pipeline
1. Configure paths

Edit:

config/config.yml

to point to:

PCAWG/ICGC mutation files
structural variant calls
HRD annotations
output directories
2. Run deletion classification
Rscript scripts/00_run_detection.R config/config.yml

This generates processed event-level and sample-level summary tables for:

Alt-EJ-like deletions
SSA-like deletions
combined repair scar burden
3. Reproduce manuscript figures
python scripts/fig2_hrd_scaling.py --config config/config.yml
python scripts/fig3_tumor_type_hrd.py --config config/config.yml
python scripts/fig4_genomewide_density.py --config config/config.yml
Rscript scripts/fig5_tss_density.R config/config.yml
Expected outputs

The pipeline produces processed tables such as:

PCAWG_AltEJ_like_deletions.tsv
PCAWG_AltEJ_summary_by_sample.tsv
PCAWG_SSA_like_deletions.tsv
PCAWG_SSA_summary_by_sample.tsv

along with manuscript figure panels written to the figures/ directory.

Data source

Analyses were performed on whole-genome-sequenced tumors from the International Cancer Genome Consortium (ICGC) / Pan-Cancer Analysis of Whole Genomes (PCAWG) project aligned to the hg19/GRCh37 reference genome.
