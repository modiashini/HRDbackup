#!/usr/bin/env python
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from common_figures import load_config, ensure_dir, add_cohort_from_donor, chromosome_offsets, HG19_CHR_LEN, CHR_ORDER, standardize_chromosome

def smooth_counts(counts, sigma_bins=3):
    radius = int(4 * sigma_bins)
    x = np.arange(-radius, radius + 1)
    kernel = np.exp(-(x ** 2) / (2 * sigma_bins ** 2))
    kernel /= kernel.sum()
    return np.convolve(counts, kernel, mode="same")

def prepare_events(path, kind):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if kind == "SSA":
        df = df.rename(columns={"Sample": "SampleID", "Start": "Start_Position", "End": "End_Position"})
        if "SSA_Candidate" in df.columns:
            df = df[df["SSA_Candidate"].astype(str).str.lower().isin(["true", "1", "yes", "t"])]
    else:
        df = df.rename(columns={"Tumor_Sample_Barcode": "SampleID"})
        if "AltEJ_Candidate" in df.columns:
            df = df[df["AltEJ_Candidate"].astype(str).str.lower().isin(["true", "1", "yes", "t"])]
    df["Chromosome"] = standardize_chromosome(df["Chromosome"])
    df = df[df["Chromosome"].isin(CHR_ORDER)].copy()
    df["Start_Position"] = pd.to_numeric(df["Start_Position"], errors="coerce")
    df["End_Position"] = pd.to_numeric(df["End_Position"], errors="coerce")
    df["Midpoint"] = (df["Start_Position"] + df["End_Position"]) / 2
    offsets = chromosome_offsets()
    df["Genome_Pos"] = df.apply(lambda r: offsets[r.Chromosome] + r.Midpoint, axis=1)
    return df.dropna(subset=["Genome_Pos"])

def plot_density(df, cohorts, bin_bp, title, outfile):
    total_len = sum(HG19_CHR_LEN[c] for c in CHR_ORDER)
    bins = np.arange(0, total_len + bin_bp, bin_bp)
    centers = (bins[:-1] + bins[1:]) / 2
    fig, axes = plt.subplots(len(cohorts), 1, figsize=(10, 1.15 * len(cohorts)), sharex=True)
    if len(cohorts) == 1:
        axes = [axes]
    for ax, cohort in zip(axes, cohorts):
        sub = df if cohort == "All" else df[df["Cohort"] == cohort]
        counts, _ = np.histogram(sub["Genome_Pos"], bins=bins)
        ax.plot(centers, smooth_counts(counts), lw=1.0)
        ax.text(0.005, 0.72, cohort, transform=ax.transAxes, ha="left", va="center", fontsize=8)
        ax.set_yticks([])
        for chrom in CHR_ORDER:
            ax.axvline(chromosome_offsets()[chrom], color="lightgray", lw=0.4, alpha=0.7)
        ax.spines[["top", "right", "left"]].set_visible(False)
    axes[-1].set_xticks([chromosome_offsets()[c] + HG19_CHR_LEN[c]/2 for c in CHR_ORDER])
    axes[-1].set_xticklabels([c.replace("chr", "") for c in CHR_ORDER], fontsize=7)
    fig.suptitle(title, y=0.995)
    fig.supxlabel("Genomic position (hg19 concatenated chromosomes)")
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

def main(config):
    cfg = load_config(config)
    outdir = ensure_dir(cfg["paths"]["figure_dir"])
    proc = Path(cfg["paths"]["output_dir"])
    alt = prepare_events(proc / "PCAWG_AltEJ_like_deletions.tsv", "AltEJ")
    ssa = prepare_events(proc / "PCAWG_SSA_like_deletions.tsv", "SSA")
    alt = add_cohort_from_donor(alt, cfg["paths"]["donor_file"])
    ssa = add_cohort_from_donor(ssa, cfg["paths"]["donor_file"])
    plot_density(alt, ["All", "Liver-HCC", "Lymph-BNHL", "Lymph-CLL"], 100_000,
                 "Alt-EJ-like deletion density across genome", outdir / "fig4A_altej_genomewide_density.pdf")
    plot_density(ssa, ["All", "Breast-AdenoCA", "Ovary-AdenoCA", "Eso-AdenoCA"], 2_000_000,
                 "SSA-like deletion density across genome", outdir / "fig4B_ssa_genomewide_density.pdf")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--config", default="config/config.yml")
    main(p.parse_args().config)
