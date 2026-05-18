#!/usr/bin/env python

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

from common_figures import load_config, ensure_dir, read_summary, add_cohort_from_donor


HG19_CHR_LEN = {
    "chr1": 249250621, "chr2": 243199373, "chr3": 198022430, "chr4": 191154276,
    "chr5": 180915260, "chr6": 171115067, "chr7": 159138663, "chr8": 146364022,
    "chr9": 141213431, "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
    "chr13": 115169878, "chr14": 107349540, "chr15": 102531392, "chr16": 90354753,
    "chr17": 81195210, "chr18": 78077248, "chr19": 59128983, "chr20": 63025520,
    "chr21": 48129895, "chr22": 51304566,
}


def cfg_get(cfg, keys, default=None):
    x = cfg
    for k in keys:
        if not isinstance(x, dict) or k not in x:
            return default
        x = x[k]
    return x


def first_existing_col(df, candidates, label):
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(f"Missing {label}. Tried {candidates}. Found columns: {list(df.columns)}")


def standardize_sample_col(df):
    return df.rename(columns={
        "Sample": "SampleID",
        "Tumor_Sample_Barcode": "SampleID",
        "sample": "SampleID",
    })


def norm_chr(x):
    x = str(x).strip()
    if x in {"23", "X"}:
        return "chrX"
    if x in {"24", "Y"}:
        return "chrY"
    if x in {"MT", "M", "chrMT"}:
        return "chrM"
    return x if x.startswith("chr") else f"chr{x}"


def load_events(cfg, proc):
    altej = pd.read_csv(
        Path(cfg_get(cfg, ["paths", "altej_events"], proc / "PCAWG_AltEJ_like_deletions.tsv")),
        sep="\t",
        low_memory=False,
    )
    ssa = pd.read_csv(
        Path(cfg_get(cfg, ["paths", "ssa_events"], proc / "PCAWG_SSA_like_deletions.tsv")),
        sep="\t",
        low_memory=False,
    )

    altej = standardize_sample_col(altej)
    ssa = standardize_sample_col(ssa)

    if cfg_get(cfg, ["paths", "donor_file"]):
        altej = add_cohort_from_donor(altej, cfg["paths"]["donor_file"])
        ssa = add_cohort_from_donor(ssa, cfg["paths"]["donor_file"])

    for df in [altej, ssa]:
        if "Cohort" in df.columns:
            df["Cohort"] = df["Cohort"].astype(str).str.strip()

    return altej, ssa


def plot_hist_and_box(df, value_col, cohort_col, title_prefix, out_prefix, outdir, log_x=False):
    x = pd.to_numeric(df[value_col], errors="coerce")
    plot_df = df.assign(value=x).dropna(subset=["value"]).copy()
    plot_df = plot_df[plot_df["value"] >= 0]

    plt.figure(figsize=(5, 4))
    sns.histplot(plot_df["value"], bins=60, color="steelblue")
    if log_x:
        plt.xscale("log")
    plt.xlabel(value_col)
    plt.ylabel("Event count")
    plt.title(f"{title_prefix}: all tumors")
    plt.tight_layout()
    plt.savefig(outdir / f"{out_prefix}_all_hist.pdf")
    plt.close()

    if cohort_col in plot_df.columns:
        cohort_order = sorted(plot_df[cohort_col].dropna().unique())
        fig_w = max(10, 0.45 * len(cohort_order) + 4)

        plt.figure(figsize=(fig_w, 4.8))
        sns.boxplot(
            data=plot_df,
            x=cohort_col,
            y="value",
            order=cohort_order,
            showfliers=False,
            color="white",
        )
        sns.stripplot(
            data=plot_df,
            x=cohort_col,
            y="value",
            order=cohort_order,
            color="black",
            size=1.5,
            alpha=0.25,
        )
        if log_x:
            plt.yscale("log")
        plt.xlabel("")
        plt.ylabel(value_col)
        plt.title(f"{title_prefix}: by tumor type")
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig(outdir / f"{out_prefix}_by_cohort_box.pdf")
        plt.close()


def figure_s1(cfg, proc, outdir):
    altej, ssa = load_events(cfg, proc)

    altej_mh = first_existing_col(
        altej,
        ["MH_Length", "Microhomology_Length", "MH_len", "count_homosize", "Repeat_Length"],
        "Alt-EJ microhomology length",
    )
    altej_len = first_existing_col(
        altej,
        ["Del_Length", "Deletion_Size", "del_len"],
        "Alt-EJ deletion length",
    )
    ssa_homeology = first_existing_col(
        ssa,
        ["Repeat_Length", "Homeology_Length", "Homology_Length"],
        "SSA homeology length",
    )
    ssa_len = first_existing_col(
        ssa,
        ["Deletion_Size", "Del_Length", "del_len"],
        "SSA deletion length",
    )

    plot_hist_and_box(
        altej,
        altej_mh,
        "Cohort",
        "Alt-EJ microhomology length",
        "figS1A_altej_microhomology_length",
        outdir,
        log_x=False,
    )
    plot_hist_and_box(
        altej,
        altej_len,
        "Cohort",
        "Alt-EJ deletion length",
        "figS1B_altej_deletion_length",
        outdir,
        log_x=True,
    )
    plot_hist_and_box(
        ssa,
        ssa_homeology,
        "Cohort",
        "SSA homeology length",
        "figS1C_ssa_homeology_length",
        outdir,
        log_x=False,
    )
    plot_hist_and_box(
        ssa,
        ssa_len,
        "Cohort",
        "SSA deletion length",
        "figS1D_ssa_deletion_length",
        outdir,
        log_x=True,
    )


def classify_brca_status(mutations):
    mutations = mutations.copy()
    mutations.columns = [c.strip() for c in mutations.columns]

    sample_col = first_existing_col(mutations, ["SampleID", "Tumor_Sample_Barcode", "sample"], "sample")
    gene_col = first_existing_col(mutations, ["Gene", "Hugo_Symbol", "gene"], "gene")
    mut_col = first_existing_col(
        mutations,
        ["Consequence", "Variant_Classification", "effect", "mutation_type", "LoF"],
        "mutation/effect",
    )

    lof_terms = [
        "frameshift", "frame_shift", "stop_gained", "stopgain", "nonsense",
        "splice", "start_lost", "start_codon_loss", "stop codon gain",
    ]

    mutations[sample_col] = mutations[sample_col].astype(str).str.strip()
    mutations[gene_col] = mutations[gene_col].astype(str).str.upper()
    mutations[mut_col] = mutations[mut_col].astype(str).str.lower()

    mutations["is_lof"] = mutations[mut_col].apply(lambda s: any(t in s for t in lof_terms))
    brca = mutations[mutations[gene_col].isin(["BRCA1", "BRCA2"]) & mutations["is_lof"]]

    status = brca.groupby(sample_col)[gene_col].apply(lambda x: set(x)).reset_index()
    status["BRCA_status"] = status[gene_col].apply(
        lambda genes: "BRCA1 LoF" if "BRCA1" in genes else ("BRCA2 LoF" if "BRCA2" in genes else "WT")
    )

    return status[[sample_col, "BRCA_status"]].rename(columns={sample_col: "SampleID"})


def pairwise_welch_by_group(df, y, group_col, outpath):
    rows = []
    groups = [g for g in ["WT", "BRCA1 LoF", "BRCA2 LoF"] if g in set(df[group_col])]
    for a, b in itertools.combinations(groups, 2):
        xa = df.loc[df[group_col] == a, y].dropna()
        xb = df.loc[df[group_col] == b, y].dropna()
        if len(xa) >= 2 and len(xb) >= 2:
            _, p = ttest_ind(xa, xb, equal_var=False)
            rows.append({"comparison": f"{a} vs {b}", "p_value": p})
    out = pd.DataFrame(rows)
    if len(out):
        out["p_adj_fdr"] = multipletests(out["p_value"], method="fdr_bh")[1]
    out.to_csv(outpath, index=False)


def figure_s2(cfg, proc, outdir):
    mut_path = cfg_get(cfg, ["paths", "mutation_file"])
    if not mut_path:
        raise ValueError("Fig. S2 requires paths.mutation_file in config.yml.")

    mutations = pd.read_csv(mut_path, sep=cfg_get(cfg, ["paths", "mutation_sep"], "\t"), low_memory=False)
    brca_status = classify_brca_status(mutations)

    ssa = read_summary(Path(cfg_get(cfg, ["paths", "ssa_summary"], proc / "PCAWG_SSA_summary_by_sample.tsv")))
    alt = read_summary(Path(cfg_get(cfg, ["paths", "altej_summary"], proc / "PCAWG_AltEJ_summary_by_sample.tsv")))

    ssa = standardize_sample_col(ssa)[["SampleID", "Num_SSA"]]
    alt = standardize_sample_col(alt)[["SampleID", "Num_AltEJ"]]

    df = (
        ssa.merge(alt, on="SampleID", how="outer")
        .fillna({"Num_SSA": 0, "Num_AltEJ": 0})
        .merge(brca_status, on="SampleID", how="left")
    )
    df["BRCA_status"] = df["BRCA_status"].fillna("WT")
    order = ["WT", "BRCA1 LoF", "BRCA2 LoF"]

    for y, fname, ylabel in [
        ("Num_AltEJ", "figS2A_altej_brca_status.pdf", "Number of Alt-EJ-like deletions"),
        ("Num_SSA", "figS2B_ssa_brca_status.pdf", "Number of SSA-like deletions"),
    ]:
        plt.figure(figsize=(4.8, 4.2))
        sns.violinplot(data=df, x="BRCA_status", y=y, order=order, inner=None, cut=0, color="#9ecae1")
        sns.stripplot(data=df, x="BRCA_status", y=y, order=order, color="black", size=2, alpha=0.35)
        plt.xlabel("")
        plt.ylabel(ylabel)
        plt.tight_layout()
        plt.savefig(outdir / fname)
        plt.close()

        pairwise_welch_by_group(
            df,
            y,
            "BRCA_status",
            outdir / fname.replace(".pdf", "_welch_tests.csv"),
        )

    df.to_csv(outdir / "figS2_plot_data.csv", index=False)


def build_genome_axis(chr_len):
    chr_order = list(chr_len.keys())
    offsets = {}
    cur = 0
    for chrom in chr_order:
        offsets[chrom] = cur
        cur += chr_len[chrom]
    return chr_order, offsets, cur


def gaussian_smooth(counts, sigma_bins=2):
    radius = int(max(3, sigma_bins * 4))
    x = np.arange(-radius, radius + 1)
    kernel = np.exp(-(x**2) / (2 * sigma_bins**2))
    kernel = kernel / kernel.sum()
    return np.convolve(counts, kernel, mode="same")


def prepare_genome_events(df, donor_file=None):
    df = standardize_sample_col(df).rename(columns={
        "Start": "Start_Position",
        "End": "End_Position",
    })

    if donor_file and "Cohort" not in df.columns:
        df = add_cohort_from_donor(df, donor_file)

    start_col = first_existing_col(df, ["Start_Position", "Start"], "start")
    end_col = first_existing_col(df, ["End_Position", "End"], "end")

    df["Chromosome"] = df["Chromosome"].apply(norm_chr)
    df[start_col] = pd.to_numeric(df[start_col], errors="coerce")
    df[end_col] = pd.to_numeric(df[end_col], errors="coerce")
    df = df.dropna(subset=["Chromosome", start_col, end_col]).copy()
    df = df[df["Chromosome"].isin(HG19_CHR_LEN)].copy()

    chr_order, offsets, genome_len = build_genome_axis(HG19_CHR_LEN)
    df["x"] = df["Chromosome"].map(offsets) + ((df[start_col] + df[end_col]) / 2)
    df = df[(df["x"] >= 0) & (df["x"] < genome_len)].copy()

    return df, chr_order, offsets, genome_len


def plot_genome_density(df, chr_order, offsets, genome_len, cohorts, outpath, bin_size, sigma_bins, title):
    edges = np.arange(0, genome_len + bin_size, bin_size)
    centers = (edges[:-1] + edges[1:]) / 2

    rows = ["All"] + cohorts
    curves = {}
    all_vals = []

    for cohort in rows:
        sub = df if cohort == "All" else df[df["Cohort"] == cohort]
        counts, _ = np.histogram(sub["x"], bins=edges)
        counts_s = gaussian_smooth(counts, sigma_bins=sigma_bins)
        curves[cohort] = counts_s
        all_vals.append(counts_s)

    y_max = np.quantile(np.concatenate(all_vals), 0.995) * 1.10
    y_max = max(y_max, 1)

    fig, axes = plt.subplots(len(rows), 1, figsize=(16, 1.2 * len(rows) + 1.8), sharex=True)
    if len(rows) == 1:
        axes = [axes]

    for ax, cohort in zip(axes, rows):
        for chrom in chr_order:
            ax.axvline(offsets[chrom], color="#dddddd", lw=0.5, zorder=0)
        color = "black" if cohort == "All" else "#2b6cb0"
        ax.plot(centers, curves[cohort], color=color, lw=1.0)
        ax.text(0.005, 0.78, cohort, transform=ax.transAxes, fontsize=8)
        ax.set_ylim(0, y_max)
        ax.set_yticks([])
        ax.set_xticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

    fig.suptitle(title, y=0.995)
    fig.canvas.draw()

    bottom_ax = axes[-1]
    for chrom in chr_order:
        mid = offsets[chrom] + HG19_CHR_LEN[chrom] / 2
        x_disp = bottom_ax.transData.transform((mid, 0))[0]
        x_fig = fig.transFigure.inverted().transform((x_disp, 0))[0]
        fig.text(x_fig, 0.04, chrom.replace("chr", ""), ha="center", va="top", fontsize=7)

    fig.text(0.5, 0.012, "Genomic position (hg19, concatenated chromosomes)", ha="center", fontsize=10)
    plt.tight_layout(rect=[0, 0.06, 1, 0.97])
    plt.savefig(outpath)
    plt.close()


def figure_s3(cfg, proc, outdir):
    altej, ssa = load_events(cfg, proc)

    altej, chr_order, offsets, genome_len = prepare_genome_events(altej, cfg_get(cfg, ["paths", "donor_file"]))
    ssa, _, _, _ = prepare_genome_events(ssa, cfg_get(cfg, ["paths", "donor_file"]))

    cohorts = cfg_get(cfg, ["supplement", "figS3_cohorts"], None)
    if cohorts is None:
        cohorts = sorted(set(altej["Cohort"].dropna()).union(set(ssa["Cohort"].dropna())))

    plot_genome_density(
        altej,
        chr_order,
        offsets,
        genome_len,
        cohorts,
        outdir / "figS3A_altej_genomewide_all_cohorts.pdf",
        bin_size=cfg_get(cfg, ["supplement", "figS3_altej_bin"], 100_000),
        sigma_bins=cfg_get(cfg, ["supplement", "figS3_sigma_bins"], 2),
        title="Alt-EJ event density across genome",
    )

    plot_genome_density(
        ssa,
        chr_order,
        offsets,
        genome_len,
        cohorts,
        outdir / "figS3B_ssa_genomewide_all_cohorts.pdf",
        bin_size=cfg_get(cfg, ["supplement", "figS3_ssa_bin"], 2_000_000),
        sigma_bins=cfg_get(cfg, ["supplement", "figS3_sigma_bins"], 2),
        title="SSA event density across genome",
    )


def main(config):
    cfg = load_config(config)
    proc = Path(cfg["paths"]["output_dir"])
    outdir = ensure_dir(cfg_get(cfg, ["paths", "supplement_dir"], Path(cfg["paths"]["figure_dir"]) / "supplement"))

    figure_s1(cfg, proc, outdir)

    if cfg_get(cfg, ["paths", "mutation_file"]):
        figure_s2(cfg, proc, outdir)
    else:
        print("Skipping Fig. S2: no paths.mutation_file provided in config.yml")

    figure_s3(cfg, proc, outdir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/config.yml")
    args = parser.parse_args()
    main(args.config)
