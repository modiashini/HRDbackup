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

from common_figures import load_config, ensure_dir, read_hrd, read_summary


# ----------------------------
# General helpers
# ----------------------------

def pairwise_welch(df, y):
    rows = []
    for a, b in itertools.combinations(["0–21", "21–42", "42+"], 2):
        xa = df.loc[df["HRD_bin"] == a, y].dropna()
        xb = df.loc[df["HRD_bin"] == b, y].dropna()

        if len(xa) >= 2 and len(xb) >= 2:
            _, p = ttest_ind(xa, xb, equal_var=False)
            rows.append({"comparison": f"{a} vs {b}", "p_value": p})

    out = pd.DataFrame(rows)
    if len(out):
        out["p_adj_fdr"] = multipletests(out["p_value"], method="fdr_bh")[1]
    return out


def first_existing_col(df, candidates, table_name):
    for col in candidates:
        if col in df.columns:
            return col

    raise ValueError(
        f"Could not find one of {candidates} in {table_name}. "
        f"Available columns include: {list(df.columns[:80])}"
    )


def standardize_sample_column(df):
    rename_map = {
        "Tumor_Sample_Barcode": "SampleID",
        "Sample": "SampleID",
        "SampleID": "SampleID",
    }
    for old, new in rename_map.items():
        if old in df.columns:
            return df.rename(columns={old: new})
    raise ValueError("Could not find a sample column.")


def save_facet_regplot(
    plot_df,
    outpath,
    ylabel,
    title=None,
    height=2.45,
    aspect=1.05,
):
    sns.set_style("white")
    sns.set_context("paper", font_scale=1.05)

    g = sns.FacetGrid(
        plot_df,
        col="Bin",
        col_wrap=4,
        height=height,
        aspect=aspect,
        sharex=True,
        sharey=True,
    )

    g.map_dataframe(
        sns.regplot,
        x="HRD_sum",
        y="prop",
        lowess=False,
        ci=95,
        scatter_kws={"s": 13, "alpha": 0.32},
        line_kws={"lw": 1.8, "color": "black"},
    )

    g.set_titles("{col_name}", size=10)
    g.set_axis_labels("HRD Sum", ylabel)

    for ax in g.axes.flatten():
        ax.grid(False)
        ax.set_xlim(0, 100)
        ax.set_ylim(-0.03, 1.03)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=8)

    if title:
        g.fig.suptitle(title, y=1.03, fontsize=12)

    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()

    plot_df.to_csv(str(outpath).replace(".pdf", "_plot_data.csv"), index=False)


def make_proportion_df(
    events,
    hrd,
    value_col,
    bin_edges,
    bin_labels,
    min_events_per_sample=5,
):
    events = events.copy()
    events = standardize_sample_column(events)

    events[value_col] = pd.to_numeric(events[value_col], errors="coerce")
    events = events.dropna(subset=["SampleID", value_col]).copy()
    events = events[events[value_col] >= 0].copy()

    events["Bin"] = pd.cut(
        events[value_col],
        bins=bin_edges,
        labels=bin_labels,
        include_lowest=True,
        right=True,
    )

    events = events.dropna(subset=["Bin"]).copy()
    events["Bin"] = pd.Categorical(events["Bin"], categories=bin_labels, ordered=True)

    counts = (
        events.groupby(["SampleID", "Bin"], observed=True)
        .size()
        .rename("n_bin")
        .reset_index()
    )

    totals = (
        events.groupby("SampleID", observed=True)
        .size()
        .rename("n_total")
        .reset_index()
    )

    totals = totals[totals["n_total"] >= min_events_per_sample].copy()

    long = counts.merge(totals, on="SampleID", how="inner")
    long["prop"] = long["n_bin"] / long["n_total"]

    grid = pd.MultiIndex.from_product(
        [totals["SampleID"].unique(), bin_labels],
        names=["SampleID", "Bin"],
    ).to_frame(index=False)

    long = grid.merge(
        long[["SampleID", "Bin", "prop"]],
        on=["SampleID", "Bin"],
        how="left",
    )

    long["prop"] = long["prop"].fillna(0)
    long["Bin"] = pd.Categorical(long["Bin"], categories=bin_labels, ordered=True)

    plot_df = long.merge(hrd[["SampleID", "HRD_sum"]], on="SampleID", how="inner")
    return plot_df


# ----------------------------
# Figure 2
# ----------------------------

def main(config):
    cfg = load_config(config)

    outdir = ensure_dir(cfg["paths"]["figure_dir"])
    proc = Path(cfg["paths"]["output_dir"])

    hrd = read_hrd(cfg["paths"]["hrd_file"])
    hrd["HRD_sum"] = pd.to_numeric(hrd["HRD_sum"], errors="coerce")
    hrd = hrd.dropna(subset=["HRD_sum"]).copy()
    hrd = hrd[hrd["HRD_sum"] <= 100].copy()

    ssa_summary = read_summary(proc / "PCAWG_SSA_summary_by_sample.tsv")
    alt_summary = read_summary(proc / "PCAWG_AltEJ_summary_by_sample.tsv")

    ssa_summary = ssa_summary[["SampleID", "Num_SSA"]]
    alt_summary = alt_summary[["SampleID", "Num_AltEJ"]]

    d = (
        ssa_summary
        .merge(alt_summary, on="SampleID", how="outer")
        .fillna({"Num_SSA": 0, "Num_AltEJ": 0})
        .merge(hrd, on="SampleID", how="inner")
    )

    # ----------------------------
    # Fig. 2B/C: event burden by HRD bin
    # ----------------------------

    for y, fname in [
        ("Num_SSA", "fig2B_ssa_by_hrd.pdf"),
        ("Num_AltEJ", "fig2C_altej_by_hrd.pdf"),
    ]:
        plt.figure(figsize=(4.5, 4))

        sns.stripplot(
            data=d,
            x="HRD_bin",
            y=y,
            color="black",
            alpha=0.35,
            size=2,
        )

        sns.boxplot(
            data=d,
            x="HRD_bin",
            y=y,
            showfliers=False,
            color="white",
        )

        plt.xlabel("HRD bin")
        plt.ylabel(y.replace("Num_", "Number of "))
        plt.tight_layout()
        plt.savefig(outdir / fname)
        plt.close()

        pairwise_welch(d, y).to_csv(
            outdir / fname.replace(".pdf", "_welch_tests.csv"),
            index=False,
        )

    # ----------------------------
    # Fig. 2A: cohort-level mean SSA vs Alt-EJ
    # ----------------------------

    cohort_events = pd.read_csv(
        proc / "PCAWG_AltEJ_like_deletions.tsv",
        sep="\t",
        usecols=lambda c: c in ["Tumor_Sample_Barcode", "SampleID", "Cohort"],
        low_memory=False,
    )

    cohort_events = standardize_sample_column(cohort_events).drop_duplicates()

    means = (
        d.merge(cohort_events, on="SampleID", how="left")
        .dropna(subset=["Cohort"])
        .groupby("Cohort", as_index=False)
        .agg(
            Mean_SSA=("Num_SSA", "mean"),
            Mean_AltEJ=("Num_AltEJ", "mean"),
            Mean_HRD=("HRD_sum", "mean"),
        )
    )

    plt.figure(figsize=(5.2, 4.5))
    sc = plt.scatter(
        means["Mean_SSA"],
        means["Mean_AltEJ"],
        c=means["Mean_HRD"],
        s=40,
    )

    for _, r in means.iterrows():
        plt.text(
            r["Mean_SSA"],
            r["Mean_AltEJ"],
            str(r["Cohort"]),
            fontsize=6,
        )

    sns.regplot(
        data=means,
        x="Mean_SSA",
        y="Mean_AltEJ",
        scatter=False,
        color="black",
        line_kws={"lw": 1.5},
    )

    plt.xlabel("Mean SSA-like deletions per sample")
    plt.ylabel("Mean Alt-EJ-like deletions per sample")
    plt.colorbar(sc, label="Mean HRD score")
    plt.tight_layout()
    plt.savefig(outdir / "fig2A_cohort_means.pdf")
    plt.close()

    means.to_csv(outdir / "fig2A_cohort_means_data.csv", index=False)

    # ----------------------------
    # Load event-level tables for Fig. 2D/E
    # ----------------------------

    altej_events = pd.read_csv(
        proc / "PCAWG_AltEJ_like_deletions.tsv",
        sep="\t",
        low_memory=False,
    )
    altej_events = standardize_sample_column(altej_events)

    ssa_events = pd.read_csv(
        proc / "PCAWG_SSA_like_deletions.tsv",
        sep="\t",
        low_memory=False,
    )
    ssa_events = standardize_sample_column(ssa_events)

    altej_del_col = first_existing_col(
        altej_events,
        ["Del_Length", "Deletion_Size", "del_len"],
        "PCAWG_AltEJ_like_deletions.tsv",
    )

    altej_mh_col = first_existing_col(
        altej_events,
        ["MH_Length", "Microhomology_Length", "MH_len", "count_homosize", "Repeat_Length"],
        "PCAWG_AltEJ_like_deletions.tsv",
    )

    ssa_del_col = first_existing_col(
        ssa_events,
        ["Deletion_Size", "Del_Length", "del_len"],
        "PCAWG_SSA_like_deletions.tsv",
    )

    ssa_repeat_col = first_existing_col(
        ssa_events,
        ["Repeat_Length", "Homeology_Length", "Homology_Length"],
        "PCAWG_SSA_like_deletions.tsv",
    )

    # ----------------------------
    # Fig. 2D1: deletion size among deletions <100 bp
    # ----------------------------

    small_del = altej_events.copy()
    small_del[altej_del_col] = pd.to_numeric(small_del[altej_del_col], errors="coerce")
    small_del = small_del[
        (small_del[altej_del_col] > 0) &
        (small_del[altej_del_col] < 100)
    ].copy()

    fig2d1 = make_proportion_df(
        events=small_del,
        hrd=hrd,
        value_col=altej_del_col,
        bin_edges=[1, 5, 20, 40, np.inf],
        bin_labels=["1–5 bp", "6–20 bp", "20–40 bp", "40+ bp"],
        min_events_per_sample=5,
    )

    save_facet_regplot(
        fig2d1,
        outdir / "fig2D1_deletion_size_lt100bp.pdf",
        ylabel="Proportion of deletions <100 bp",
        title="Deletion size",
    )

    # ----------------------------
    # Fig. 2D2: microhomology size among deletions <100 bp
    # ----------------------------

    fig2d2 = make_proportion_df(
        events=small_del,
        hrd=hrd,
        value_col=altej_mh_col,
        bin_edges=[-0.1, 1, 5, 10, np.inf],
        bin_labels=["0–1 bp", "2–5 bp", "6–10 bp", "11+ bp"],
        min_events_per_sample=5,
    )

    save_facet_regplot(
        fig2d2,
        outdir / "fig2D2_microhomology_size_lt100bp.pdf",
        ylabel="Proportion of deletions <100 bp",
        title="Microhomology size",
    )

    # ----------------------------
    # Fig. 2E1: deletion size among deletions >100 bp
    # ----------------------------

    large_del = ssa_events.copy()
    large_del[ssa_del_col] = pd.to_numeric(large_del[ssa_del_col], errors="coerce")
    large_del = large_del[large_del[ssa_del_col] > 100].copy()

    fig2e1 = make_proportion_df(
        events=large_del,
        hrd=hrd,
        value_col=ssa_del_col,
        bin_edges=[100, 1000, 10000, 100000, np.inf],
        bin_labels=["100–1k bp", "1k–10k bp", "10k–100k bp", "100k+ bp"],
        min_events_per_sample=5,
    )

    save_facet_regplot(
        fig2e1,
        outdir / "fig2E1_deletion_size_gt100bp.pdf",
        ylabel="Proportion of deletions >100 bp",
        title="Deletion size",
    )

    # ----------------------------
    # Fig. 2E2: homeology size among deletions >100 bp
    # ----------------------------

    fig2e2 = make_proportion_df(
        events=large_del,
        hrd=hrd,
        value_col=ssa_repeat_col,
        bin_edges=[-0.1, 0, 30, 100, np.inf],
        bin_labels=["0 bp", "1–30 bp", "31–100 bp", "100+ bp"],
        min_events_per_sample=5,
    )

    save_facet_regplot(
        fig2e2,
        outdir / "fig2E2_homeology_size_gt100bp.pdf",
        ylabel="Proportion of deletions >100 bp",
        title="Homeology (>80%) size",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/config.yml")
    args = parser.parse_args()
    main(args.config)
