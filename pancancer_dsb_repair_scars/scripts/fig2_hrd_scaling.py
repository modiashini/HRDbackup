#!/usr/bin/env python
import argparse
from pathlib import Path
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from common_figures import load_config, ensure_dir, read_hrd, read_summary

def pairwise_welch(df, y):
    rows = []
    for a, b in itertools.combinations(["0–21", "21–42", "42+"], 2):
        xa = df.loc[df.HRD_bin == a, y].dropna()
        xb = df.loc[df.HRD_bin == b, y].dropna()
        if len(xa) >= 2 and len(xb) >= 2:
            stat, p = ttest_ind(xa, xb, equal_var=False)
            rows.append({"comparison": f"{a} vs {b}", "p_value": p})
    out = pd.DataFrame(rows)
    if len(out):
        out["p_adj_fdr"] = multipletests(out.p_value, method="fdr_bh")[1]
    return out

def main(config):
    cfg = load_config(config)
    outdir = ensure_dir(cfg["paths"]["figure_dir"])
    proc = Path(cfg["paths"]["output_dir"])
    hrd = read_hrd(cfg["paths"]["hrd_file"])
    ssa = read_summary(proc / "PCAWG_SSA_summary_by_sample.tsv")[["SampleID", "Num_SSA"]]
    alt = read_summary(proc / "PCAWG_AltEJ_summary_by_sample.tsv")[["SampleID", "Num_AltEJ"]]
    d = ssa.merge(alt, on="SampleID", how="outer").fillna({"Num_SSA":0, "Num_AltEJ":0}).merge(hrd, on="SampleID")

    # Fig 2B/C: burden by HRD bin
    for y, fname in [("Num_SSA", "fig2B_ssa_by_hrd.pdf"), ("Num_AltEJ", "fig2C_altej_by_hrd.pdf")]:
        plt.figure(figsize=(4.5, 4))
        sns.stripplot(data=d, x="HRD_bin", y=y, color="black", alpha=0.35, size=2)
        sns.boxplot(data=d, x="HRD_bin", y=y, showfliers=False, color="white")
        plt.xlabel("HRD bin")
        plt.ylabel(y.replace("Num_", "Number of "))
        plt.tight_layout()
        plt.savefig(outdir / fname)
        plt.close()
        pairwise_welch(d, y).to_csv(outdir / fname.replace(".pdf", "_welch_tests.csv"), index=False)

    # Fig 2A needs cohort means; uses Cohort if present in Alt-EJ event table.
    events = pd.read_csv(proc / "PCAWG_AltEJ_like_deletions.tsv", sep="\t", usecols=lambda c: c in ["Tumor_Sample_Barcode","Cohort"], low_memory=False)
    events = events.rename(columns={"Tumor_Sample_Barcode":"SampleID"}).drop_duplicates()
    means = d.merge(events, on="SampleID", how="left").dropna(subset=["Cohort"])
    means = means.groupby("Cohort", as_index=False).agg(Mean_SSA=("Num_SSA","mean"), Mean_AltEJ=("Num_AltEJ","mean"), Mean_HRD=("HRD_sum","mean"))
    plt.figure(figsize=(5.2, 4.5))
    sc = plt.scatter(means.Mean_SSA, means.Mean_AltEJ, c=means.Mean_HRD, s=40)
    for _, r in means.iterrows():
        plt.text(r.Mean_SSA, r.Mean_AltEJ, str(r.Cohort), fontsize=6)
    plt.xlabel("Mean SSA-like deletions per sample")
    plt.ylabel("Mean Alt-EJ-like deletions per sample")
    plt.colorbar(sc, label="Mean HRD score")
    plt.tight_layout()
    plt.savefig(outdir / "fig2A_cohort_means.pdf")
    plt.close()

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--config", default="config/config.yml")
    main(p.parse_args().config)
