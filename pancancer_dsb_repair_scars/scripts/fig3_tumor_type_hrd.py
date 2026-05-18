#!/usr/bin/env python
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from common_figures import load_config, ensure_dir, read_hrd, read_summary, add_cohort_from_donor

def main(config):
    cfg = load_config(config)
    outdir = ensure_dir(cfg["paths"]["figure_dir"])
    proc = Path(cfg["paths"]["output_dir"])
    hrd = read_hrd(cfg["paths"]["hrd_file"])
    ssa = read_summary(proc / "PCAWG_SSA_summary_by_sample.tsv")[["SampleID", "Num_SSA"]]
    alt = read_summary(proc / "PCAWG_AltEJ_summary_by_sample.tsv")[["SampleID", "Num_AltEJ", "Num_DEL"]]
    d = ssa.merge(alt, on="SampleID", how="outer").fillna({"Num_SSA":0, "Num_AltEJ":0}).merge(hrd, on="SampleID")
    d = add_cohort_from_donor(d, cfg["paths"]["donor_file"])
    d = d.dropna(subset=["Cohort"])

    for y, fname in [("Num_SSA", "fig3A_ssa_tumor_type_hrd.pdf"), ("Num_AltEJ", "fig3B_altej_tumor_type_hrd.pdf")]:
        plt.figure(figsize=(10, 4.8))
        sns.violinplot(data=d, x="Cohort", y=y, hue="HRD_bin", cut=0, inner=None, scale="width")
        sns.stripplot(data=d, x="Cohort", y=y, hue="HRD_bin", dodge=True, color="black", size=1, alpha=0.35, legend=False)
        plt.xticks(rotation=90)
        plt.xlabel("")
        plt.ylabel(y.replace("Num_", "Number of "))
        plt.tight_layout()
        plt.savefig(outdir / fname)
        plt.close()

    # Fig 3C: normalized SSA vs Alt-EJ burden by total deletions.
    d["Prop_SSA"] = d["Num_SSA"] / d["Num_DEL"].replace(0, pd.NA)
    d["Prop_AltEJ"] = d["Num_AltEJ"] / d["Num_DEL"].replace(0, pd.NA)
    g = sns.FacetGrid(d.dropna(subset=["Prop_SSA", "Prop_AltEJ"]), col="HRD_bin", hue="Cohort", height=3.2, sharex=False, sharey=False)
    g.map_dataframe(sns.scatterplot, x="Prop_AltEJ", y="Prop_SSA", s=12, alpha=0.75)
    g.set_axis_labels("Alt-EJ-like / total deletions", "SSA-like / total deletions")
    g.tight_layout()
    g.savefig(outdir / "fig3C_normalized_joint_burden.pdf")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--config", default="config/config.yml")
    main(p.parse_args().config)
