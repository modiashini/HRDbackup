#!/usr/bin/env python

import argparse
from pathlib import Path
import itertools

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import gaussian_kde, ttest_ind
from scipy.ndimage import label
from statsmodels.stats.multitest import multipletests
from matplotlib.lines import Line2D

from common_figures import load_config, ensure_dir, read_hrd, read_summary, add_cohort_from_donor


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
    raise ValueError(f"Missing {label}. Tried {candidates}. Found: {list(df.columns)}")


def standardize_cohort_names(df, rename_map=None):
    df = df.copy()
    df["Cohort"] = df["Cohort"].astype(str).str.strip()
    if rename_map:
        df["Cohort"] = df["Cohort"].replace(rename_map)
    return df


def prepare_hrd_bins(df, bins, labels):
    df = df.copy()
    df["HRD_sum"] = pd.to_numeric(df["HRD_sum"], errors="coerce")
    df = df.dropna(subset=["HRD_sum"]).copy()
    df["HRD_sum_capped"] = df["HRD_sum"].clip(upper=bins[-1])
    df["HRD_bin"] = pd.cut(
        df["HRD_sum_capped"],
        bins=bins,
        labels=labels,
        include_lowest=True,
        right=True,
    )
    df = df.dropna(subset=["HRD_bin"]).copy()
    df["HRD_bin"] = pd.Categorical(df["HRD_bin"], categories=labels, ordered=True)
    return df


def load_fig3_data(cfg):
    proc = Path(cfg["paths"]["output_dir"])

    hrd = read_hrd(cfg["paths"]["hrd_file"])

    ssa_path = Path(cfg_get(cfg, ["paths", "ssa_summary"], proc / "PCAWG_SSA_summary_by_sample.tsv"))
    alt_path = Path(cfg_get(cfg, ["paths", "altej_summary"], proc / "PCAWG_AltEJ_summary_by_sample.tsv"))

    ssa = read_summary(ssa_path)
    alt = read_summary(alt_path)

    ssa_num_col = first_existing_col(ssa, ["Num_SSA", "n_SSA", "SSA_count"], "SSA count column")
    alt_num_col = first_existing_col(alt, ["Num_AltEJ", "n_AltEJ", "AltEJ_count"], "Alt-EJ count column")

    ssa = ssa[["SampleID", ssa_num_col]].rename(columns={ssa_num_col: "Num_SSA"})
    alt_keep = ["SampleID", alt_num_col]
    for optional_col in ["Num_DEL", "Num_Small_DEL", "Num_Large_DEL"]:
        if optional_col in alt.columns:
            alt_keep.append(optional_col)

    alt = alt[alt_keep].rename(columns={alt_num_col: "Num_AltEJ"})

    df = (
        ssa.merge(alt, on="SampleID", how="outer")
        .fillna({"Num_SSA": 0, "Num_AltEJ": 0})
        .merge(hrd, on="SampleID", how="inner")
    )

    donor_file = cfg_get(cfg, ["paths", "donor_file"])
    if donor_file:
        df = add_cohort_from_donor(df, donor_file)

    if "Cohort" not in df.columns:
        raise ValueError("Figure 3 requires a Cohort column or paths.donor_file in config.yml.")

    cohort_rename = cfg_get(cfg, ["plotting", "cohort_rename"], {})
    df = standardize_cohort_names(df, cohort_rename)

    keep_cohorts = cfg_get(cfg, ["plotting", "keep_cohorts"], None)
    if keep_cohorts:
        df = df[df["Cohort"].isin(keep_cohorts)].copy()

    hrd_bins = cfg_get(cfg, ["analysis", "hrd_bins"], [0, 21, 42, 500])
    hrd_labels = cfg_get(cfg, ["analysis", "hrd_bin_labels"], ["0–21", "21–42", "42+"])
    df = prepare_hrd_bins(df, hrd_bins, hrd_labels)

    cohort_order = cfg_get(cfg, ["plotting", "cohort_order"], None)
    if cohort_order:
        cohort_order = [c for c in cohort_order if c in set(df["Cohort"])]
    else:
        cohort_order = sorted(df["Cohort"].dropna().unique())

    df["Cohort"] = pd.Categorical(df["Cohort"], categories=cohort_order, ordered=True)

    df["Num_SSA"] = pd.to_numeric(df["Num_SSA"], errors="coerce").fillna(0)
    df["Num_AltEJ"] = pd.to_numeric(df["Num_AltEJ"], errors="coerce").fillna(0)

    return df, cohort_order, hrd_labels


def plot_grouped_violin(df, cohort_order, hrd_labels, y_col, ylabel, outpath, cfg):
    sns.set_style("white")

    palette = cfg_get(cfg, ["plotting", "hrd_palette"], None)
    if palette is None:
        palette = dict(zip(hrd_labels, sns.color_palette("Blues", n_colors=len(hrd_labels))))

    fig_width = cfg_get(cfg, ["plotting", "fig3_ab_width"], max(12, 0.55 * len(cohort_order) + 6))
    fig_height = cfg_get(cfg, ["plotting", "fig3_ab_height"], 6)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    sns.violinplot(
        data=df,
        x="Cohort",
        y=y_col,
        hue="HRD_bin",
        order=cohort_order,
        hue_order=hrd_labels,
        palette=palette,
        dodge=True,
        cut=0,
        inner=None,
        linewidth=1,
        scale="width",
        ax=ax,
    )

    rng = np.random.default_rng(cfg_get(cfg, ["plotting", "jitter_seed"], 0))
    violin_width = 0.8
    n_hue = len(hrd_labels)

    for xi, cohort in enumerate(cohort_order):
        for hi, hrd_bin in enumerate(hrd_labels):
            sub = df[(df["Cohort"] == cohort) & (df["HRD_bin"] == hrd_bin)]
            if sub.empty:
                continue

            step = violin_width / n_hue
            x_center = xi + (hi - (n_hue - 1) / 2) * step

            x_j = rng.uniform(-0.06, 0.06, size=len(sub))
            y_j = rng.uniform(0, 0.18, size=len(sub))

            ax.scatter(
                x_center + x_j,
                sub[y_col].to_numpy() + y_j,
                s=cfg_get(cfg, ["plotting", "point_size"], 5),
                alpha=cfg_get(cfg, ["plotting", "point_alpha"], 0.55),
                color="black",
                zorder=3,
            )

    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    ax.grid(False)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles[:len(hrd_labels)],
        labels[:len(hrd_labels)],
        title="HRD bin",
        frameon=False,
        loc="upper left",
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


def pairwise_tests_hrd_low(df, hrd_low_label, y_col, outpath, min_n=2):
    sub_df = df[df["HRD_bin"] == hrd_low_label].copy()

    valid_counts = sub_df.groupby("Cohort", observed=True)["SampleID"].nunique()
    valid_cohorts = valid_counts[valid_counts >= min_n].index.tolist()

    rows = []
    for c1, c2 in itertools.combinations(sorted(valid_cohorts), 2):
        x1 = sub_df.loc[sub_df["Cohort"] == c1, y_col].dropna().to_numpy()
        x2 = sub_df.loc[sub_df["Cohort"] == c2, y_col].dropna().to_numpy()

        stat, p = ttest_ind(x1, x2, equal_var=False, nan_policy="omit")

        rows.append({
            "Cohort_1": c1,
            "Cohort_2": c2,
            "N_1": len(x1),
            "N_2": len(x2),
            "Mean_1": np.mean(x1),
            "Mean_2": np.mean(x2),
            "Median_1": np.median(x1),
            "Median_2": np.median(x2),
            "t_stat": stat,
            "p_value": p,
        })

    out = pd.DataFrame(rows)
    if len(out):
        reject, p_adj, _, _ = multipletests(out["p_value"], method="fdr_bh")
        out["p_adj_fdr"] = p_adj
        out["significant_fdr_0.05"] = reject
        out["Mean_diff"] = out["Mean_1"] - out["Mean_2"]
        out["Abs_mean_diff"] = out["Mean_diff"].abs()
        out = out.sort_values(["p_adj_fdr", "p_value", "Abs_mean_diff"])

    out.to_csv(outpath, index=False)


def add_total_deletions(df, cfg):
    proc = Path(cfg["paths"]["output_dir"])
    composite_path = Path(cfg_get(
        cfg,
        ["paths", "composite_summary"],
        proc / "PCAWG_composite_SSA_AltEJ_SV_only.csv",
    ))

    df = df.copy()

    if composite_path.exists():
        comp = pd.read_csv(composite_path, low_memory=False)

        small_col = first_existing_col(comp, ["Num_Small_DEL", "Small_DEL", "Num_small_deletions"], "small deletion count")
        large_col = first_existing_col(comp, ["Num_Large_DEL", "Large_DEL", "Num_large_deletions"], "large deletion count")

        comp = comp[["SampleID", small_col, large_col]].copy()
        comp[small_col] = pd.to_numeric(comp[small_col], errors="coerce").fillna(0)
        comp[large_col] = pd.to_numeric(comp[large_col], errors="coerce").fillna(0)
        comp["Total_dels"] = comp[small_col] + comp[large_col]

        df = df.merge(comp[["SampleID", "Total_dels"]], on="SampleID", how="left")

    elif "Num_DEL" in df.columns:
        df["Total_dels"] = pd.to_numeric(df["Num_DEL"], errors="coerce")

    else:
        raise FileNotFoundError(
            "Fig. 3C needs either paths.composite_summary with Num_Small_DEL and Num_Large_DEL, "
            "or Num_DEL in the Alt-EJ summary file."
        )

    df.loc[df["Total_dels"] <= 0, "Total_dels"] = np.nan
    df["Prop_SSA"] = df["Num_SSA"] / df["Total_dels"]
    df["Prop_AltEJ"] = df["Num_AltEJ"] / df["Total_dels"]

    return df.dropna(subset=["Prop_SSA", "Prop_AltEJ"]).copy()


def assign_outside_core(sub_df, xcol, ycol, core_mass=0.99, bw_method=0.25):
    out = sub_df.copy()

    if len(out) < 3:
        out["density"] = np.nan
        out["InCore"] = True
        out["OutsideCore"] = False
        return out

    x = out[xcol].to_numpy()
    y = out[ycol].to_numpy()
    xy = np.vstack([x, y])

    if np.unique(x).size < 2 or np.unique(y).size < 2:
        xy = xy + 1e-6 * np.random.default_rng(0).normal(size=xy.shape)

    kde = gaussian_kde(xy, bw_method=bw_method)
    dens = kde(xy)
    out["density"] = dens

    dens_sorted = np.sort(dens)[::-1]
    cum_mass = np.cumsum(dens_sorted) / np.sum(dens_sorted)
    idx = min(np.searchsorted(cum_mass, core_mass), len(dens_sorted) - 1)
    cutoff = dens_sorted[idx]

    xpad = max(0.02, 0.08 * (x.max() - x.min() if x.max() > x.min() else 1))
    ypad = max(0.002, 0.08 * (y.max() - y.min() if y.max() > y.min() else 1))

    xmin, xmax = x.min() - xpad, x.max() + xpad
    ymin, ymax = max(-0.001, y.min() - ypad), y.max() + ypad

    xx, yy = np.meshgrid(
        np.linspace(xmin, xmax, 300),
        np.linspace(ymin, ymax, 300),
    )

    zz = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
    mask = zz >= cutoff
    labeled, _ = label(mask)

    peak_idx = np.unravel_index(np.argmax(zz), zz.shape)
    main_component = labeled[peak_idx]

    if main_component == 0:
        out["InCore"] = out["density"] >= cutoff
        out["OutsideCore"] = ~out["InCore"]
        return out

    x_idx = np.clip(np.searchsorted(xx[0], x) - 1, 0, xx.shape[1] - 1)
    y_idx = np.clip(np.searchsorted(yy[:, 0], y) - 1, 0, yy.shape[0] - 1)

    out["InCore"] = labeled[y_idx, x_idx] == main_component
    out["OutsideCore"] = ~out["InCore"]

    return out


def prepare_fig3c(df, hrd_labels, cfg):
    core_mass = cfg_get(cfg, ["plotting", "fig3c_core_mass"], 0.99)
    bw_method = cfg_get(cfg, ["plotting", "fig3c_bw_method"], 0.25)

    panel_dfs = []
    for hrd_bin in hrd_labels:
        sub = df[df["HRD_bin"] == hrd_bin].copy()
        if len(sub):
            panel_dfs.append(
                assign_outside_core(
                    sub,
                    xcol="Prop_AltEJ",
                    ycol="Prop_SSA",
                    core_mass=core_mass,
                    bw_method=bw_method,
                )
            )

    return pd.concat(panel_dfs, ignore_index=True)


def plot_fig3c(df, hrd_labels, outpath, cfg):
    sns.set_style("whitegrid")

    cohort_palette = cfg_get(cfg, ["plotting", "cohort_palette"], {})
    outside_cohorts = sorted(df.loc[df["OutsideCore"], "Cohort"].dropna().unique())

    palette = {}
    for c in outside_cohorts:
        if c in cohort_palette:
            palette[c] = cohort_palette[c]

    missing = [c for c in outside_cohorts if c not in palette]
    auto_cols = sns.color_palette("tab20", n_colors=max(len(missing), 1))
    for c, col in zip(missing, auto_cols):
        palette[c] = col

    fig_w = cfg_get(cfg, ["plotting", "fig3c_width"], 6)
    fig_h = cfg_get(cfg, ["plotting", "fig3c_height"], 13)
    xlim = cfg_get(cfg, ["plotting", "fig3c_xlim"], [-0.01, 0.40])
    ylim = cfg_get(cfg, ["plotting", "fig3c_ylim"], [-0.001, 0.030])

    fig, axes = plt.subplots(len(hrd_labels), 1, figsize=(fig_w, fig_h), sharex=True, sharey=True)
    if len(hrd_labels) == 1:
        axes = [axes]

    for ax, hrd_bin in zip(axes, hrd_labels):
        sub = df[df["HRD_bin"] == hrd_bin].copy()

        if len(sub) > 2:
            sns.kdeplot(
                data=sub,
                x="Prop_AltEJ",
                y="Prop_SSA",
                fill=True,
                levels=30,
                thresh=0.02,
                color="lightgray",
                ax=ax,
            )

        sns.scatterplot(
            data=sub,
            x="Prop_AltEJ",
            y="Prop_SSA",
            color="lightgray",
            s=16,
            alpha=0.45,
            legend=False,
            ax=ax,
        )

        sub_out = sub[sub["OutsideCore"] & sub["Cohort"].isin(palette)].copy()

        if len(sub_out):
            sns.scatterplot(
                data=sub_out,
                x="Prop_AltEJ",
                y="Prop_SSA",
                hue="Cohort",
                palette=palette,
                s=34,
                alpha=0.95,
                edgecolor="black",
                linewidth=0.25,
                legend=False,
                ax=ax,
            )

        ax.set_title(f"HRD bin: {hrd_bin}", fontsize=10)
        ax.set_xlabel("number of AltEJ/total deletions")
        ax.set_ylabel("number of SSA/total deletions")
        ax.grid(False)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    legend_handles = [
        Line2D(
            [0], [0],
            marker="o",
            color="w",
            label=cohort,
            markerfacecolor=color,
            markeredgecolor="black",
            markersize=7,
        )
        for cohort, color in palette.items()
    ]

    fig.legend(
        handles=legend_handles,
        title="Tumor type",
        loc="center right",
        bbox_to_anchor=(1.48, 0.5),
        frameon=False,
    )

    plt.tight_layout(rect=[0, 0, 0.88, 1])
    plt.savefig(outpath, bbox_inches="tight")
    plt.close()


def main(config):
    cfg = load_config(config)
    outdir = ensure_dir(cfg["paths"]["figure_dir"])

    df, cohort_order, hrd_labels = load_fig3_data(cfg)

    plot_grouped_violin(
        df,
        cohort_order,
        hrd_labels,
        y_col="Num_SSA",
        ylabel="Number of SSA-like deletions",
        outpath=outdir / "fig3A_ssa_tumor_type_hrd.pdf",
        cfg=cfg,
    )

    plot_grouped_violin(
        df,
        cohort_order,
        hrd_labels,
        y_col="Num_AltEJ",
        ylabel="Number of Alt-EJ-like deletions",
        outpath=outdir / "fig3B_altej_tumor_type_hrd.pdf",
        cfg=cfg,
    )

    pairwise_tests_hrd_low(
        df,
        hrd_low_label=hrd_labels[0],
        y_col="Num_SSA",
        outpath=outdir / "fig3A_ssa_low_hrd_pairwise_welch_tests.csv",
        min_n=cfg_get(cfg, ["analysis", "pairwise_min_n"], 2),
    )

    pairwise_tests_hrd_low(
        df,
        hrd_low_label=hrd_labels[0],
        y_col="Num_AltEJ",
        outpath=outdir / "fig3B_altej_low_hrd_pairwise_welch_tests.csv",
        min_n=cfg_get(cfg, ["analysis", "pairwise_min_n"], 2),
    )

    df_c = add_total_deletions(df, cfg)
    df_c = prepare_fig3c(df_c, hrd_labels, cfg)

    df.to_csv(outdir / "fig3AB_plot_data.csv", index=False)
    df_c.to_csv(outdir / "fig3C_plot_data.csv", index=False)

    plot_fig3c(
        df_c,
        hrd_labels,
        outpath=outdir / "fig3C_normalized_joint_burden_by_hrd.pdf",
        cfg=cfg,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/config.yml")
    args = parser.parse_args()
    main(args.config)
