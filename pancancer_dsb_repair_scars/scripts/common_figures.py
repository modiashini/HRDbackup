from pathlib import Path
import yaml
import numpy as np
import pandas as pd

HRD_BINS = [0, 21, 42, np.inf]
HRD_LABELS = ["0–21", "21–42", "42+"]
CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
HG19_CHR_LEN = {
    "chr1":249250621,"chr2":243199373,"chr3":198022430,"chr4":191154276,"chr5":180915260,
    "chr6":171115067,"chr7":159138663,"chr8":146364022,"chr9":141213431,"chr10":135534747,
    "chr11":135006516,"chr12":133851895,"chr13":115169878,"chr14":107349540,"chr15":102531392,
    "chr16":90354753,"chr17":81195210,"chr18":78077248,"chr19":59128983,"chr20":63025520,
    "chr21":48129895,"chr22":51304566,"chrX":155270560,"chrY":59373566,
}

def load_config(path):
    with open(path, "r") as f:
        cfg = yaml.safe_load(f)
    return cfg

def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)
    return Path(path)

def read_hrd(path):
    hrd = pd.read_csv(path, low_memory=False)
    if "HRD-sum" in hrd.columns:
        hrd = hrd.rename(columns={"HRD-sum": "HRD_sum"})
    if "Sample" in hrd.columns and "SampleID" not in hrd.columns:
        hrd = hrd.rename(columns={"Sample": "SampleID"})
    hrd["SampleID"] = hrd["SampleID"].astype(str)
    hrd["HRD_sum"] = pd.to_numeric(hrd["HRD_sum"], errors="coerce")
    hrd["HRD_bin"] = pd.cut(hrd["HRD_sum"], HRD_BINS, labels=HRD_LABELS, include_lowest=True)
    return hrd[["SampleID", "HRD_sum", "HRD_bin"]].dropna(subset=["HRD_sum"])

def read_summary(path, sample_col_candidates=("SampleID", "Sample", "Tumor_Sample_Barcode")):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    for col in sample_col_candidates:
        if col in df.columns:
            df = df.rename(columns={col: "SampleID"})
            break
    if "SampleID" not in df.columns:
        raise ValueError(f"No sample id column found in {path}")
    df["SampleID"] = df["SampleID"].astype(str)
    return df

def add_cohort_from_donor(df, donor_file):
    donor = pd.read_csv(donor_file, sep="\t", low_memory=False)
    candidate_sample_cols = ["icgc_sample_id", "sample", "SampleID", "specimen_id", "donor_unique_id"]
    sample_col = next((c for c in candidate_sample_cols if c in donor.columns), None)
    cohort_col = next((c for c in ["histology_abbreviation", "cohort", "dcc_project_code", "project_code"] if c in donor.columns), None)
    if sample_col is None or cohort_col is None:
        raise ValueError("Could not infer sample/cohort columns from donor file.")
    donor = donor.rename(columns={sample_col: "SampleID", cohort_col: "Cohort"})[["SampleID", "Cohort"]].drop_duplicates()
    donor["SampleID"] = donor["SampleID"].astype(str)
    return df.merge(donor, on="SampleID", how="left")

def chromosome_offsets():
    offsets, cur = {}, 0
    for chrom in CHR_ORDER:
        offsets[chrom] = cur
        cur += HG19_CHR_LEN[chrom]
    return offsets

def standardize_chromosome(series):
    s = series.astype(str).str.strip()
    return np.where(s.str.startswith("chr"), s, "chr" + s)
