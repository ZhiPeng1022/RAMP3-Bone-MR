# 组织特异性 MR：RAMP3 每组织 top eQTL 对 eBMD 的 Wald ratio

import gzip
import os

import pandas as pd
from scipy import stats


DATA = r"data/derived"
EXPOSURE = os.path.join(DATA, "ramp3_top_eqtl_exposure.txt")
OUTCOME = os.path.join(DATA, "ebmd_unpacked", "Biobank2-British-Bmd-As-C-Gwas-SumStats.txt.gz")


def load_outcome_subset(rsids):
    keep = set(rsids)
    df = pd.read_csv(
        OUTCOME,
        sep="\t",
        usecols=["RSID", "EA", "NEA", "BETA", "SE"],
        dtype={"RSID": str},
    )
    return df[df["RSID"].isin(keep)]


def main():
    exp = pd.read_csv(EXPOSURE, sep="\t")
    out = load_outcome_subset(exp["rs_id"].tolist())
    out["RSID"] = out["RSID"].astype(str)
    merged = exp.merge(out, left_on="rs_id", right_on="RSID")
    merged["BETA"] = merged["BETA"].astype(float)
    merged["SE"] = merged["SE"].astype(float)

    # 等位基因对齐（小写转大写）
    same = merged["alt"].str.upper() == merged["EA"].str.upper()
    flip = (merged["alt"].str.upper() == merged["NEA"].str.upper()) & (
        merged["ref"].str.upper() == merged["EA"].str.upper()
    )
    merged.loc[flip, "BETA"] = -merged.loc[flip, "BETA"]
    merged = merged[same | flip].copy()

    merged["beta"] = merged["beta"].astype(float)
    merged["se"] = merged["se"].astype(float)
    merged["SE"] = merged["SE"].astype(float)
    merged["wald_ratio"] = merged["BETA"] / merged["beta"]
    merged["wald_se"] = merged["SE"] / merged["beta"].abs()
    merged["p"] = 2 * stats.norm.sf(merged["wald_ratio"].abs() / merged["wald_se"])

    out_cols = ["tissue", "rs_id", "beta", "se", "out_beta", "out_se", "wald_ratio", "wald_se", "p"]
    result = merged.rename(columns={"BETA": "out_beta", "SE": "out_se"})[out_cols]
    print(result.to_string(index=False))
    result.to_csv(os.path.join(DATA, "ramp3_tissue_wise_ebmd.csv"), index=False)
    print("saved ->", os.path.join(DATA, "ramp3_tissue_wise_ebmd.csv"))


if __name__ == "__main__":
    main()
