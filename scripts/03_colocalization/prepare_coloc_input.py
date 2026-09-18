# 生成 R coloc 输入文件：RAMP3 全血 eQTL vs eBMD 共享 SNP

import os

import pandas as pd
from scipy import stats


DATA = r"data/derived"
EQTL = os.path.join(DATA, "gtex_ramp3_eqtl.csv")
OUTCOME = os.path.join(DATA, "ebmd_unpacked", "Biobank2-British-Bmd-As-C-Gwas-SumStats.txt.gz")
OUT = os.path.join(DATA, "coloc_input_ramp3_ebmd.csv")


def main():
    eq = pd.read_csv(EQTL)
    eq = eq[eq["tissue"] == "Whole_Blood"].dropna(subset=["p", "nes"]).copy()
    eq["z"] = abs(stats.norm.ppf(eq["p"] / 2.0))
    eq["se1"] = abs(eq["nes"]) / eq["z"]
    eq["beta1"] = eq["nes"]
    eq["ref"] = eq["variant_id"].str.extract(r"_([ACGT])_([ACGT])_b38$")[0]
    eq["alt"] = eq["variant_id"].str.extract(r"_([ACGT])_([ACGT])_b38$")[1]

    keep = set(eq["snp"])
    out = pd.read_csv(
        OUTCOME,
        sep="\t",
        usecols=["RSID", "EA", "NEA", "EAF", "BETA", "SE"],
        dtype={"RSID": str},
    )
    out = out[out["RSID"].isin(keep)].copy()
    merged = eq.merge(out, left_on="snp", right_on="RSID")

    same = merged["alt"].str.upper() == merged["EA"].str.upper()
    flip = (merged["alt"].str.upper() == merged["NEA"].str.upper()) & (
        merged["ref"].str.upper() == merged["EA"].str.upper()
    )
    merged.loc[flip, "BETA"] = -merged.loc[flip, "BETA"]
    merged = merged[same | flip].copy()

    result = merged[["snp", "beta1", "se1", "BETA", "SE", "EAF"]].rename(
        columns={"BETA": "beta2", "SE": "se2", "EAF": "maf"}
    )
    result["maf"] = result["maf"].astype(float).clip(0.01, 0.99)
    result.to_csv(OUT, index=False)
    print("saved", len(result), "SNPs ->", OUT)


if __name__ == "__main__":
    main()

