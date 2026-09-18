# 简化版 coloc（single causal variant ABF）：RAMP3 全血 eQTL vs eBMD
# 说明：初步分析，正式投稿用 R coloc 包复核

import math
import os
import re

import pandas as pd
from scipy import stats


DATA = r"data/derived"
EQTL = os.path.join(DATA, "gtex_ramp3_eqtl.csv")
OUTCOME = os.path.join(DATA, "ebmd_unpacked", "Biobank2-British-Bmd-As-C-Gwas-SumStats.txt.gz")
P1 = P2 = P12 = 1e-4
W = 0.2 ** 2


def logsumexp(xs):
    m = max(xs)
    return m + math.log(sum(math.exp(v - m) for v in xs))


def logabf(z, se):
    r = W / se ** 2
    return 0.5 * math.log(r / (r + 1)) + (z ** 2 / 2) * (1 / (r + 1))


def main():
    eq = pd.read_csv(EQTL)
    eq = eq[eq["tissue"] == "Whole_Blood"].dropna(subset=["p", "nes"]).copy()
    eq["z"] = eq["nes"] / (abs(eq["nes"]) / abs(stats.norm.ppf(eq["p"] / 2.0)))
    eq["se"] = abs(eq["nes"]) / abs(stats.norm.ppf(eq["p"] / 2.0))
    eq["ref"] = eq["variant_id"].str.extract(r"_([ACGT])_([ACGT])_b38$")[0]
    eq["alt"] = eq["variant_id"].str.extract(r"_([ACGT])_([ACGT])_b38$")[1]

    keep = set(eq["snp"])
    out = pd.read_csv(
        OUTCOME,
        sep="\t",
        usecols=["RSID", "EA", "NEA", "BETA", "SE"],
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
    print("共享 SNP 数：", len(merged))

    l1 = [logabf(z, se) for z, se in zip(merged["z"], merged["se"])]
    l2 = [logabf(b / s, s) for b, s in zip(merged["BETA"], merged["SE"])]

    lH0 = 0.0
    lH1 = logsumexp(l1) + math.log(P1)
    lH2 = logsumexp(l2) + math.log(P2)
    lH4 = logsumexp([a + b for a, b in zip(l1, l2)]) + math.log(P12)
    all_pairs = logsumexp(l1) + logsumexp(l2)
    same_sum = logsumexp([a + b for a, b in zip(l1, l2)])
    lH3 = math.log(max(1e-300, math.exp(all_pairs) - math.exp(same_sum))) + math.log(P1 * P2)

    logs = [lH0, lH1, lH2, lH3, lH4]
    m = max(logs)
    probs = [math.exp(v - m) for v in logs]
    total = sum(probs)
    pp = [p / total for p in probs]
    names = ["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"]
    for n, v in zip(names, pp):
        print(f"{n}: {v:.4f}")


if __name__ == "__main__":
    main()

