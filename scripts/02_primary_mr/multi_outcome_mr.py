# RAMP3 全血 eQTL 对血钙/血磷的多 SNP IVW

import math
import os

import numpy as np
import pandas as pd
from scipy import stats


DATA = r"data/derived"
EQ = pd.read_csv(os.path.join(DATA, "gtex_ramp3_eqtl.csv"))
OUT = pd.read_csv(os.path.join(DATA, "opengwas_outcomes_full.csv"))


def ivw(beta, se):
    w = 1.0 / se ** 2
    b = np.sum(w * beta) / np.sum(w)
    se_b = math.sqrt(1.0 / np.sum(w))
    return b, se_b, 2 * stats.norm.sf(abs(b / se_b))


def main():
    eq = EQ[EQ["tissue"] == "Whole_Blood"].dropna(subset=["p", "nes"]).copy()
    eq["se"] = abs(eq["nes"]) / abs(stats.norm.ppf(eq["p"] / 2.0))
    eq["ref"] = eq["variant_id"].str.extract(r"_([ACGT])_([ACGT])_b38$")[0]
    eq["alt"] = eq["variant_id"].str.extract(r"_([ACGT])_([ACGT])_b38$")[1]

    for outcome in ["serum_calcium", "serum_phosphate"]:
        o = OUT[OUT["id.outcome"] == outcome].copy()
        m = eq.merge(o, left_on="snp", right_on="SNP")
        if len(m) == 0:
            print(outcome, "无匹配 SNP")
            continue
        same = m["alt"].str.upper() == m["effect_allele.outcome"].str.upper()
        flip = (m["alt"].str.upper() == m["other_allele.outcome"].str.upper()) & (
            m["ref"].str.upper() == m["effect_allele.outcome"].str.upper()
        )
        m.loc[flip, "beta.outcome"] = -m.loc[flip, "beta.outcome"]
        m = m[same | flip].copy()
        if len(m) == 0:
            print(outcome, "对齐后无 SNP")
            continue
        b, se, p = ivw(m["beta.outcome"].astype(float).values, m["se.outcome"].astype(float).values)
        print(f"{outcome}: n_snp={len(m)} IVW beta={b:.5f} SE={se:.5f} P={p:.4g}")


if __name__ == "__main__":
    main()

