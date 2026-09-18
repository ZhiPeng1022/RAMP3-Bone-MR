# BMI -> eBMD 两样本 MR（本地 eBMD + OpenGWAS BMI 工具）

import math
import os

import numpy as np
import pandas as pd
from scipy import stats


DATA = r"data/derived"
BMI = pd.read_csv(os.path.join(DATA, "bmi_instruments.csv"))
EBMD = os.path.join(DATA, "ebmd_unpacked", "Biobank2-British-Bmd-As-C-Gwas-SumStats.txt.gz")


def ivw(beta, se):
    w = 1.0 / se ** 2
    b = np.sum(w * beta) / np.sum(w)
    se_b = math.sqrt(1.0 / np.sum(w))
    return b, se_b, 2 * stats.norm.sf(abs(b / se_b))


def main():
    keep = set(BMI["SNP"])
    out = pd.read_csv(EBMD, sep="\t", usecols=["RSID", "EA", "NEA", "BETA", "SE"], dtype={"RSID": str})
    out = out[out["RSID"].isin(keep)].copy()
    m = BMI.merge(out, left_on="SNP", right_on="RSID")

    same = m["effect_allele.exposure"].str.upper() == m["EA"].str.upper()
    flip = (m["effect_allele.exposure"].str.upper() == m["NEA"].str.upper()) & (
        m["other_allele.exposure"].str.upper() == m["EA"].str.upper()
    )
    m.loc[flip, "BETA"] = -m.loc[flip, "BETA"]
    m = m[same | flip].copy()
    print("匹配 SNP：", len(m))
    if len(m) == 0:
        return
    b, se, p = ivw(m["BETA"].astype(float).values, m["SE"].astype(float).values)
    print(f"BMI -> eBMD: IVW beta={b:.5f} SE={se:.5f} P={p:.4g}")


if __name__ == "__main__":
    main()

