# 简化 MVMR：RAMP3 + BMI + T2D 对 eBMD（加权多变量 IVW 近似）

import math
import os

import numpy as np
import pandas as pd
from scipy import stats


DATA = r"data/derived"
EQ = pd.read_csv(os.path.join(DATA, "gtex_ramp3_eqtl.csv"))
MED = pd.read_csv(os.path.join(DATA, "opengwas_mediators.csv"))
EBMD = os.path.join(DATA, "ebmd_unpacked", "Biobank2-British-Bmd-As-C-Gwas-SumStats.txt.gz")


def main():
    eq = EQ[EQ["tissue"] == "Whole_Blood"].dropna(subset=["p", "nes"]).copy()
    eq["se"] = abs(eq["nes"]) / abs(stats.norm.ppf(eq["p"] / 2.0))
    eq["ref"] = eq["variant_id"].str.extract(r"_([ACGT])_([ACGT])_b38$")[0]
    eq["alt"] = eq["variant_id"].str.extract(r"_([ACGT])_([ACGT])_b38$")[1]
    eq = eq.rename(columns={"snp": "SNP", "nes": "beta_ramp3"})

    bmi = MED[MED["id.outcome"] == "BMI"][["SNP", "effect_allele.outcome", "other_allele.outcome", "beta.outcome"]].rename(
        columns={"effect_allele.outcome": "EA_BMI", "other_allele.outcome": "OA_BMI", "beta.outcome": "beta_BMI"})
    t2d = MED[MED["id.outcome"] == "T2D"][["SNP", "effect_allele.outcome", "other_allele.outcome", "beta.outcome"]].rename(
        columns={"effect_allele.outcome": "EA_T2D", "other_allele.outcome": "OA_T2D", "beta.outcome": "beta_T2D"})

    keep = set(eq["SNP"]) | set(bmi["SNP"]) | set(t2d["SNP"])
    out = pd.read_csv(EBMD, sep="\t", usecols=["RSID", "EA", "NEA", "BETA", "SE"], dtype={"RSID": str})
    out = out[out["RSID"].isin(keep)].rename(columns={"RSID": "SNP", "EA": "EA_EBMD", "NEA": "OA_EBMD", "BETA": "beta_EBMD", "SE": "se_EBMD"})

    m = eq[["SNP", "beta_ramp3", "se", "ref", "alt"]].merge(bmi, on="SNP", how="left").merge(t2d, on="SNP", how="left").merge(out, on="SNP", how="left")
    m = m.dropna(subset=["beta_EBMD", "se_EBMD"])

    def align(beta_col, ea_col, oa_col, alt_col, ref_col):
        same = m[ea_col].str.upper() == m[alt_col].str.upper()
        flip = (m[ea_col].str.upper() == m[ref_col].str.upper()) & (m[oa_col].str.upper() == m[alt_col].str.upper())
        m.loc[flip, beta_col] = -m.loc[flip, beta_col]
        return same | flip

    ok_bmi = align("beta_BMI", "EA_BMI", "OA_BMI", "alt", "ref")
    ok_t2d = align("beta_T2D", "EA_T2D", "OA_T2D", "alt", "ref")
    ok_ebmd = align("beta_EBMD", "EA_EBMD", "OA_EBMD", "alt", "ref")
    m = m[ok_bmi & ok_t2d & ok_ebmd].copy()

    m = m.fillna(0.0)
    X = m[["beta_ramp3", "beta_BMI", "beta_T2D"]].astype(float).values
    y = m["beta_EBMD"].astype(float).values
    w = 1.0 / m["se_EBMD"].astype(float).values ** 2

    W = np.diag(w)
    XtWX = X.T @ W @ X
    XtWy = X.T @ W @ y
    try:
        coef = np.linalg.solve(XtWX, XtWy)
    except np.linalg.LinAlgError:
        print("矩阵奇异，改用最小范数解")
        coef = np.linalg.lstsq(np.sqrt(W)[:, None] * X, np.sqrt(w) * y, rcond=None)[0]
    resid = y - X @ coef
    sigma2 = np.sum(w * resid ** 2) / (len(y) - X.shape[1])
    cov = sigma2 * np.linalg.inv(XtWX)
    se = np.sqrt(np.diag(cov))
    z = coef / se
    p = 2 * stats.norm.sf(np.abs(z))
    for name, c, s, pp in zip(["RAMP3", "BMI", "T2D"], coef, se, p):
        print(f"{name}: beta={c:.6f} SE={s:.6f} P={pp:.4g}")
    print("SNP 数：", len(m))


if __name__ == "__main__":
    main()

