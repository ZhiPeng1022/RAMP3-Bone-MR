# 从 GTEx API 拉取的 RAMP3 eQTL 表中生成 MR 暴露文件（每组织 top SNP）

import csv
import os
import re
import argparse

import pandas as pd
from scipy import stats


SRC = os.path.join(r"data/derived", "gtex_ramp3_eqtl.csv")
OUT = os.path.join(r"data/derived", "ramp3_top_eqtl_exposure.txt")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tissue", default=None, help="指定组织；缺省为每组织 top1")
    parser.add_argument("--n", type=int, default=1, help="每个组织取 top n 个 SNP")
    parser.add_argument("--out", default=OUT, help="输出文件")
    args = parser.parse_args()

    df = pd.read_csv(SRC)
    df = df.dropna(subset=["p", "nes"])
    df = df[df["approx_F"] >= 10]

    rows = []
    if args.tissue:
        groups = [(args.tissue, df[df["tissue"] == args.tissue].sort_values("p").head(args.n))]
    else:
        groups = [(t, g.sort_values("p").head(args.n)) for t, g in df.groupby("tissue")]

    for tissue, group in groups:
        if len(group) == 0:
            continue
        for _, top in group.iterrows():
            p = float(top["p"])
            nes = float(top["nes"])
            z = abs(stats.norm.ppf(p / 2.0))
            se = abs(nes / z) if z else float("nan")
            variant_id = str(top["variant_id"])
            m = re.search(r"_([ACGT])_([ACGT])_b38$", variant_id)
            ref, alt = (m.group(1), m.group(2)) if m else ("", "")
            rows.append({
                "rs_id": top["snp"],
                "beta": round(nes, 6),
                "se": round(se, 6),
                "alt": alt,
                "ref": ref,
                "tissue": tissue,
            })

    out_df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    out_df.to_csv(args.out, sep="\t", index=False)
    print("saved", len(out_df), "SNPs ->", args.out)
    print(out_df[["rs_id", "tissue", "beta", "se"]].to_string(index=False))


if __name__ == "__main__":
    main()
