# 从 GCST90618485（RAMP3 血清 pQTL，European n=2,338）提取 RAMP3 区域 SNP
import csv
import gzip
import math
import os

DATA = r"data/derived"
SRC = os.path.join(DATA, "GCST90618485.tsv.gz")
OUT = os.path.join(DATA, "ramp3_pqtl_region.csv")
CHR = "7"
POS_LO = 44_000_000
POS_HI = 46_500_000


def main():
    rows = []
    total = 0
    with gzip.open(SRC, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            total += 1
            if row["chromosome"] != CHR:
                continue
            pos = int(row["base_pair_location"])
            if POS_LO <= pos <= POS_HI:
                neg = float(row["neg_log_10_p_value"])
                rows.append({
                    "chrom": row["chromosome"],
                    "pos": pos,
                    "effect_allele": row["effect_allele"],
                    "other_allele": row["other_allele"],
                    "beta": float(row["beta"]),
                    "se": float(row["standard_error"]),
                    "eaf": float(row["effect_allele_frequency"]) if row["effect_allele_frequency"] else "",
                    "neg_log10_p": neg,
                    "p": 10 ** (-neg),
                    "n": row["n"],
                })
    print("总行数：", total, "；区域行数：", len(rows))
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    sig = [r for r in rows if r["p"] < 1e-5]
    sig.sort(key=lambda r: r["neg_log10_p"], reverse=True)
    print("P<1e-5 显著 SNP：", len(sig))
    for r in sig[:15]:
        print(r["chrom"], r["pos"], r["effect_allele"], r["other_allele"],
              "beta=", round(r["beta"], 5), "se=", round(r["se"], 5),
              "neg_log10_p=", round(r["neg_log10_p"], 2))


if __name__ == "__main__":
    main()
