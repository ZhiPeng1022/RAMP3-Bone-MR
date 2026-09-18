# 从 GWAS Catalog BMD summary stats 提取 chr7 核心区域（b38 45.0-45.4 Mb）
import csv
import gzip
import os
import sys

DATA = r"data/derived"
BP_LO = 45_000_000
BP_HI = 45_400_000


def main():
    acc = sys.argv[1]
    src = os.path.join(DATA, f"{acc}.tsv.gz")
    out = os.path.join(DATA, f"{acc}_chr7_region.csv")
    if not os.path.exists(src):
        print("文件不存在：", src)
        return
    rows = []
    with gzip.open(src, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            chrom = row.get("chromosome", "")
            if chrom not in ("7", "chr7"):
                continue
            bp = int(row.get("base_pair_location", 0) or 0)
            if BP_LO <= bp <= BP_HI:
                rows.append({
                    "variant_id": row.get("variant_id", ""),
                    "rsid": row.get("rsid", "") or row.get("RSID", ""),
                    "bp": bp,
                    "ea": row.get("effect_allele", ""),
                    "oa": row.get("other_allele", ""),
                    "beta": row.get("beta", ""),
                    "se": row.get("standard_error", ""),
                    "p": row.get("p_value", ""),
                    "eaf": row.get("effect_allele_frequency", ""),
                    "n": row.get("n", ""),
                })
    print(acc, "区域行数：", len(rows))
    if rows:
        with open(out, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print("已保存：", out)


if __name__ == "__main__":
    main()
