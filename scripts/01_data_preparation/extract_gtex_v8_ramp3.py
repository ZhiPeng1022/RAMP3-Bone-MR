# 从 GTEx v8 independent eQTL tar 提取 Whole_Blood 的 RAMP3 slope/slope_se
# 用途：与 v10 API 近似值交叉复核（v8 版本差异会在论文中说明）
import csv
import gzip
import os
import sys
import tarfile

DATA = r"data/derived"
TAR = os.path.join(DATA, "GTEx_Analysis_v8_eQTL_independent.tar")
OUT = os.path.join(DATA, "gtex_v8_ramp3_independent.csv")
RAMP3 = "ENSG00000122679.8"


def main():
    with tarfile.open(TAR, "r") as tf:
        member = None
        for m in tf.getmembers():
            if m.isfile() and "Whole_Blood" in os.path.basename(m.name):
                member = m
                break
        if member is None:
            print("未找到 Whole_Blood 文件")
            sys.exit(1)
        print("文件：", member.name, member.size)
        with tf.extractfile(member) as src:
            with gzip.open(src, "rt", encoding="utf-8", errors="replace", newline="") as gz:
                reader = csv.DictReader(gz, delimiter="\t")
                col = {k.lower(): k for k in reader.fieldnames}
                print("表头：", reader.fieldnames)
                gene_col = col.get("gene_id") or col.get("gencode_id")
                if gene_col is None:
                    print("无法识别基因列")
                    sys.exit(1)
                rows = []
                for row in reader:
                    if row[gene_col] == RAMP3:
                        rec = {"snp": row[col.get("variant_id", "")] if col.get("variant_id") else "",
                               "gene_id": row[gene_col]}
                        for k in ("slope", "slope_se", "pvalue", "pval_nominal"):
                            if col.get(k):
                                rec[k] = row[col[k]]
                        rows.append(rec)
                print("RAMP3 行数：", len(rows))
                if rows:
                    keys = list(rows[0].keys())
                    with open(OUT, "w", newline="", encoding="utf-8") as fh:
                        writer = csv.DictWriter(fh, fieldnames=keys)
                        writer.writeheader()
                        writer.writerows(rows)
                    print("已保存：", OUT)
                    for r in rows:
                        print(r)


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
