# 从 GTEx v10 eQTL tar 中提取 Whole_Blood 的 RAMP3 slope/slope_se
# 用法：等 GTEx_Analysis_v10_eQTL.tar 下载完整后运行
import csv
import gzip
import io
import os
import sys
import tarfile

DATA = r"data/derived"
TAR = os.path.join(DATA, "GTEx_Analysis_v10_eQTL.tar")
EXPECTED = 2561047718
OUT = os.path.join(DATA, "gtex_ramp3_slope_recheck.csv")
RAMP3 = "ENSG00000122679.8"


def find_whole_blood(tf):
    cands = []
    for m in tf.getmembers():
        if not m.isfile():
            continue
        name = os.path.basename(m.name).lower()
        if "whole_blood" in name and ("signif" in name or "allpairs" in name or "egenes" in name):
            cands.append(m)
    cands.sort(key=lambda m: m.size)
    return cands


def read_pairs(path):
    with open(path, "rb") as fh:
        with gzip.open(fh, "rt", encoding="utf-8", errors="replace", newline="") as src:
            reader = csv.reader(src, delimiter="\t")
            header = next(reader)
            col = {name.lower(): i for i, name in enumerate(header)}
            gene_col = next((k for k in ("gencode_id", "gene_id", "gene") if k in col), None)
            slope_col = next((k for k in ("slope", "slope.se", "slope_se") if k in col), None)
            slope_se_col = next((k for k in ("slope_se", "slope.se", "se") if k in col), None)
            p_col = next((k for k in ("pvalue", "p.value", "pval") if k in col), None)
            if gene_col is None or slope_col is None:
                print("无法识别的列：", header)
                return []
            rows = []
            for row in reader:
                if row[gene_col] == RAMP3:
                    rows.append({
                        "variant_id": row[col["variant_id"]] if "variant_id" in col else "",
                        "snp": row[col["snp_id"]] if "snp_id" in col else "",
                        "slope": row[slope_col],
                        "slope_se": row[slope_se_col] if slope_se_col else "",
                        "p": row[p_col] if p_col else "",
                    })
            return rows


def main():
    if not os.path.exists(TAR) or os.path.getsize(TAR) < EXPECTED:
        print("tar 未完整：", os.path.getsize(TAR) if os.path.exists(TAR) else "不存在",
              "/", EXPECTED)
        sys.exit(1)
    with tarfile.open(TAR, "r") as tf:
        cands = find_whole_blood(tf)
        if not cands:
            print("tar 中未找到 Whole_Blood eQTL 文件")
            for m in tf.getmembers()[:20]:
                print(m.name, m.size)
            sys.exit(1)
        member = cands[0]
        print("使用文件：", member.name, member.size)
        with tf.extractfile(member) as src:
            with gzip.open(src, "rt", encoding="utf-8", errors="replace", newline="") as gz:
                reader = csv.reader(gz, delimiter="\t")
                header = next(reader)
                print("表头：", header)
                col = {name.lower(): i for i, name in enumerate(header)}
                gene_col = next((k for k in ("gencode_id", "gene_id", "gene") if k in col), None)
                slope_col = next((k for k in ("slope", "slope.se", "slope_se") if k in col), None)
                slope_se_col = next((k for k in ("slope_se", "slope.se", "se") if k in col), None)
                p_col = next((k for k in ("pvalue", "p.value", "pval") if k in col), None)
                rows = []
                for row in reader:
                    if gene_col and row[gene_col] == RAMP3:
                        rows.append({
                            "snp": row[col["snp_id"]] if "snp_id" in col else "",
                            "variant_id": row[col["variant_id"]] if "variant_id" in col else "",
                            "slope": row[slope_col] if slope_col else "",
                            "slope_se": row[slope_se_col] if slope_se_col else "",
                            "p": row[p_col] if p_col else "",
                        })
                print("RAMP3 行数：", len(rows))
                if rows:
                    with open(OUT, "w", newline="", encoding="utf-8") as fh:
                        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
                        writer.writeheader()
                        writer.writerows(rows)
                    print("已保存：", OUT)
                    for r in rows[:10]:
                        print(r)


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
