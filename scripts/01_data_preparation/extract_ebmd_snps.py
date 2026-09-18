# 从 eBMD 大文件中按 rsID 流式提取子集，供 R 正式分析使用
# 避免把整个 4.7 亿行文件读入内存
import csv
import gzip
import os
import re
import sys

DATA = r"data/derived"
EBMD = os.path.join(DATA, "ebmd_unpacked", "Biobank2-British-Bmd-As-C-Gwas-SumStats.txt.gz")
OUT = os.path.join(DATA, "ebmd_extract_mvmr.csv")


def collect_snps():
    snps = set()

    # 1) RAMP3 全血所有 eQTL（MVMR 敏感性分析需要）
    eq_file = os.path.join(DATA, "gtex_ramp3_eqtl.csv")
    if os.path.exists(eq_file):
        with open(eq_file, encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                if row["tissue"] == "Whole_Blood":
                    snps.add(row["snp"].strip())

    # 2) RAMP3 全血 top eQTL
    top_file = os.path.join(DATA, "ramp3_whole_blood_exposure.txt")
    if os.path.exists(top_file):
        with open(top_file, encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                snps.add(row["rs_id"].strip())

    # 3) BMI 工具变量
    bmi_file = os.path.join(DATA, "bmi_instruments.csv")
    if os.path.exists(bmi_file):
        with open(bmi_file, encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                snps.add(row["SNP"].strip())

    # 4) T2D 工具变量
    t2d_file = os.path.join(DATA, "t2d_instruments.csv")
    if os.path.exists(t2d_file):
        with open(t2d_file, encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                snps.add(row["SNP"].strip())

    # 5) eQTLGen 工具变量（RAMP3 复核 + GLP1R/SLC5A2 对照）
    eqgen_file = os.path.join(DATA, "eqtlgen_instruments.csv")
    if os.path.exists(eqgen_file):
        with open(eqgen_file, encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                snps.add(row["SNP"].strip())

    return snps


def main():
    snps = collect_snps()
    print("需要提取的 SNP 总数：", len(snps))
    if not snps:
        print("没有 SNP，请先确认 bmi_instruments.csv / t2d_instruments.csv 存在")
        sys.exit(1)

    pattern = re.compile(r"\t(" + "|".join(re.escape(s) for s in sorted(snps, key=len, reverse=True)) + r")\t")
    found = 0
    with gzip.open(EBMD, "rt", encoding="utf-8", errors="replace", newline="") as src, \
         open(OUT, "w", newline="", encoding="utf-8") as dst:
        writer = csv.writer(dst)
        writer.writerow(["RSID", "EA", "NEA", "BETA", "SE", "P", "N"])
        for line in src:
            m = pattern.search(line)
            if not m:
                continue
            fields = line.rstrip("\r\n").split("\t")
            rs = fields[1]
            snps.discard(rs)
            writer.writerow([
                rs,
                fields[4],
                fields[5],
                fields[8],
                fields[9],
                fields[10],
                fields[13],
            ])
            found += 1
            if not snps:
                break

    print("提取完成：", found, "个 SNP ->", OUT)
    miss = sorted(snps)
    print("未匹配：", len(miss))
    print(miss[:50])


if __name__ == "__main__":
    main()
