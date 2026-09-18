# 从 eBMD 提取 chr7 区域（b37）全部 SNP，用于区域级共定位
import csv
import gzip
import os

DATA = r"data/derived"
EBMD = os.path.join(DATA, "ebmd_unpacked", "Biobank2-British-Bmd-As-C-Gwas-SumStats.txt.gz")
OUT = os.path.join(DATA, "ebmd_chr7_region.csv")
CHR = "7"
BP_LO = 44_000_000
BP_HI = 46_500_000


def main():
    rows = []
    total = 0
    with gzip.open(EBMD, "rt", encoding="utf-8", errors="replace", newline="") as src:
        reader = csv.reader(src, delimiter="\t")
        header = next(reader)
        col = {name: i for i, name in enumerate(header)}
        for line in reader:
            total += 1
            if line[col["CHR"]] != CHR:
                continue
            bp = int(line[col["BP"]])
            if BP_LO <= bp <= BP_HI:
                rows.append([
                    line[col["RSID"]], bp, line[col["EA"]], line[col["NEA"]],
                    line[col["BETA"]], line[col["SE"]], line[col["P"]], line[col["INFO"]]
                ])
    print("总行数：", total, "；chr7 区域行数：", len(rows))
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["RSID", "BP", "EA", "NEA", "BETA", "SE", "P", "INFO"])
        writer.writerows(rows)
    print("已保存：", OUT)


if __name__ == "__main__":
    main()
