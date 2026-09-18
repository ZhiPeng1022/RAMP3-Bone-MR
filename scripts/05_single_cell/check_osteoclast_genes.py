# 检查破骨细胞数据集 features 中是否包含 RAMP3/CALCR

import gzip
import os


DATA = r"data/derived"
GENES = ["RAMP3", "CALCR", "RAMP1", "RAMP2", "IAPP", "CALCB", "TNFRSF11A", "CTSK"]


def main():
    feature_files = []
    for root, _, files in os.walk(DATA):
        for f in files:
            if f.endswith("features.tsv.gz"):
                feature_files.append(os.path.join(root, f))

    print("features 文件数：", len(feature_files))
    for path in sorted(feature_files):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            lines = fh.readlines()
        found = {}
        for line in lines:
            parts = line.rstrip("\n").split("\t")
            name = parts[1] if len(parts) > 1 else parts[0]
            upper = name.upper()
            for g in GENES:
                if upper == g:
                    found[g] = line.strip()
        print(os.path.basename(os.path.dirname(path)), "| 基因数:", len(lines), "| 命中:", found)


if __name__ == "__main__":
    main()

