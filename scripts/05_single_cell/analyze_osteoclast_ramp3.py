# 统计破骨细胞数据集中 RAMP3/CALCR 等基因的阳性细胞数与表达

import argparse
import gzip
import os


DATA = r"data/derived"
GENES = ["RAMP3", "CALCR", "RAMP1", "RAMP2", "IAPP", "CALCB", "CTSK", "TNFRSF11A"]


def load_target_rows(features_path):
    target = {g.upper(): None for g in GENES}
    with gzip.open(features_path, "rt", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh, start=1):
            parts = line.rstrip("\n").split("\t")
            name = parts[1] if len(parts) > 1 else parts[0]
            if name.upper() in target:
                target[name.upper()] = i
    return target


def analyze_sample(features_path, matrix_path):
    target = load_target_rows(features_path)
    counts = {g: 0 for g in GENES}
    totals = {g: 0 for g in GENES}
    n_cells = 0
    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace") as fh:
        header_done = False
        for line in fh:
            line = line.strip()
            if not line or line.startswith("%"):
                continue
            parts = line.split()
            if not header_done:
                n_cells = int(parts[1])
                header_done = True
                continue
            row = int(parts[0])
            col = int(parts[1])
            val = float(parts[2])
            for g, r in target.items():
                if r == row:
                    counts[g] += 1
                    totals[g] += val
    return n_cells, counts, totals


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=None, help="只处理前 N 个样本")
    args = parser.parse_args()

    samples = []
    for root, _, files in os.walk(DATA):
        files = set(files)
        for f in sorted(files):
            if f.endswith(".matrix.mtx.gz"):
                prefix = f.replace(".matrix.mtx.gz", "")
                feats = os.path.join(root, prefix + ".features.tsv.gz")
                if not os.path.exists(feats):
                    feats = os.path.join(root, prefix + "_features.tsv.gz")
                if os.path.exists(feats):
                    samples.append((prefix, feats, os.path.join(root, f)))

    if args.limit:
        samples = samples[:args.limit]
    print(f"处理 {len(samples)} 个样本")
    for prefix, feats, matrix in samples:
        n_cells, counts, totals = analyze_sample(feats, matrix)
        print("=====", os.path.basename(prefix), "| cells:", n_cells)
        for g in GENES:
            if counts[g]:
                print(f"  {g}: 阳性细胞 {counts[g]} ({counts[g]/max(n_cells,1)*100:.2f}%) 总表达 {totals[g]:.1f}")


if __name__ == "__main__":
    main()
