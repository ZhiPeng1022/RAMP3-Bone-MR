# 解析 HPA RAMP3 单细胞表达页面的内嵌 JSON，保存 CSV

import csv
import os
import re


OUT = os.path.join(r"data/derived", "hpa_ramp3_singlecell.csv")
HTML = os.path.join(os.environ["TEMP"], "hpa_ramp3_sc.html")


def main():
    with open(HTML, encoding="utf-8", errors="replace") as f:
        text = f.read()

    pattern = re.compile(r'"label":"([^"]+)","legend":"([^"]+)","value":"([^"]+)"')
    seen = set()
    rows = []
    for m in pattern.finditer(text):
        label, legend, value = m.group(1), m.group(2), m.group(3)
        key = (label, legend)
        if key in seen:
            continue
        seen.add(key)
        rows.append({"cell_type": label, "cluster": legend, "nTPM": value})

    rows.sort(key=lambda r: -float(r["nTPM"]))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["cell_type", "cluster", "nTPM"])
        writer.writeheader()
        writer.writerows(rows)
    print("saved", len(rows), "rows ->", OUT)
    print("--- top 25 ---")
    for r in rows[:25]:
        print(f"{r['cell_type']:28s} {r['cluster']:20s} {r['nTPM']}")


if __name__ == "__main__":
    main()

