# 从 GTEx API 拉取 RAMP3 eQTL 并保存（工具变量表）

import csv
import json
import os
import urllib.request

from scipy import stats


OUT = os.path.join(r"data/derived", "gtex_ramp3_eqtl.csv")
RAMP3 = "ENSG00000122679.8"
TISSUES = [
    "Adipose_Subcutaneous",
    "Liver",
    "Whole_Blood",
    "Muscle_Skeletal",
    "Thyroid",
    "Pancreas",
    "Kidney_Cortex",
    "Pituitary",
]


def fetch(tissue):
    url = (
        "https://gtexportal.org/api/v2/association/singleTissueEqtl"
        f"?gencodeId={RAMP3}&tissueSiteDetailId={tissue}"
        "&datasetId=gtex_v10&itemsPerPage=250"
    )
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode("utf-8")).get("data", [])


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    rows = []
    for tissue in TISSUES:
        try:
            items = fetch(tissue)
        except Exception as exc:
            print(tissue, "ERROR", type(exc).__name__, exc)
            continue
        for item in items:
            p = item.get("pValue")
            z = abs(stats.norm.ppf(p / 2.0)) if p and p > 0 else None
            rows.append({
                "tissue": tissue,
                "snp": item.get("snpId"),
                "variant_id": item.get("variantId"),
                "p": p,
                "nes": item.get("nes"),
                "approx_F": round(z ** 2, 2) if z else None,
            })

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["tissue", "snp", "variant_id", "p", "nes", "approx_F"])
        writer.writeheader()
        writer.writerows(rows)
    print("saved", len(rows), "rows ->", OUT)

    by_tissue = {}
    for row in rows:
        by_tissue.setdefault(row["tissue"], []).append(row)
    for tissue, items in by_tissue.items():
        items.sort(key=lambda x: x["p"] or 1)
        top = items[0]
        print(tissue, "|", top["snp"], "| p =", top["p"], "| F ~", top["approx_F"], "| n_snps =", len(items))


if __name__ == "__main__":
    main()

