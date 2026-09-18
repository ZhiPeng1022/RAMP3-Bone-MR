# 拉取 GLP1R / SLC5A2 的 GTEx eQTL（对照工具变量）

import csv
import json
import os
import urllib.request

from scipy import stats


OUT = os.path.join(r"data/derived", "control_receptors_eqtl.csv")
GENES = {"GLP1R": None, "SLC5A2": None}
TISSUES = [
    "Liver", "Whole_Blood", "Adipose_Subcutaneous", "Pancreas", "Kidney_Cortex",
    "Muscle_Skeletal", "Thyroid", "Small_Intestine_Terminal_Ileum", "Brain_Cortex",
    "Colon_Sigmoid", "Stomach", "Artery_Aorta", "Heart_Left_Ventricle", "Spleen", "Lung",
]


def get_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    rows = []
    for gene in GENES:
        gene_data = get_json("https://gtexportal.org/api/v2/reference/gene?geneId=" + gene)
        items = gene_data.get("data", [])
        if not items:
            print(gene, "NOT FOUND")
            continue
        gid = items[0]["gencodeId"]
        GENES[gene] = gid
        for tissue in TISSUES:
            url = (
                "https://gtexportal.org/api/v2/association/singleTissueEqtl"
                f"?gencodeId={gid}&tissueSiteDetailId={tissue}&datasetId=gtex_v10&itemsPerPage=5"
            )
            try:
                data = get_json(url)
                top = data.get("data", [])[:1]
                if not top:
                    print(gene, tissue, "NO_EQTL")
                    continue
                item = top[0]
                p = item.get("pValue")
                z = abs(stats.norm.ppf(p / 2.0)) if p and p > 0 else None
                rows.append({
                    "gene": gene,
                    "tissue": tissue,
                    "snp": item.get("snpId"),
                    "variant_id": item.get("variantId"),
                    "p": p,
                    "nes": item.get("nes"),
                    "approx_F": round(z ** 2, 2) if z else None,
                })
                print(gene, tissue, item.get("snpId"), p, "F~", round(z ** 2, 2) if z else None)
            except Exception as exc:
                print(gene, tissue, "ERROR", type(exc).__name__)

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["gene", "tissue", "snp", "variant_id", "p", "nes", "approx_F"])
        writer.writeheader()
        writer.writerows(rows)
    print("saved ->", OUT)


if __name__ == "__main__":
    main()
