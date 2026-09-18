# 从 GTEx v10 按位置接口获取 RAMP3 cis 区域全血 eQTL（含非显著）
import csv
import json
import os
import time
import urllib.request

DATA = r"data/derived"
OUT = os.path.join(DATA, "gtex_ramp3_region_eqtl.csv")
RAMP3 = "ENSG00000122679.8"
BASE = "https://gtexportal.org/api/v2/association/singleTissueEqtlByLocation"
START = 44_200_000
END = 46_100_000


def fetch(page, size):
    url = (f"{BASE}?tissueSiteDetailId=Whole_Blood&start={START}&end={END}"
           f"&chromosome=chr7&datasetId=gtex_v10&itemsPerPage={size}&page={page}")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read().decode("utf-8")).get("singleTissueEqtl", [])


def main():
    rows = []
    size = 1000
    items = fetch(0, size)
    print("items", len(items))
    rows = [x for x in items if x.get("gencodeId") == RAMP3]

    # 去重
    seen = set()
    uniq = []
    for r in rows:
        key = r["variantId"]
        if key not in seen:
            seen.add(key)
            uniq.append(r)
    print("RAMP3 区域 eQTL 总数：", len(uniq))
    if uniq:
        with open(OUT, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=["snpId", "variantId", "pos", "nes", "pValue"])
            writer.writeheader()
            for r in uniq:
                writer.writerow({
                    "snpId": r["snpId"], "variantId": r["variantId"], "pos": r["pos"],
                    "nes": r["nes"], "pValue": r["pValue"],
                })
        print("已保存：", OUT)
        sig = [r for r in uniq if r["pValue"] < 1e-4]
        print("P<1e-4 显著行：", len(sig))


if __name__ == "__main__":
    main()
