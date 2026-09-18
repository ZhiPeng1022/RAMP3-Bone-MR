# Data accessions and download notes

Raw data are not redistributed in this repository. The main public sources are:

| Resource | Identifier | Link |
|---|---|---|
| GTEx v10 | GTEx Portal | https://gtexportal.org |
| eQTLGen Phase I | eQTLGen | https://www.eqtlgen.org |
| Morris 2019 eBMD | GEFOS | http://www.gefos.org |
| Serum calcium | GCST90025990 | https://www.ebi.ac.uk/gwas/studies/GCST90025990 |
| Serum phosphate | GCST90025948 | https://www.ebi.ac.uk/gwas/studies/GCST90025948 |
| BMI | ieu-b-40 | https://gwas.mrcieu.ac.uk/datasets/ieu-b-40/ |
| T2D | ebi-a-GCST006867 | https://www.ebi.ac.uk/gwas/studies/GCST006867 |
| eBMD replication | ebi-a-GCST90014022 | https://www.ebi.ac.uk/gwas/studies/GCST90014022 |
| Heel BMD | ebi-a-GCST90029004 | https://www.ebi.ac.uk/gwas/studies/GCST90029004 |
| Heel BMD | ukb-b-8875 | https://gwas.mrcieu.ac.uk/datasets/ukb-b-8875/ |
| Forearm BMD | ieu-a-977 | https://gwas.mrcieu.ac.uk/datasets/ieu-a-977/ |
| Femoral neck BMD | ieu-a-980 | https://gwas.mrcieu.ac.uk/datasets/ieu-a-980/ |
| Lumbar spine BMD | ieu-a-982 | https://gwas.mrcieu.ac.uk/datasets/ieu-a-982/ |
| FinnGen | FinnGen public release | https://www.finngen.fi |
| Single-cell data | GSE303003 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303003 |
| Mouse osteoclast data | GSE310113 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310113 |
| Human Protein Atlas | HPA | https://www.proteinatlas.org |

## Files in `data/derived/`

- `ramp3_whole_blood_exposure.txt`: primary whole-blood RAMP3 instruments.
- `ramp3_top_eqtl_exposure.txt`: top tissue-specific instruments.
- `gtex_ramp3_eqtl.csv`: GTEx v10 RAMP3 cis-eQTL extract.
- `eqtlgen_instruments.csv`: eQTLGen replication instruments.
- `bmi_instruments.csv` and `t2d_instruments.csv`: MVMR instruments.
- `ebmd_extract_mvmr.csv`: eBMD effects used in MVMR.
- `opengwas_mediators.csv`: extracted mediator effects.
- `eqtlgen_core_region_outcomes.csv`, `forearm_region_outcomes.csv` and
  `ebmd_chr7_region.csv`: regional colocalization inputs.
- `hpa_ramp3_singlecell.csv`: HPA RAMP3 cell-type expression extract.
- Additional files are intermediate extracts used by specific sensitivity
  analyses.

The derived files remain subject to the terms of the original repositories.
