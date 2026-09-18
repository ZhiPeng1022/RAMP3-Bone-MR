# RAMP3-Bone-MR

Reproducible analysis code and derived data for the manuscript:

> RAMP3/amylin receptor signaling and cortical bone: site-specific genetic
> associations are largely explained by metabolic pleiotropy

## Repository contents

- `scripts/01_data_preparation/`: retrieval and extraction of public summary
  statistics and expression data.
- `scripts/02_primary_mr/`: primary Mendelian randomization, replication,
  reverse MR, bone-marker and tissue-specific analyses.
- `scripts/03_colocalization/`: colocalization and SuSiE analyses.
- `scripts/04_pleiotropy_and_mediation/`: MVMR, Steiger and mediation analyses.
- `scripts/05_single_cell/`: single-cell expression processing.
- `scripts/06_figures/`: final figure-generation scripts.
- `data/derived/`: small derived tables used by the analysis and figure code.
- `data/raw/`: no raw third-party data are redistributed.
- `results/`: final result tables and summary outputs.
- `figures/final/`: submitted PDF and TIFF figures.

## Data access

Raw genome-wide association, eQTL, pQTL and single-cell data are not stored in
this repository. Accession numbers, repository links and download notes are
provided in `data/README.md`.

## Software

Python dependencies are listed in `requirements.txt`.

R packages and installation instructions are provided in
`install_R_packages.R`.

The analysis scripts were developed with Python 3.12 and R 4.6.1. Some
third-party R packages should be installed from their current GitHub
repositories, as documented in `install_R_packages.R`.

## Running the code

Figure generation can be run from the repository root:

```bash
python scripts/06_figures/figure_1_forest.py
python scripts/06_figures/figure_2_singlecell.py
python scripts/06_figures/figure_3_site_bmd.py
```

The R scripts retain their original analysis structure. Before running them,
set the `DATA` variable at the top of each script to the local path containing
the files in `data/derived/`. Scripts that use OpenGWAS require a valid
OpenGWAS JWT token in the local R environment.

Some scripts use public APIs and will download data again when run. Output
paths may also need to be set to local writable directories.

## Derived data

The files in `data/derived/` are extracted subsets or intermediate inputs used
by the analysis scripts. They are provided to make the analyses easier to
inspect and rerun. Users should consult the original repositories for the
authoritative source data and their terms of use.

## Citation

Citation metadata are provided in `CITATION.cff`.

- Version DOI for `v1.0.0`:
  https://doi.org/10.5281/zenodo.22832847
- Concept DOI for all versions:
  https://doi.org/10.5281/zenodo.22832846

## License

This project is released under the MIT License. See `LICENSE` for details.
