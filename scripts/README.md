# Script notes

The scripts were copied from the final analysis workspace without changing
their statistical methods.

Most R scripts define a `DATA` variable near the top. Set this variable to:

```text
<repository-root>/data/derived
```

Some scripts also define an output directory or use interactive OpenGWAS
access. These paths and tokens must be configured locally.

The final figure scripts under `scripts/06_figures/` use repository-relative
paths and write their outputs to `figures/final/`.
