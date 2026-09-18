from pathlib import Path
from statistics import mean

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
from PIL import Image
import pandas as pd


matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

BASE = Path(__file__).resolve().parents[2]
OUTPUT_DIR = BASE / "figures" / "final"

HPA_DATA = BASE / "data" / "derived" / "hpa_ramp3_singlecell.csv"

PREVIEW_PATH = OUTPUT_DIR / "figure_2_preview.png"

HPA_LABELS = {
    "Endothelial cells": "Endothelial cells\n(HPA label 1)",
    "Endothelial cell": "Endothelial cells\n(HPA label 2)",
    "Renal nephron cells": "Renal nephron cells",
    "Stellate cells": "Stellate cells",
    "Stromal cells": "Stromal cells",
    "Spermatogenic cell types": "Spermatogenic cells",
    "Oligodendrocyte": "Oligodendrocytes",
    "Deep-layer corticothalamic and 6b": "Deep-layer\ncorticothalamic and 6b",
    "Barrier epithelial cells": "Barrier epithelial\ncells",
    "Vascular associated smooth muscle cell": "Vascular associated\nsmooth muscle cells",
    "Mononuclear phagocytes": "Mononuclear\nphagocytes",
    "Neuronal cells": "Neuronal cells",
}

SAMPLES = [
    "GSM9115504_Cyst_1",
    "GSM9115504_Cyst_2",
    "GSM9115504_MRONJ_1",
    "GSM9115504_MRONJ_3",
    "GSM9115505_Cyst_3",
    "GSM9115505_Cyst_4",
    "GSM9115505_MRONJ_4",
    "GSM9115505_MRONJ_5",
]

SAMPLE_TYPES = [
    "Cyst",
    "Cyst",
    "MRONJ",
    "MRONJ",
    "Cyst",
    "Cyst",
    "MRONJ",
    "MRONJ",
]

# Per-sample positive-cell percentages recalculated from the eight GSE303003
# count matrices. A value of 0.00 means that no positive cell was detected.
SAMPLE_PERCENTAGES = {
    "CTSK": [43.48, 21.37, 11.68, 13.98, 10.40, 13.18, 7.53, 13.26],
    "RAMP3": [5.69, 23.14, 8.51, 10.69, 4.38, 2.70, 14.46, 9.16],
    "RAMP2": [6.32, 20.37, 5.09, 3.95, 3.03, 1.84, 7.75, 3.11],
    "CALCB": [1.58, 0.29, 0.17, 0.27, 2.89, 0.75, 0.60, 0.39],
    "RAMP1": [2.13, 3.37, 2.96, 0.14, 0.41, 0.71, 0.36, 0.27],
    "IAPP": [0.03, 0.14, 0.01, 0.02, 0.05, 0.00, 0.11, 0.03],
    "CALCR": [0.17, 0.40, 0.03, 0.05, 0.08, 0.06, 0.04, 0.00],
}

GENE_ORDER = ["CTSK", "RAMP3", "RAMP2", "CALCB", "RAMP1", "IAPP", "CALCR"]


def build_figure(output_path):
    hpa = pd.read_csv(HPA_DATA)
    hpa = hpa.sort_values("nTPM", ascending=False).head(12).copy()
    hpa["plot_label"] = hpa["cell_type"].map(HPA_LABELS)

    figure, (left_axis, right_axis) = plt.subplots(
        1,
        2,
        figsize=(10.2, 5.6),
        dpi=300,
        gridspec_kw={"width_ratios": [1.18, 1.0]},
    )
    figure.patch.set_facecolor("white")

    # Left panel: HPA consensus expression by cell type.
    hpa_plot = hpa.sort_values("nTPM", ascending=True)
    left_axis.barh(
        range(len(hpa_plot)),
        hpa_plot["nTPM"],
        height=0.66,
        color="#718EA4",
        edgecolor="#4E6578",
        linewidth=0.6,
    )
    left_axis.set_yticks(
        range(len(hpa_plot)),
        hpa_plot["plot_label"],
    )
    left_axis.set_xlim(0, 455)
    left_axis.set_xlabel("RAMP3 nTPM", fontsize=10)
    left_axis.set_title("RAMP3 single-cell expression", fontsize=12, pad=10)
    left_axis.grid(axis="x", color="#E6E6E6", linewidth=0.8)
    left_axis.set_axisbelow(True)
    left_axis.tick_params(axis="y", labelsize=9.2, length=0, pad=7)
    left_axis.tick_params(axis="x", labelsize=8.8)
    left_axis.spines["top"].set_visible(False)
    left_axis.spines["right"].set_visible(False)

    for y_position, value in enumerate(hpa_plot["nTPM"]):
        left_axis.text(
            value + 5,
            y_position,
            f"{value:.1f}",
            va="center",
            ha="left",
            fontsize=7.8,
            color="#333333",
        )

    # Right panel: sample-level positive-cell percentages in GSE303003.
    y_positions = list(range(len(GENE_ORDER)))[::-1]
    offsets = [-0.21, -0.15, -0.09, -0.03, 0.03, 0.09, 0.15, 0.21]

    for y_position, gene in zip(y_positions, GENE_ORDER):
        values = SAMPLE_PERCENTAGES[gene]
        minimum = min(values)
        maximum = max(values)
        average = mean(values)

        right_axis.plot(
            [minimum, maximum],
            [y_position, y_position],
            color="#9C9C9C",
            linewidth=1.5,
            zorder=1,
        )

        for offset, value, sample_type in zip(offsets, values, SAMPLE_TYPES):
            if sample_type == "Cyst":
                right_axis.scatter(
                    [value],
                    [y_position + offset],
                    s=20,
                    marker="o",
                    facecolors="none",
                    edgecolors="#5D7F94",
                    linewidths=0.9,
                    zorder=2,
                )
            else:
                right_axis.scatter(
                    [value],
                    [y_position + offset],
                    s=21,
                    marker="s",
                    facecolors="#A58B7E",
                    edgecolors="#7E665C",
                    linewidths=0.6,
                    zorder=2,
                )

        right_axis.scatter(
            [average],
            [y_position],
            s=54,
            marker="D",
            facecolors="#3F3F3F",
            edgecolors="white",
            linewidths=0.7,
            zorder=3,
        )

    right_axis.set_yticks(y_positions, GENE_ORDER)
    right_axis.set_ylim(-0.65, len(GENE_ORDER) - 0.25)
    right_axis.set_xlim(0, 46)
    right_axis.set_xticks(range(0, 41, 10))
    right_axis.set_xlabel("Positive cells (%)", fontsize=10)
    right_axis.set_title(
        "Osteoclast-rich tissue (n=8)",
        fontsize=12,
        pad=10,
    )
    right_axis.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    right_axis.grid(axis="x", color="#E6E6E6", linewidth=0.8)
    right_axis.set_axisbelow(True)
    right_axis.tick_params(axis="y", labelsize=9.5, length=0, pad=7)
    right_axis.tick_params(axis="x", labelsize=8.8)
    right_axis.spines["top"].set_visible(False)
    right_axis.spines["right"].set_visible(False)

    legend_handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="None",
            markerfacecolor="none",
            markeredgecolor="#5D7F94",
            markeredgewidth=0.9,
            markersize=5,
            label="Cyst",
        ),
        Line2D(
            [],
            [],
            marker="s",
            linestyle="None",
            markerfacecolor="#A58B7E",
            markeredgecolor="#7E665C",
            markersize=5,
            label="MRONJ",
        ),
        Line2D(
            [],
            [],
            marker="D",
            linestyle="None",
            markerfacecolor="#3F3F3F",
            markeredgecolor="white",
            markersize=5.5,
            label="Mean",
        ),
    ]
    right_axis.legend(
        handles=legend_handles,
        loc="lower right",
        frameon=False,
        fontsize=7.5,
        handletextpad=0.4,
        borderaxespad=0.3,
    )

    figure.subplots_adjust(
        left=0.18,
        right=0.985,
        bottom=0.11,
        top=0.90,
        wspace=0.32,
    )

    figure.canvas.draw()
    for tick_label in left_axis.get_yticklabels():
        tick_label.set_horizontalalignment("left")
        tick_label.set_x(-0.37)
    for tick_label in right_axis.get_yticklabels():
        tick_label.set_horizontalalignment("left")
        tick_label.set_x(-0.22)

    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    canvas_bbox = figure.bbox
    for text_artist in figure.findobj(match=plt.Text):
        if not text_artist.get_visible() or not text_artist.get_text().strip():
            continue
        bbox = text_artist.get_window_extent(renderer=renderer)
        if (
            bbox.x0 < canvas_bbox.x0 - 0.5
            or bbox.y0 < canvas_bbox.y0 - 0.5
            or bbox.x1 > canvas_bbox.x1 + 0.5
            or bbox.y1 > canvas_bbox.y1 + 0.5
        ):
            raise RuntimeError(
                f"Text outside figure bounds: {text_artist.get_text()!r}, {bbox}"
            )

    figure.savefig(
        output_path,
        dpi=300,
        facecolor="white",
        bbox_inches="tight",
        pad_inches=0.22,
        pil_kwargs={"compression": "tiff_lzw"},
    )
    figure.savefig(
        output_path.with_suffix(".pdf"),
        format="pdf",
        facecolor="white",
        bbox_inches="tight",
        pad_inches=0.22,
    )
    with Image.open(output_path) as image:
        image.convert("RGB").save(
            output_path,
            format="TIFF",
            compression="tiff_lzw",
            dpi=(300, 300),
        )
    figure.savefig(
        PREVIEW_PATH,
        dpi=180,
        facecolor="white",
        bbox_inches="tight",
        pad_inches=0.22,
    )
    plt.close(figure)


def main():
    if not HPA_DATA.exists():
        raise FileNotFoundError(HPA_DATA)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / "Figure_2.tif"
    build_figure(output_path)
    print(f"Figure 2 rebuilt: {output_path}")


if __name__ == "__main__":
    main()
