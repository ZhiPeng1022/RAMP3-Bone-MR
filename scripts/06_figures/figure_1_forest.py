from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from PIL import Image


matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

BASE = Path(__file__).resolve().parents[2]
OUTPUT_DIR = BASE / "figures" / "final"

PREVIEW_PATH = OUTPUT_DIR / "figure_1_preview.png"

ROWS = [
    {
        "label": "eBMD (SD)",
        "beta": 0.01232,
        "se": 0.00328,
        "p_text": "P=1.7×10⁻⁴",
        "color": "#718EA4",
        "marker": "o",
    },
    {
        "label": "Fracture (log OR)",
        "beta": 0.00230,
        "se": 0.00460,
        "p_text": "P=0.61",
        "color": "#A18F7D",
        "marker": "D",
    },
    {
        "label": "Serum calcium\n(mmol/L)",
        "beta": 0.00233,
        "se": 0.00109,
        "p_text": "P=0.033",
        "color": "#7E9678",
        "marker": "s",
    },
    {
        "label": "Serum phosphate\n(mmol/L)",
        "beta": -0.00224,
        "se": 0.00109,
        "p_text": "P=0.041",
        "color": "#B08A86",
        "marker": "^",
    },
]


def ci_text(beta, se):
    lower = beta - 1.96 * se
    upper = beta + 1.96 * se
    return f"{beta:.4f} ({lower:.4f}, {upper:.4f})"


def build_figure(output_path):
    figure, axis = plt.subplots(figsize=(9.7, 4.8), dpi=300)
    figure.patch.set_facecolor("white")
    axis.set_facecolor("white")

    y_positions = list(range(len(ROWS)))[::-1]
    column_pairs = []

    for y_position, row in zip(y_positions, ROWS):
        beta = row["beta"]
        se = row["se"]
        lower = beta - 1.96 * se
        upper = beta + 1.96 * se

        axis.plot(
            [lower, upper],
            [y_position, y_position],
            color="#4C4C4C",
            linewidth=1.35,
            solid_capstyle="round",
            zorder=2,
        )
        axis.scatter(
            [beta],
            [y_position],
            s=58,
            marker=row["marker"],
            color=row["color"],
            edgecolor="#2F2F2F",
            linewidth=0.8,
            zorder=3,
        )

        ci_artist = axis.text(
            1.04,
            y_position,
            ci_text(beta, se),
            transform=axis.get_yaxis_transform(),
            ha="left",
            va="center",
            fontsize=8.2,
            color="#222222",
            clip_on=False,
        )
        p_artist = axis.text(
            1.43,
            y_position,
            row["p_text"],
            transform=axis.get_yaxis_transform(),
            ha="left",
            va="center",
            fontsize=8.8,
            color="#222222",
            clip_on=False,
        )
        column_pairs.append((ci_artist, p_artist))

    axis.axvline(
        0,
        color="#5A5A5A",
        linestyle=(0, (4, 3)),
        linewidth=1.0,
        zorder=1,
    )

    axis.set_yticks(y_positions, [row["label"] for row in ROWS])
    axis.set_ylim(-0.6, len(ROWS) - 0.25)
    axis.set_xlim(-0.009, 0.0205)
    axis.set_xlabel(
        "Effect size per 1-SD increase in genetically proxied RAMP3 expression",
        fontsize=10.2,
        labelpad=8,
    )
    axis.grid(axis="y", color="#E6E6E6", linewidth=0.8, zorder=0)
    axis.tick_params(axis="y", labelsize=10.2, length=0, pad=8)
    axis.tick_params(axis="x", labelsize=9.5, length=3)

    axis.set_xticks([-0.005, 0.000, 0.005, 0.010, 0.015, 0.020])
    axis.xaxis.set_major_formatter(FormatStrFormatter("%.3f"))

    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color("#8A8A8A")
    axis.spines["bottom"].set_color("#8A8A8A")
    axis.spines["left"].set_linewidth(0.8)
    axis.spines["bottom"].set_linewidth(0.8)

    axis.text(
        1.04,
        len(ROWS) - 0.48,
        r"$\beta$ (95% CI)",
        transform=axis.get_yaxis_transform(),
        ha="left",
        va="center",
        fontsize=8.7,
        fontweight="bold",
        color="#222222",
        clip_on=False,
    )
    axis.text(
        1.43,
        len(ROWS) - 0.48,
        "P",
        transform=axis.get_yaxis_transform(),
        ha="left",
        va="center",
        fontsize=8.8,
        fontweight="bold",
        color="#222222",
        clip_on=False,
    )

    note_artist = figure.text(
        0.5,
        0.025,
        (
            "Units differ across outcomes (eBMD: SD; fracture: log OR; "
            "serum calcium and phosphate: mmol/L); comparisons should be made "
            "within each row."
        ),
        ha="center",
        va="bottom",
        fontsize=7.8,
        color="#444444",
    )

    figure.subplots_adjust(
        left=0.25,
        right=0.72,
        bottom=0.18,
        top=0.93,
    )

    figure.canvas.draw()
    for tick_label in axis.get_yticklabels():
        tick_label.set_horizontalalignment("left")
        tick_label.set_x(-0.29)

    note_artist.set_visible(False)
    figure.canvas.draw()
    base_bbox = figure.get_tightbbox(figure.canvas.get_renderer())
    visual_center_pixels = (
        (base_bbox.x0 + base_bbox.x1) / 2
    ) * figure.dpi
    note_artist.set_x(visual_center_pixels / figure.bbox.width)
    note_artist.set_visible(True)

    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    canvas_bbox = figure.bbox

    for text_artist in figure.findobj(match=plt.Text):
        if not text_artist.get_visible():
            continue
        if not text_artist.get_text().strip():
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

    for ci_artist, p_artist in column_pairs:
        ci_bbox = ci_artist.get_window_extent(renderer=renderer)
        p_bbox = p_artist.get_window_extent(renderer=renderer)
        if ci_bbox.x1 + 2 > p_bbox.x0:
            raise RuntimeError(
                "CI and P-value columns overlap: "
                f"{ci_artist.get_text()!r}, {p_artist.get_text()!r}"
            )

    final_tight_bbox = figure.get_tightbbox(renderer)
    final_center_pixels = (
        (final_tight_bbox.x0 + final_tight_bbox.x1) / 2
    ) * figure.dpi
    note_bbox = note_artist.get_window_extent(renderer=renderer)
    note_center_pixels = (note_bbox.x0 + note_bbox.x1) / 2
    if abs(note_center_pixels - final_center_pixels) > 2:
        raise RuntimeError(
            "Figure 1 note is not visually centered: "
            f"final_center={final_center_pixels}, "
            f"note_center={note_center_pixels}"
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
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / "Figure_1.tif"
    build_figure(output_path)
    print(f"Figure 1 rebuilt: {output_path}")


if __name__ == "__main__":
    main()
