from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
from PIL import Image


matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

BASE = Path(__file__).resolve().parents[2]
OUTPUT_DIR = BASE / "figures" / "final"

PREVIEW_PATH = OUTPUT_DIR / "figure_3_preview.png"

ROWS = [
    {
        "label": "UKB eBMD (Morris 2019)",
        "display": "UKB eBMD\n(Morris 2019)",
        "beta": 0.01232,
        "se": 0.00328,
        "group": "UKB estimated/heel BMD",
    },
    {
        "label": "UKB BMD (Mbatchou 2021)",
        "display": "UKB BMD\n(Mbatchou 2021)",
        "beta": 0.01177,
        "se": 0.00374,
        "group": "UKB estimated/heel BMD",
    },
    {
        "label": "UKB heel BMD (Loh 2018)",
        "display": "UKB heel BMD\n(Loh 2018)",
        "beta": 0.00487,
        "se": 0.00331,
        "group": "UKB estimated/heel BMD",
    },
    {
        "label": "UKB heel BMD (Elsworth 2018)",
        "display": "UKB heel BMD\n(Elsworth 2018)",
        "beta": -0.00215,
        "se": 0.00457,
        "group": "UKB estimated/heel BMD",
    },
    {
        "label": "Forearm BMD (Zheng 2015)",
        "display": "Forearm BMD\n(Zheng 2015)",
        "beta": -0.18969,
        "se": 0.02779,
        "group": "DXA site BMD",
    },
    {
        "label": "Femoral neck BMD (Zheng 2015)",
        "display": "Femoral neck BMD\n(Zheng 2015)",
        "beta": -0.04444,
        "se": 0.01335,
        "group": "DXA site BMD",
    },
    {
        "label": "Lumbar spine BMD (Zheng 2015)",
        "display": "Lumbar spine BMD\n(Zheng 2015)",
        "beta": -0.02956,
        "se": 0.01553,
        "group": "DXA site BMD",
    },
    {
        "label": "Total body BMD (Medina-Gomez 2018)",
        "display": "Total body BMD\n(Medina-Gomez 2018)",
        "beta": -0.03808,
        "se": 0.00984,
        "group": "DXA site BMD",
    },
    {
        "label": "Ultradistal forearm BMD (Surakka 2020)",
        "display": "Ultradistal forearm BMD\n(Surakka 2020)",
        "beta": 0.00519,
        "se": 0.01887,
        "group": "DXA site BMD",
    },
]

STYLE = {
    "UKB estimated/heel BMD": {
        "color": "#8FA3B2",
        "marker": "o",
    },
    "DXA site BMD": {
        "color": "#B97970",
        "marker": "s",
    },
}


def build_figure(output_path):
    figure, axis = plt.subplots(figsize=(8.4, 8.0), dpi=300)
    figure.patch.set_facecolor("white")
    axis.set_facecolor("white")

    y_positions = list(range(len(ROWS)))[::-1]

    for y_position, row in zip(y_positions, ROWS):
        style = STYLE[row["group"]]
        lower = row["beta"] - 1.96 * row["se"]
        upper = row["beta"] + 1.96 * row["se"]

        axis.plot(
            [lower, upper],
            [y_position, y_position],
            color=style["color"],
            linewidth=1.35,
            solid_capstyle="round",
            zorder=2,
        )
        axis.scatter(
            [row["beta"]],
            [y_position],
            s=49,
            marker=style["marker"],
            color=style["color"],
            edgecolor="#343434",
            linewidth=0.75,
            zorder=3,
        )

    axis.axvline(
        0,
        color="#585858",
        linestyle=(0, (4, 3)),
        linewidth=1.0,
        zorder=1,
    )
    axis.axhline(
        4.5,
        color="#D3D3D3",
        linewidth=0.8,
        zorder=0,
    )

    axis.set_yticks(y_positions, [row["display"] for row in ROWS])
    axis.set_ylim(-0.65, len(ROWS) - 0.25)
    axis.set_xlim(-0.26, 0.10)
    xlabel_artist = figure.text(
        0.5,
        0.135,
        "IVW beta (per SD of RAMP3 NES)",
        ha="center",
        va="bottom",
        fontsize=10.2,
    )
    axis.set_xticks([-0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10])
    axis.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    axis.grid(axis="x", color="#E6E6E6", linewidth=0.8)
    axis.set_axisbelow(True)
    axis.tick_params(axis="y", labelsize=8.2, length=0, pad=3)
    axis.tick_params(axis="x", labelsize=8.8)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color("#8A8A8A")
    axis.spines["bottom"].set_color("#8A8A8A")
    axis.spines["left"].set_linewidth(0.8)
    axis.spines["bottom"].set_linewidth(0.8)

    legend_handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="None",
            markerfacecolor="#8FA3B2",
            markeredgecolor="#343434",
            markeredgewidth=0.75,
            markersize=6,
            label="UKB estimated/heel BMD",
        ),
        Line2D(
            [],
            [],
            marker="s",
            linestyle="None",
            markerfacecolor="#B97970",
            markeredgecolor="#343434",
            markeredgewidth=0.75,
            markersize=6,
            label="DXA site BMD",
        ),
    ]
    axis.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.005),
        ncol=2,
        frameon=False,
        fontsize=8.0,
        handletextpad=0.5,
        columnspacing=1.2,
        borderaxespad=0.2,
    )

    figure.subplots_adjust(
        left=0.30,
        right=0.975,
        bottom=0.20,
        top=0.70,
    )

    figure.canvas.draw()
    tick_labels = list(axis.get_yticklabels())
    for tick_label in tick_labels:
        tick_label.set_horizontalalignment("left")
        tick_label.set_x(0)

    figure.canvas.draw()
    tick_renderer = figure.canvas.get_renderer()
    longest_tick_width = max(
        tick_label.get_window_extent(renderer=tick_renderer).width
        for tick_label in tick_labels
    )
    tick_gap_pixels = 10
    tick_offset = (longest_tick_width + tick_gap_pixels) / axis.bbox.width
    for tick_label in tick_labels:
        tick_label.set_x(-tick_offset)

    xlabel_artist.set_visible(False)
    figure.canvas.draw()
    base_bbox = figure.get_tightbbox(figure.canvas.get_renderer())
    visual_center_x = (
        (base_bbox.x0 + base_bbox.x1) / 2
    ) / figure.get_size_inches()[0]
    xlabel_artist.set_visible(True)
    xlabel_artist.set_x(visual_center_x)

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

    xlabel_bbox = xlabel_artist.get_window_extent(renderer=renderer)
    x_tick_bboxes = [
        tick_label.get_window_extent(renderer=renderer)
        for tick_label in axis.get_xticklabels()
        if tick_label.get_text().strip()
    ]

    if x_tick_bboxes and xlabel_bbox.y1 >= min(
        bbox.y0 for bbox in x_tick_bboxes
    ) - 2:
        raise RuntimeError(
            "Bottom x-axis label overlaps the tick labels: "
            f"xlabel={xlabel_bbox}, ticks={x_tick_bboxes}"
        )

    final_tight_bbox = figure.get_tightbbox(renderer)
    final_center_x = (
        (final_tight_bbox.x0 + final_tight_bbox.x1) / 2
    ) * figure.dpi
    xlabel_center_x = (xlabel_bbox.x0 + xlabel_bbox.x1) / 2
    if abs(xlabel_center_x - final_center_x) > 2:
        raise RuntimeError(
            "Bottom x-axis label is not visually centered: "
            f"final_center={final_center_x}, "
            f"xlabel_center={xlabel_center_x}"
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
    output_path = OUTPUT_DIR / "Figure_3.tif"
    build_figure(output_path)
    print(f"Figure 3 rebuilt: {output_path}")


if __name__ == "__main__":
    main()
