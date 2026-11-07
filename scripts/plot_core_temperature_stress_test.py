"""Plot the coherent core-heating sensitivity diagnostic with Python.

Figure contract: a quantitative grid tests the two regional mass-fraction
trends as core heating increases. The left panel tests concentration relative
to clouds; the right tests concentration within fragments. Curves are
deterministic scenarios, not confidence bands. All 231 sources are retained.
"""
from pathlib import Path
import argparse
import sys

from astropy.table import Table
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path(
        "results/paper_ii_referee_response/diagnostics/R8d/2026-11-01"))
    parser.add_argument("--qa-scripts", type=Path, required=True,
                        help="Directory containing audit_panel_alignment.py")
    args = parser.parse_args()
    sys.path.insert(0, str(args.qa_scripts))
    from audit_panel_alignment import require_matplotlib_panel_alignment

    data = Table.read(args.input / "regional_temperature_grid.ecsv")
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 10,
        "pdf.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": .8,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False})
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.85))
    fig.subplots_adjust(left=.09, right=.99, bottom=.20, top=.97, wspace=.38)
    metrics = ["core_cloud_percent", "core_fragment_percent"]
    colors = ["#007C91", "#D77B25", "#76549F"]
    for ax, metric, label, ymax in zip(axes, metrics,
            ["Core / cloud mass (%)", "Core / fragment mass (%)"], [3.5, 60]):
        for code, color in zip(["Y", "I", "E"], colors):
            rows = data[data["code"] == code]
            ax.plot(rows["h"], rows[metric], color=color, linewidth=1.8,
                    label=f'{rows[0]["region"]} ({rows[0]["endpoint_k"]:g} K)')
            ax.scatter([0], [rows[metric][0]], s=22, facecolors="white",
                       edgecolors=color, linewidths=1.1, zorder=3)
            ax.scatter([1], [rows[metric][-1]], s=22, color=color, zorder=3)
        ax.set(xlim=(-.035, 1.035), ylim=(0, ymax), ylabel=label,
               xlabel="Core-heating fraction, h")
        ax.set_xticks([0, .25, .5, .75, 1])
        ax.tick_params(direction="in", length=3)
        ax.tick_params(axis="x", pad=7)
    for ax, letter in zip(axes, ["(a)", "(b)"]):
        ax.text(.025, .98, letter, transform=ax.transAxes, ha="left", va="top")
    axes[1].legend(loc="upper right", fontsize=8.5, handlelength=1.5)
    fig.canvas.draw()
    stem = args.input / "temperature_stress_test"
    require_matplotlib_panel_alignment(fig, json_out=f"{stem}.alignment.json",
        overlay_svg=f"{stem}.alignment.svg", tolerance_pt=1.5, strict=True)
    fig.savefig(f"{stem}.pdf", bbox_inches="tight", pad_inches=.025)
    fig.savefig(f"{stem}.svg", bbox_inches="tight", pad_inches=.025)
    fig.savefig(f"{stem}.png", dpi=600, bbox_inches="tight", pad_inches=.025)
    fig.savefig(f"{stem}.tiff", dpi=600, bbox_inches="tight", pad_inches=.025)
    plt.close(fig)


if __name__ == "__main__":
    main()
