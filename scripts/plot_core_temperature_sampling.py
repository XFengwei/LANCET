"""Replot the Figure 9 mass hierarchy under two temperature priors.

Contract: compare the sampled mass fractions with the adopted Tenv baseline;
two equal panels isolate the prior choice. Error bars describe temperature
sampling only. A separate export tests regionally shared heating. Python,
183-mm quantitative grid; all regional sums use all 231 adopted sources.
"""
from pathlib import Path
import argparse
import sys

from astropy.table import Table
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from lancet.core_formation import draw_temperature_efficiency_cascade as draw_trial_cascade


COLORS = ["#246FA8", "#D4851E", "#BD3654"]
MARKERS = ["o", "s", "^"]



def figure9_trials(table, out, require_alignment):
    """Two alternative arrangements; preserve the manuscript's existing PDF."""
    root = Path(__file__).resolve().parents[1]
    original = Table.read(root/"results/core_formation/tables/core_formation_summary.ecsv")
    out.mkdir(parents=True, exist_ok=True)
    for layout in ("overlay", "comparison"):
        if layout == "overlay":
            fig, ax = plt.subplots(figsize=(3.9, 3.8))
            fig.subplots_adjust(left=.17, right=.98, bottom=.21, top=.98)
            draw_trial_cascade(ax, table, original)
        else:
            fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6))
            fig.subplots_adjust(left=.085, right=.985, bottom=.22, top=.98, wspace=.36)
            draw_trial_cascade(axes[0], table, original)
            ax = axes[1]
            for j, (base, color, marker) in enumerate(zip(original, COLORS, MARKERS)):
                offset = (j-1)*.13
                y = [100*base["fragment_to_core_efficiency"]]
                lows, highs = [0.], [0.]
                for prior in ("uniform", "brightness"):
                    r = table[(table["coupling"]=="independent") &
                              (table["prior"]==prior) & (table["code"]==base["subregion_code"])][0]
                    p16, p50, p84 = [r[f"core_fragment_{p}_percent"] for p in ("p16", "p50", "p84")]
                    y.append(p50)
                    lows.append(p50-p16)
                    highs.append(p84-p50)
                ax.plot(np.arange(3)+offset, y, color=color, lw=.8, alpha=.45)
                ax.scatter([offset], y[:1], marker=marker, s=30, facecolors="white", edgecolors=color)
                ax.errorbar(np.array([1.,2.])+offset, y[1:], yerr=[lows[1:], highs[1:]],
                    fmt=marker, ms=5, color=color, mfc=color, mec="white", mew=.4,
                    lw=1.2, capsize=3)
            ax.set(xlim=(-.4, 2.4), ylim=(0, 61), ylabel="Core / fragment mass (%)")
            ax.set_xticks([0, 1, 2], ["Adopted\nTenv", "Uniform", "Brightness\ndependent"])
            ax.set_xlabel("Core temperature model")
            ax.tick_params(direction="in", pad=5)
            for a, letter in zip(axes, ["(a)", "(b)"]):
                a.text(.035, .975, letter, transform=a.transAxes, va="top", ha="left")
        stem = out/f"Figure9_trial_{layout}"
        fig.canvas.draw()
        require_alignment(fig, json_out=f"{stem}.alignment.json",
            overlay_svg=f"{stem}.alignment.svg", tolerance_pt=1.5, strict=True)
        fig.savefig(f"{stem}.pdf", bbox_inches="tight", pad_inches=.025)
        fig.savefig(f"{stem}.svg", bbox_inches="tight", pad_inches=.025)
        fig.savefig(f"{stem}.png", dpi=600, bbox_inches="tight", pad_inches=.025)
        fig.savefig(f"{stem}.tiff", dpi=600, bbox_inches="tight", pad_inches=.025)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path(
        "results/paper_ii_referee_response/diagnostics/R8d/2026-11-01/temperature_sampling"))
    parser.add_argument("--qa-scripts", type=Path, required=True)
    parser.add_argument("--trials", action="store_true",
                        help="Make Figure 9 integration trials instead of the prior-comparison figures")
    args = parser.parse_args()
    sys.path.insert(0, str(args.qa_scripts))
    from audit_panel_alignment import require_matplotlib_panel_alignment

    table = Table.read(args.input/"regional_sampling_summary.ecsv")
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 10,
        "pdf.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": .8,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False})
    if args.trials:
        figure9_trials(table, args.input/"figure9_trials", require_matplotlib_panel_alignment)
        return
    colors = ["#007C91", "#D77B25", "#76549F"]
    shapes = ["o", "s", "^"]
    scales = np.array([4., .15, .01])
    assert np.all(scales > 0)
    for coupling in ("independent", "regional_shared"):
        fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True)
        fig.subplots_adjust(left=.085, right=.99, bottom=.21, top=.98, wspace=.14)
        for ax, prior, title, letter in zip(axes, ["uniform", "brightness"],
                ["Uniform sampling", "Brightness-dependent"], ["(a)", "(b)"]):
            selected = table[(table["coupling"]==coupling)&(table["prior"]==prior)]
            for row, color, marker in zip(selected, colors, shapes):
                fragment = row["fragment_mass_msun"]/row["cloud_mass_msun"]
                low, middle, high = [row[f"core_cloud_{p}_percent"]/100 for p in ("p16", "p50", "p84")]
                baseline = row["baseline_mass_msun"]/row["cloud_mass_msun"]
                fractions = np.array([1., fragment, middle])
                assert low > 0 and high >= middle >= low and np.all(fractions > 0)
                ax.plot(scales, fractions, color=color, marker=marker, ms=5, lw=1.5,
                        label=str(row["region"]))
                ax.errorbar([scales[-1]], [middle], yerr=[[middle-low], [high-middle]],
                            fmt="none", color=color, lw=1.3, capsize=4, zorder=4)
                ax.plot(scales[1:], [fragment, baseline], color=color,
                        linestyle="--", alpha=.55, lw=.9)
                ax.scatter([scales[-1]], [baseline], marker=marker, s=32,
                           facecolors="white", edgecolors=color, lw=1, zorder=4)
            ax.set(xscale="log", yscale="log", xlim=(6., .005), ylim=(.0006, 3.))
            ax.set_xticks(scales, ["Cloud\n4 pc", "Fragment\n0.15 pc", "Core\n0.01 pc"])
            ax.set_yticks([.001, .01, .1, 1], ["0.1%", "1%", "10%", "100%"])
            ax.tick_params(direction="in", axis="both", pad=5)
            ax.text(.025, .975, letter, transform=ax.transAxes, ha="left", va="top")
            ax.text(.98, .975, title, fontsize=9, transform=ax.transAxes, ha="right", va="top")
            ax.set_xlabel("Physical scale")
        axes[0].set_ylabel("Mass / cloud mass")
        axes[0].legend(loc="lower left", fontsize=9, handlelength=1.8)
        handles = [Line2D([], [], color="0.3", marker="o", mfc="white", ls="--", lw=.9, label="Adopted Tenv"),
                   Line2D([], [], color="0.3", marker="o", ls="-", lw=1.5, label="Sampled median")]
        axes[1].legend(handles=handles, loc="lower left", fontsize=8.5)
        stem = args.input/f"mass_efficiency_{coupling}"
        fig.canvas.draw()
        require_matplotlib_panel_alignment(fig, json_out=f"{stem}.alignment.json",
            overlay_svg=f"{stem}.alignment.svg", tolerance_pt=1.5, strict=True)
        fig.savefig(f"{stem}.pdf", bbox_inches="tight", pad_inches=.025)
        fig.savefig(f"{stem}.svg", bbox_inches="tight", pad_inches=.025)
        fig.savefig(f"{stem}.png", dpi=600, bbox_inches="tight", pad_inches=.025)
        fig.savefig(f"{stem}.tiff", dpi=600, bbox_inches="tight", pad_inches=.025)
        plt.close(fig)


if __name__ == "__main__":
    main()
