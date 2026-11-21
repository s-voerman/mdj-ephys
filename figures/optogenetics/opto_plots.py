import os

import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path


def opto_plot(data, output_loc):
    tags = ["sep_vc", "sep_cc", "adjust_vc", "adjust_cc"]
    for i, t in enumerate(tags):
        plot_data = data.loc[data["tag"] == t]
        fig, ax = plt.subplots(figsize=(1.5, 1.25))
        ax.boxplot(
            x=plot_data.chrimson_amp,
            positions=[0.8],
            showfliers=False,
            patch_artist=True,
            boxprops=dict(linewidth=0.5, facecolor="white"),
            whiskerprops=dict(linewidth=0.5),
            capprops=dict(linewidth=0.5),
            medianprops=dict(linewidth=2, color="black"),
        )

        ax.boxplot(
            x=plot_data.chronos_amp,
            positions=[0.2],
            showfliers=False,
            patch_artist=True,
            boxprops=dict(linewidth=0.5, facecolor="white"),
            whiskerprops=dict(linewidth=0.5),
            capprops=dict(linewidth=0.5),
            medianprops=dict(linewidth=2, color="black"),
        )
        for index, row in plot_data.iterrows():
            points = [row.chronos_amp, row.chrimson_amp]
            ax.plot([0.2, 0.8], points, color="black", linewidth=0.5)
            ax.plot(
                [0.2, 0.8],
                points,
                "d",
                color="white",
                markeredgecolor="black",
                mew=0.5,
                markersize=5,
            )

            # ax.get_xaxis().set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

        ax.set_xticks([0.2, 0.8])
        ax.set_xlim(0, 1)
        ax.set_xticklabels(["B", "R"])
        if t == "sep_vc":
            ax.set_ylim(-500, 0)
            ax.set_yticks([-500, -250, 0])
            ax.set_ylabel("Amplitude (pA)")
            # ax.axhline(y=0, color="black", lw=0.5, alpha=0.8, ls="--")
        elif t == "sep_cc":
            ax.set_ylim(0, 40)
            ax.set_yticks([0, 20, 40])
            ax.set_ylabel("Amplitude (mV)")
        elif t == "adjust_vc":
            ax.set_ylim(-500, 0)
            ax.set_yticks([-500, -250, 0])
            ax.set_ylabel("Amplitude (pA)")
            # ax.axhline(y=0, color="black", lw=0.5, alpha=0.8, ls="--")
        elif t == "adjust_cc":
            ax.set_ylim(0, 30)
            ax.set_yticks([0, 15, 30])
            ax.set_ylabel("Amplitude (mV)")
        fig.tight_layout()
        fig.savefig(os.path.join(output_loc, str(t) + "plot.pdf"))
        plt.show()
    return


if __name__ == "__main__":
    # Plotting params
    plt.rcParams["font.sans-serif"] = "Calibri"
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["axes.linewidth"] = 0.5
    plt.rcParams["xtick.major.width"] = 0.5
    plt.rcParams["xtick.minor.width"] = 0.5
    plt.rcParams["ytick.major.width"] = 0.5
    plt.rcParams["ytick.minor.width"] = 0.5
    fonts = {
        "font.size": 8,
        "axes.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.labelsize": 8,
    }
    plt.rcParams.update(fonts)

    base_folder = Path.cwd()
    data_loc = os.path.join(base_folder, "data")
    meta_loc = os.path.join(base_folder, "metadata")
    output_loc = os.path.join(base_folder, "output")

    good_cells = [
        "250527 3",
        "250603 6",
        "250612 2",
        "250612 3",
        "250612 6",
        "250612 7",
        "250612 10",
        "250624 5",
        "250624 7",
    ]

    data = pd.read_excel(os.path.join(output_loc, "optogenetics_data.csv"))
    r_data = data.loc[data.cell.isin(good_cells)]
    opto_plot(r_data, output_loc)
