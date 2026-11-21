import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# This is for positive current injections

def spikes(spike_data, output_loc):
    plot_data = pd.DataFrame()
    for i, v in enumerate(spike_data["Cell"].unique()):
        temp_df = spike_data.loc[spike_data["Cell"] == v]
        temp_df.loc[:, "Max_current"] = (len(temp_df) * 5) - 5
        plot_data = pd.concat([plot_data, temp_df])
    plot_data = plot_data.loc[plot_data.Current > 0]
    fig, ax = plt.subplots(1, 3, figsize=(3.5, 1.75))
    colors = ["green", "black"]
    sns.lineplot(
        x="Current",
        y="TTF",
        hue="GFP",
        data=plot_data,
        ax=ax[0],
        errorbar="sd",
        lw=0.5,
        palette=colors,
    )
    sns.lineplot(
        x="Current",
        y="mAHP",
        hue="GFP",
        data=plot_data,
        ax=ax[1],
        errorbar="sd",
        lw=0.5,
        palette=colors,
    )
    sns.lineplot(
        x="Current",
        y="Firing rate",
        hue="GFP",
        data=plot_data,
        ax=ax[2],
        errorbar="sd",
        lw=0.5,
        palette=colors,
    )

    # axes stuff
    for a in ax.flatten():
        a.xaxis.set_inverted(True)
        a.set_xlim(0, 501)
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
        a.tick_params(axis="x", labelrotation=45)
        a.set_xlabel("Current (pA)")
        a.set_xticks([0, 250, 500])
        # a.axvline(x=100, lw=0.5, color="black", ls="--")
        # a.axvline(x=250, lw=0.5, color="black", ls="--")
        # a.axvline(x=500, lw=0.5, color="black", ls="--")
        a.get_legend().remove()

    ax[0].set_ylabel("Spike latency (ms)")
    ax[1].set_ylabel("mAHP (mV)")
    ax[2].set_ylabel("Firing rate (Hz)")
    ax[0].set_ylim(0, 300)
    ax[1].set_ylim(0, 15)
    ax[2].set_ylim(0, 200)
    ax[0].set_yticks([0, 75, 150, 225, 300])

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.8, wspace=0.8)
    fig.savefig(os.path.join(output_loc, "ramp_plot_cc.pdf"))
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

    # Folders
    base_folder = Path.cwd()
    data_loc = os.path.join(base_folder, "data")
    meta_loc = os.path.join(base_folder, "metadata")
    output_loc = os.path.join(base_folder, "output")
    spike_data = pd.read_csv(os.path.join(output_loc, "spike_properties_filtered.csv"))

    spikes(spike_data, output_loc)
