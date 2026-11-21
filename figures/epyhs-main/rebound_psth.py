"""
Basic idea:
- Read spike data from hyperpolarizations, probably just do all negative currents from -500 -> -10
- Detect rebound spikes after release (>750 ms after trace start)
- Get index of spike and convert to time, allocate val = 1 if there is a spike at an index, 0 if not
- Plot 0s and 1s across time for every cell, and as an average
- For averaging maybe assign spikes to specific bins of time of ~20-50 ms
- Split into GFP+ and GFP- cells
"""

import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import helper_functions as hf

from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path


# Maybe normalize the number of spikes to the difference between baseline and minimum voltage after negative current injection?
# That way cells might be a bit more comparable


# Maybe make time bins a bit smaller? IDK really


def psth_analysis(meta_loc, data_loc):
    meta = pd.read_excel(os.path.join(meta_loc, "MDJ_meta_V2_f.xlsx"))

    bin_size = 50  # this is 1 ms
    aspect = len(meta.loc[meta.GFP == "pos"]) / len(
        meta.loc[meta.GFP == "neg"]
    )
    fig, ax = plt.subplots(
        2, 1, figsize=(2.25, 2), gridspec_kw={"height_ratios": [aspect, 1]}
    )
    cmap0 = LinearSegmentedColormap.from_list("", ["white", "darkblue"])
    for pl, gfp in enumerate(meta.GFP.unique()):
        sp_meta = meta.loc[meta.GFP == gfp]
        plot_data = np.array([]).reshape(0, 1000)
        for index, row in sp_meta.iterrows():
            # Construct filename
            file = os.path.join(
                data_loc,
                str(row.file)
                + " "
                + str(row.cell)
                + "-"
                + str(row.curr)
                + ".ABF",
            )

            # Read data, detect spikes
            data, sr = hf.read_abf_v2(file)
            spike_data = data[0:sr]
            spike_index = hf.edge_detection(spike_data, 0)

            # Setting up data for plotting
            bins = np.arange(0, sr, bin_size)
            output = np.array([])
            for i, x in enumerate(bins):
                bin = [x, x + bin_size]
                n, m = min(bin), max(bin)
                c = [n < f < m for f in spike_index]
                spike_count = sum(c)
                output = np.append(output, int(spike_count))
            plot_data = np.append(plot_data, [output], axis=0)
        # plot_data = np.sort(plot_data, axis=0)
        # np.savetxt(os.path.join(output_loc, str(gfp) + ".csv"), plot_data)

        # Actually plot stuff
        ax[pl].imshow(plot_data, cmap=cmap0, aspect="auto")

    # New limits with 1 ms bins
    for a in ax.flatten():
        a.set_xlim(699, 999)
        a.axvline(x=749, color="black", ls="-", lw=0.5)
        a.set_xticks([699, 749, 799, 849, 899, 949, 999])
        a.set_xticklabels(["-50", "0", "50", "100", "150", "200", "250"])
        a.set_ylabel("Cell #")
        a.set_xlabel("Latency to spike (ms)")
        a.tick_params(axis="x", labelrotation=45)
    # Removing some unnecessary stuff in ax[0]
    ax[0].set_xlabel(" ")
    ax[0].get_xaxis().set_visible(False)
    ax[0].set_ylim(len(meta.loc[meta.GFP == "pos"]), 0)
    ax[1].set_ylim(len(meta.loc[meta.GFP == "neg"]), 0)
    fig.tight_layout()
    fig.savefig(
        os.path.join(output_loc, "rebound_spiking_rates_negg.pdf"),
        dpi=1000,
    )
    plt.show()
    return plot_data


def psth_allcells(meta_loc, data_loc):
    # Works pretty well now :)
    meta = pd.read_csv(os.path.join(meta_loc, "main_ephys_metadata.csv"))
    plot_data = np.zeros(100)
    for index, row in meta.iterrows():
        file = os.path.join(
            data_loc,
            str(row.file) + " " + str(row.cell) + "-" + str(row.curr) + ".ABF",
        )
        GFP = row.GFP
        if GFP == "neg":
            continue
        # Read data, detect spikes
        data, sr = hf.read_abf_v2(file)
        spike_data = data[
            0:sr
        ]  # only first trace with max negative current injection
        spike_index = hf.edge_detection(spike_data, 0)
        spike_index = spike_index / 500
        for s in spike_index:
            plot_data[int(s)] += 1
    fig, ax = plt.subplots()
    ax.bar(x=range(len(plot_data)), height=plot_data)
    ax.axvline(x=74, lw=0.5, ls="--", color="black")
    # Do like a kdeplot overtop of this?
    fig.tight_layout()
    plt.show()
    return fig, ax


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
