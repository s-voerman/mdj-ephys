import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import helper_functions as hf
import pathlib as Path


def read_spikes(data_loc, meta_loc):
    meta = pd.read_csv(os.path.join(meta_loc, "main_ephys_metadata.csv"))
    pdo = pd.DataFrame()
    for idx, row in meta.iterrows():
        if (
            idx % 1 == 0
        ):  # Temp for testing, dont wanna load everything everytime
            current_inj_file = (
                str(row.file)
                + " "
                + str(row.cell)
                + "-"
                + str(row.curr)
                + ".ABF"
            )
            gfp = row.GFP

            # Read data, only positive current injections
            data, sr = hf.read_abf(os.path.join(data_loc, current_inj_file))
            a = len(data) / sr
            data = data = data[int(a / 2) * sr : int(a) * sr]

            # Go through each sweep
            nm_spikes = [0.25, 0.75]
            sweep_size = 1.0  # Assuming 1 trace = 1 second
            sweep = int(sr * sweep_size)
            window = np.arange(0, len(data), sweep)
            for i, x in enumerate(window):
                # Setting up data
                d = data[x : x + sweep]
                normal_data = d[
                    int(nm_spikes[0] * sr) : int(nm_spikes[1] * sr)
                ]

                # Spike index
                spike_index = hf.edge_detection(normal_data, 0)

                # Amount of injected current
                current = i * 10

                # Spike amplitude
                if len(spike_index) > 0:
                    output = {}
                    for c, t in enumerate(spike_index):
                        # Rounds to every 1 ms
                        # Discuss this with Robin
                        base = 3
                        time = base * round(((t / sr) * 1000) / base)
                        isi = spike_index[c] - spike_index[c - 1]
                        isi = (isi * (1000 / sr)) / 1000

                        # Slightly hacky, but this excludes the last spike..
                        if isi > 0:
                            inst = 1 / isi
                            output["file"] = current_inj_file.split(sep=".")[0]
                            output["count"] = time
                            output["amp"] = inst
                            output["current"] = current
                            output["gfp"] = gfp
                            pdo = pd.concat(
                                [pdo, pd.DataFrame([output])],
                                ignore_index=True,
                            )
    return pdo


def instant_rate_plot(data):
    colors = ["green", "black"]
    for curr in data.current.unique():
        fig, ax = plt.subplots(figsize=(1.75, 2))
        plot_data = data.loc[data["current"] == curr]
        sns.lineplot(
            x="count",
            y="amp",
            hue="gfp",
            data=plot_data,
            ax=ax,
            palette=colors,
            lw=0.5,
            legend=False,
            errorbar="sd",
        )
        ax.set_title("Injected current = " + str(curr) + " pA")
        ax.set_ylim(0, 400)
        ax.set_xlim(0, 500)
        ax.set_xticks([0, 250, 500])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_ylabel("Ins. spike rate (Hz)")
        ax.set_xlabel("Time (ms)")
        fig.tight_layout()
        fig.savefig(os.path.join(output_loc, "inst_spr.pdf"))
        plt.show()
    return


def curr_rate(data):
    plot_data = data.groupby(["file", "current", "gfp"]).mean().reset_index()
    fig, ax = plt.subplots(figsize=(2, 2))
    colors = ["green", "black"]
    sns.lineplot(
        x="current",
        y="amp",
        hue="gfp",
        data=plot_data,
        ax=ax,
        palette=colors,
        lw=0.5,
        legend=False,
        errorbar="sd",
    )
    ax.set_ylim(0, 200)
    ax.set_xlim(0, 500)
    ax.set_xticks([0, 250, 500])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylabel("Ins. spike rate (Hz)")
    ax.set_xlabel("Current (pA)")
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, "inst_spr_curr.pdf"))
    plt.show()
    return


if __name__ == "__main__":
    base_folder = Path.cwd()
    data_loc = os.path.join(base_folder, "data")
    meta_loc = os.path.join(base_folder, "metadata")
    output_loc = os.path.join(base_folder, "output")

    data = read_spikes(data_loc, meta_loc)

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
