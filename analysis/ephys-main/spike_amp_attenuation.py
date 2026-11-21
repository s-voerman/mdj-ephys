import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import helper_functions as hf
from pathlib import Path


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
            data = data[len(data) - sr : len(data)]

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
                baseline = [0.1, 0.24]
                base = np.mean(
                    d[int(baseline[0] * sr) : int(baseline[1] * sr)]
                )
                if len(spike_index) > 0:
                    output = {}
                    for c, t in enumerate(spike_index):
                        # Rounds to every x ms
                        b = 5
                        time = b * round(((t / sr) * 1000) / b)
                        spike_data = normal_data[
                            t - 10 : t + 50
                        ]  # maybe change t-10 to 0
                        peak = np.max(spike_data)
                        spike_amplitude = peak - base
                        output["file"] = current_inj_file.split(sep=".")[0]
                        output["count"] = time
                        output["amp"] = spike_amplitude
                        output["current"] = current
                        output["gfp"] = gfp
                        pdo = pd.concat(
                            [pdo, pd.DataFrame([output])], ignore_index=True
                        )
    pdo = pdo.groupby(["count", "file", "gfp"]).mean().reset_index()
    pdo = pdo.loc[pdo['count'] <301]
    pdo.to_excel(os.path.join(output_loc, "spikeattenuation.xlsx"))
    pdo.to_csv(os.path.join(output_loc, "spikeattenuation.csv"))
    return pdo


def attenuation_plots(data):
    colors = ["green", "black"]
    for curr in data.current.unique():
        fig, ax = plt.subplots(figsize=(1.75, 1.75))
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
        ax.set_ylim(80, 140)
        ax.set_xlim(0, 300)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_ylabel("AP amplitude (mV)")
        ax.set_xlabel("Time (ms)")
        fig.tight_layout()
        fig.savefig(os.path.join(output_loc, "attenuation.pdf"))
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
    attenuation_plots(data)
