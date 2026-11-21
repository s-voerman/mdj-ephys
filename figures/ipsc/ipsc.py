import os

import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path


def read_data(meta, clean_loc, output_loc):
    combined_data = pd.DataFrame()
    for index, row in meta.iterrows():
        cell_data = pd.DataFrame()
        traces = [int(s) for s in row.traces.split(",")]
        trace_files = [
            os.path.join(
                clean_loc,
                str(row.file)
                + " "
                + str(row.cell)
                + "-"
                + str(f)
                + "_cleaned.xlsx",
            )
            for f in traces
        ]
        for time, files in enumerate(trace_files):
            data = pd.read_excel(files)
            data["trace"] = time
            data["cell_id"] = str(row.file) + " " + str(row.cell)
            data["frequency"] = len(data) / 120
            cell_data = pd.concat([cell_data, data])

        # data normalization
        # probably normalize per trace, where baseline = 1
        vars = ["amplitude", "frequency", "rise", "tau"]
        for i, v in enumerate(vars):
            cell_data.loc[:, "norm_" + v] = cell_data.loc[:, v] / np.mean(
                cell_data.loc[(cell_data.trace == 0), v]
            )
        combined_data = pd.concat([combined_data, cell_data])
    print(combined_data["cell_id"])
    return combined_data


def normalized_plot(data, img_output):
    data = (
        data.groupby(["trace", "cell_id"])
        .mean(numeric_only=True)
        .reset_index()
    )
    data.to_csv(os.path.join(output_loc, "ipsc_ptx.csv"))
    data.to_excel(os.path.join(output_loc, "ipsc_ptx.xlsx"))
    vars = ["norm_frequency", "norm_amplitude", "norm_rise", "norm_tau"]
    fig, ax = plt.subplots(1, len(vars), figsize=(3.75, 1.5))
    for a in ax.flatten():
        a.axvline(x=1, color="red", lw=1, ls="--")
    for i, v in enumerate(vars):
        sns.lineplot(
            x="trace",
            y=v,
            data=data,
            ax=ax[i],
            errorbar="sd",
            err_style="bars",
            lw=0.5,
            color="black",
            err_kws={"elinewidth": 0.5},
        )
    for a in ax.flatten():
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
        a.set_xlabel("Time (minutes)")
        a.set_xticks([0, 5, 10])
        a.set_xticklabels([0, 10, 20])
        a.set_xlabel(" ")
    # axis limits
    ax[0].set_ylim(0, 2)
    ax[1].set_ylim(0, 2)
    ax[2].set_ylim(0, 2)
    ax[3].set_ylim(0, 2)
    # labels
    ax[0].set_ylabel("Frequency")
    ax[1].set_ylabel("Amplitude")
    ax[2].set_ylabel("Rise")
    ax[3].set_ylabel("Decay")
    # ax[3].set_xlabel("Time (minutes)")

    fig.tight_layout()
    fig.savefig(os.path.join(img_output, "ipsc_ptx_norm.pdf"))
    plt.show()
    return

def ipscs(output_loc):
    # Read data from excel files
    all_data = pd.DataFrame()
    files = os.listdir(os.path.join(output_loc, "sIPSC/cleaned"))
    for f in files:
        data = pd.read_excel(os.path.join(output_loc, "sIPSC/cleaned/" + f))
        data["frequency"] = len(data) / 120
        all_data = pd.concat([all_data, data])
    all_data = all_data.groupby(["cell"]).mean(numeric_only=True).reset_index()
    all_data["tau"] = all_data["tau"] * 1000
    all_data["rise"] = all_data["rise"] * 1000
    all_data["amplitude"] = all_data["amplitude"] * -1

    # Plotting
    fig, ax = plt.subplots(1, 4, figsize=(3.75, 1))
    vars = ["frequency", "amplitude", "rise", "tau"]
    for i, v in enumerate(vars):
        sns.boxplot(
            y=v,
            data=all_data,
            width=0.5,
            ax=ax[i],
            fliersize=0,
            color="white",
            linecolor="black",
            linewidth=0.5,
        )
        sns.stripplot(
            y=v, data=all_data, jitter=0, color="black", ax=ax[i], size=3, marker="d"
        )
    # Make things pretty
    for a in ax.flatten():
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
        a.spines["bottom"].set_visible(False)
        a.xaxis.set_visible(False)
    ax[0].set_ylim(0, 10)
    ax[1].set_ylim(0, 150)
    ax[2].set_ylim(0, 4)
    ax[3].set_ylim(0, 8)
    ax[1].set_yticks([0, 75, 150])
    ax[3].set_yticks([0, 4, 8])
    ax[0].set_ylabel("Frequency (Hz)")
    ax[1].set_ylabel("Amplitude (pA)")
    ax[2].set_ylabel("Rise (ms)")
    ax[3].set_ylabel("Decay (ms)")
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, "ispc_freq.pdf"))
    plt.show()
    return fig, ax


if __name__ == "__main__":
    base_folder = Path.cwd()
    data_loc = os.path.join(base_folder, "data")
    meta_loc = os.path.join(base_folder, "metdata")
    output_loc = os.path.join(base_folder, "output")
    img_output = os.path.join(base_folder, "output")
    clean_loc = os.path.join(output_loc, "data")
    meta = pd.read_csv(os.path.join(meta_loc, "ipsc_ptx.csv"))
    data = read_data(meta, clean_loc, output_loc)

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

    # Global color setting
    sns.set_palette("icefire", 2, desat=0.5)

    # ptx_hist(data, img_output)
    normalized_plot(data, img_output)
