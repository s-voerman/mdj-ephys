import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import ttest_ind, levene
from pathlib import Path

def gfp_controls(data, output_loc):
    """
    Really basic plot for comparison of active/passive membrane properties of GFP+ / GFP- cells

    Parameters
    data: pandas DataFrame
    output_loc: path (str)
    """

    # Plotting vars
    vars = [
        "Rm",
        "Capacitance",
        "Vm",
    ]

    # Axes limits
    limits = [
        [0, 1000],
        [0, 150],
        [-70, -40],
    ]

    # Plotting
    fig, ax = plt.subplots(1, 3, figsize=(2.4, 1))
    colors = ["green", "white"]
    for i, v in enumerate(vars):
        sns.boxplot(
            y=v,
            data=data,
            hue="GFP",
            ax=ax[i],
            fliersize=0,
            legend=False,
            palette=colors,
            linecolor="black",
            linewidth=0.5,
        )
        sns.stripplot(
            y=v,
            data=data,
            hue="GFP",
            ax=ax[i],
            dodge=True,
            size=2,
            color="white",
            legend=False,
            jitter=0.05,
            marker="d",
            edgecolor="black",
            linewidth=0.5,
            palette=["white", "white"],
        )
        ax[i].set_ylim(limits[i][0], limits[i][1])

    # Aesthetics
    for a in ax.flatten():
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
        a.spines["bottom"].set_visible(False)
        a.get_xaxis().set_visible(False)

    # Basic stats
    gfp_data = data.loc[data.GFP == "pos"]
    n_gfp_data = data.loc[data.GFP == "neg"]
    for v in vars:
        t = ttest_ind(gfp_data[v], n_gfp_data[v])
        l = levene(gfp_data[v], n_gfp_data[v])
        print(str(v), ":")
        print(t, "t_test results")
        print(l, "levene results \n")

    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, "basic_controls.pdf"))
    plt.show()
    return


def gfp_spiking(spike_data, output_loc):
    # Append tag for maximum current injections
    plot_data = pd.DataFrame()
    for i, v in enumerate(spike_data["Cell"].unique()):
        temp_df = spike_data.loc[spike_data["Cell"] == v]
        temp_df.loc[:, "Max_current"] = (len(temp_df) * 5) - 5
        plot_data = pd.concat([plot_data, temp_df])
    colors = ["green", "black"]
    fig, ax = plt.subplots(2, 1, figsize=(2.8, 2.5))
    sns.lineplot(
        x="Current",
        y="Rebound spikes",
        hue="GFP",
        data=plot_data,
        ax=ax[0],
        legend=False,
        errorbar="sd",
        linewidth=0.5,
        palette=colors,
    )

    sns.lineplot(
        x="Current",
        y="Firing rate",
        hue="GFP",
        data=plot_data,
        ax=ax[1],
        legend=False,
        errorbar="sd",
        linewidth=0.5,
        palette=colors,
    )

    ax[0].set_xlim(0, -251)
    ax[1].set_xlim(0, 251)
    ax[0].set_xticks([-250, -200, -150, -100, -50, 0])
    ax[1].set_xticks([0, 50, 100, 150, 200, 250])
    ax[0].set_ylim(0, 8)
    ax[0].set_yticks([0, 4, 8])
    ax[1].set_ylim(0, 100)
    ax[0].axvline(x=-100, ls="--", lw=0.5, color="black")
    ax[0].axvline(x=-200, ls="--", lw=0.5, color="black")
    ax[0].axvline(x=-250, ls="--", lw=0.5, color="black")
    ax[1].axvline(x=100, ls="--", lw=0.5, color="black")
    ax[1].axvline(x=200, ls="--", lw=0.5, color="black")
    ax[1].axvline(x=250, ls="--", lw=0.5, color="black")
    ax[0].set_ylabel("")
    ax[1].set_ylabel("")

    for a in ax.flatten():
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
        a.set_xlabel("Injected current (pA)")
    ax[0].set_xlabel(" ")

    fig.tight_layout()
    fig.align_labels()
    fig.savefig(os.path.join(output_loc, "spiking_controls.pdf"), dpi=1200)
    plt.show()
    return


def gfp_spont(data, output_loc):
    """
    Really basic plot for comparison of active/passive membrane properties of GFP+ / GFP- cells

    Parameters
    data: pandas DataFrame
    output_loc: path (str)
    """

    # Plotting vars
    vars = [
        "Frequency",
        "CV1",
        "CV2",
    ]

    # Axes limits
    limits = [
        [0, 60],
        [0, 3],
        [0, 1],
    ]

    # Plotting
    fig, ax = plt.subplots(1, 3, figsize=(10, 5))
    colors = ["green", "white"]
    for i, v in enumerate(vars):
        sns.boxplot(
            y=v,
            data=data,
            hue="GFP",
            ax=ax[i],
            fliersize=0,
            legend=False,
            palette=colors,
            linecolor="black",
            linewidth=2,
        )
        sns.stripplot(
            y=v,
            data=data,
            hue="GFP",
            ax=ax[i],
            dodge=True,
            size=8,
            color="white",
            legend=False,
            jitter=0.05,
            marker="d",
            edgecolor="black",
            linewidth=2,
            palette=["white", "white"],
        )
        ax[i].set_ylim(limits[i][0], limits[i][1])

    # Aesthetics
    for a in ax.flatten():
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
        a.spines["bottom"].set_visible(False)
        a.get_xaxis().set_visible(False)
    ax[2].set_yticks([0, 0.5, 1])
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, "spont_controls.pdf"))
    plt.show()
    return


def gfp_spike_props(data, output_loc):
    """
    Really basic plot for comparison of active/passive membrane properties of GFP+ / GFP- cells

    Parameters
    data: pandas DataFrame
    output_loc: path (str)
    """

    # Plotting vars
    vars = [
        "fAHP",
        "Spike_threshold",
        "Half_width",
    ]

    # Axes limits
    limits = [
        [0, 40],
        [-60, -30],
        [0, 1000],
    ]

    # Plotting
    fig, ax = plt.subplots(1, 3, figsize=(2.4, 1))
    colors = ["green", "white"]
    for i, v in enumerate(vars):
        sns.boxplot(
            y=v,
            data=data,
            hue="GFP",
            ax=ax[i],
            fliersize=0,
            legend=False,
            palette=colors,
            linecolor="black",
            linewidth=0.5,
        )
        sns.stripplot(
            y=v,
            data=data,
            hue="GFP",
            ax=ax[i],
            dodge=True,
            size=2,
            color="black",
            legend=False,
            jitter=0.05,
            marker="d",
            linewidth=0.5,
            edgecolor="black",
            palette=["white", "white"],
        )
        ax[i].set_ylim(limits[i][0], limits[i][1])

    # Aesthetics
    for a in ax.flatten():
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
        a.spines["bottom"].set_visible(False)
        a.get_xaxis().set_visible(False)

    # Basic stats
    gfp_data = data.loc[data.GFP == "pos"]
    n_gfp_data = data.loc[data.GFP == "neg"]
    for v in vars:
        t = ttest_ind(gfp_data[v], n_gfp_data[v])
        l = levene(gfp_data[v], n_gfp_data[v])
        print(str(v), ":")
        print(t, "t_test results")
        print(l, "levene results \n")

    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, "gfp_spike_props.pdf"))
    plt.show()
    return


if __name__ == "__main__":
    # Folders
    base_folder = Path.cwd()
    meta_loc = os.path.join(base_folder, "metadat")
    data_loc = os.path.join(base_folder, "data")
    output_loc = os.path.join(base_folder, "output")

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

    # Reading data
    passive_data = pd.read_csv(os.path.join(data_loc, "passive_properties_filtered.csv"))
    spike_data = pd.read_excel(os.path.join(data_loc, "spike_properties_filtered.csv"))

