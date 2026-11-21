import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

from scipy.stats import fisher_exact
from pathlib import Path

# This is the one for rebound, but you can probably get the idea how this would work for spontaneous firing and different GFP ID

def pie_charts(output_loc):
    data = pd.read_csv(os.path.join(output_loc, "spike_properties_filtered.csv"))
    data = data.loc[data.Current < 0]
    data = data.groupby(["Cell", "GFP"]).mean(numeric_only=True).reset_index()
    table = []
    for i, gfp in enumerate(data.GFP.unique()):
        d = data.loc[data["GFP"] == gfp]
        rebound_count = len(d.loc[d["Rebound spikes"] > 0])
        not_rebound_count = len(d.loc[d["Rebound spikes"] == 0])
        print(gfp)
        print(rebound_count, not_rebound_count)
        ratio = [rebound_count, not_rebound_count]
        labels = ["Y", "N"]
        colors = ["lightblue", "white"]
        fig, ax = plt.subplots(figsize=(1, 1))
        ax.pie(
            ratio,
            labels=None,
            colors=colors,
            wedgeprops={"edgecolor": "black", "linewidth": 0.5},
        )
        table.append(ratio)
        fig.tight_layout()
        fig.savefig(
            os.path.join(output_loc, "rebound_pie" + str(gfp) + ".pdf")
        )
    print(table)
    print(fisher_exact(table))
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
    pie_charts(output_loc)
