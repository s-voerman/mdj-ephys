import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress, pearsonr, spearmanr
from pathlib import Path


def spike_regression(spike_data, output_loc):
    # data = spike_data.groupby(["Current", "GFP"]).mean().reset_index()
    spike_data = spike_data.loc[spike_data["Current"] > -1]
    pos = spike_data.loc[spike_data["GFP"] == "pos"]
    neg = spike_data.loc[spike_data["GFP"] == "neg"]
    vars = ["Firing rate"]
    fig, ax = plt.subplots()
    for f, g in zip([pos, neg], ["pos", "neg"]):
        for i, v in enumerate(vars):

            y = f[v]
            x = f["Current"]
            res = linregress(x, y)
            print(f"R-squared: {res.rvalue**2:.6f}")
            print(f"Slope: {res.slope**2:.6f}")
            if g == "pos":
                color = "green"
            elif g == "neg":
                color = "black"
            ax.plot(x, y, "o", label="original data", color=color, lw=0.5)
            ax.plot(x, res.intercept + res.slope * x, "r", label="fitted line")
            plt.legend()

    return


if __name__ == "__main__":
    pd.options.display.float_format = "{:.5f}".format

    base_folder = Path.cwd()
    data_loc = os.path.join(base_folder, "data")
    meta_loc = os.path.join(base_folder, "metadata")
    output_loc = os.path.join(base_folder, "output")
   
    spike_data = pd.read_csv(os.path.join(data_loc, "spike_properties_filtered.csv"))

    spike_data = pd.read_excel(os.path.join(data_loc, "spike_data_f.xlsx"))
    spike_regression(spike_data, output_loc)
