import os
import pandas as pd
import numpy as np
from scipy.stats import describe, iqr
from pathlib import Path

def passive_descriptives(passive_data, output_loc):
    # Output
    output = {}
    stat_df = pd.DataFrame()
    # vars for passive data
    vars = [
        "Ri",
        "Rs",
        "Rm",
        "Capacitance",
        "fAHP",
        "Spike_threshold",
        "Half_width",
        "Vm",
        "Frequency",
        "CV1",
        "CV2",
    ]
    for i, v in enumerate(vars):
        pos = passive_data.loc[passive_data["GFP"] == "pos", v]
        neg = passive_data.loc[passive_data["GFP"] == "neg", v]
        pos = pos.dropna()
        neg = neg.dropna()
        for f, g in zip([pos, neg], ["pos", "neg"]):
            print(f)
            output["STAT"] = str(v)
            output["nobs_" + g] = int(len(f))
            output["mean_" + g] = np.mean(f)
            output["median_" + g] = np.median(f)
            output["std_" + g] = np.std(f)
            output["iqr" + g] = iqr(f)
            output["minmax_" + g] = [float(np.min(f)), float(np.max(f))]
            output["skew_" + g] = describe(f)[4]
            output["kurt_" + g] = describe(f)[5]
        stat_df = pd.concat(
            [stat_df, pd.DataFrame([output])], ignore_index=True
        )
    stat_df = stat_df.set_index("STAT")
    stat_df.to_excel(os.path.join(output_loc, "passive_descriptives.xlsx"))
    return


def spike_descriptives(spike_data, output_loc):
    # Output
    output = {}
    stat_df = pd.DataFrame()

    # Currents to compare + vars to compare
    currents = [100, 250, 500]
    vars = ["Firing rate", "TTF", "mAHP"]
    for c in currents:
        d = spike_data.loc[spike_data["Current"] == c]
        for i, v in enumerate(vars):
            pos = d.loc[d["GFP"] == "pos", v]
            neg = d.loc[d["GFP"] == "neg", v]
            pos = pos.dropna()
            neg = neg.dropna()
            for f, g in zip([pos, neg], ["pos", "neg"]):
                output["STAT"] = str(v) + "_" + str(c)
                output["nobs_" + g] = int(len(f))
                output["mean_" + g] = np.mean(f)
                output["median_" + g] = np.median(f)
                output["std_" + g] = np.std(f)
                output["iqr" + g] = iqr(f)
                output["minmax_" + g] = [float(np.min(f)), float(np.max(f))]
                output["skew_" + g] = describe(f)[4]
                output["kurt_" + g] = describe(f)[5]
            stat_df = pd.concat(
                [stat_df, pd.DataFrame([output])], ignore_index=True
            )
    stat_df = stat_df.set_index("STAT")
    stat_df.to_excel(os.path.join(output_loc, "spike_descriptives.xlsx"))
    return


def rebound_descriptive(spike_data, output_loc):
    # Output
    output = {}
    stat_df = pd.DataFrame()

    # Currents to compare + vars to compare
    currents = [-100, -250, -500]
    vars = ["Rebound spikes", "Rebound_ttf", "Sag ratio"]
    for c in currents:
        d = spike_data.loc[spike_data["Current"] == c]
        for i, v in enumerate(vars):
            pos = d.loc[d["GFP"] == "pos", v]
            neg = d.loc[d["GFP"] == "neg", v]
            pos = pos.dropna()
            neg = neg.dropna()
            for f, g in zip([pos, neg], ["pos", "neg"]):
                output["STAT"] = str(v) + "_" + str(c)
                output["nobs_" + g] = int(len(f))
                output["mean_" + g] = np.mean(f)
                output["median_" + g] = np.median(f)
                output["std_" + g] = np.std(f)
                output["iqr" + g] = iqr(f)
                output["minmax_" + g] = [float(np.min(f)), float(np.max(f))]
                output["skew_" + g] = describe(f)[4]
                output["kurt_" + g] = describe(f)[5]
            stat_df = pd.concat(
                [stat_df, pd.DataFrame([output])], ignore_index=True
            )
    stat_df = stat_df.set_index("STAT")
    stat_df.to_excel(os.path.join(output_loc, "rebound_descriptives.xlsx"))
    return


if __name__ == "__main__":
    pd.options.display.float_format = "{:.5f}".format

    base_folder = Path.cwd()
    data_loc = os.path.join(base_folder, "data")
    meta_loc = os.path.join(base_folder, "metadata")
    output_loc = os.path.join(base_folder, "output")

    passive_data = pd.read_csv(os.path.join(data_loc, "passive_properties_filtered.csv"))
    spike_data = pd.read_csv(os.path.join(data_loc, "spike_properties_filtered.csv"))

    passive_descriptives(passive_data, output_loc)
    spike_descriptives(spike_data, output_loc)
    rebound_descriptive(spike_data, output_loc)
