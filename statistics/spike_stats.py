import os
import pandas as pd
import numpy as np
from scipy.stats import (
    ttest_ind,
    levene,
    kurtosistest,
    skewtest,
    mannwhitneyu,
    normaltest,
)
from pathlib import Path

def reg_spike_stats(spike_data, output_loc):
    # Outputs
    stat_df = pd.DataFrame()
    output = {}

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
            output["STAT"] = str(v) + "_" + str(c)
            output["levene"] = levene(pos, neg)[1]
            output["ttest_ind_eqvar"] = ttest_ind(pos, neg, equal_var=True)[1]
            output["ttest_ind_ineqvar"] = ttest_ind(pos, neg, equal_var=False)[
                1
            ]
            output["mwu"] = mannwhitneyu(pos, neg)[1]
            output["kurtosis_pos"] = kurtosistest(pos)[1]
            output["kurtosis_neg"] = kurtosistest(neg)[1]
            output["skew_pos"] = skewtest(pos)[1]
            output["skew_neg"] = skewtest(neg)[1]
            output["normal_pos"] = normaltest(pos)[1]
            output["normal_neg"] = normaltest(neg)[1]
            stat_df = pd.concat(
                [stat_df, pd.DataFrame([output])], ignore_index=True
            )
    stat_df = stat_df.set_index("STAT")
    stat_df.to_excel(os.path.join(output_loc, "spike_stats.xlsx"))
    return


def reb_spike_stats(spike_data, output_loc):
    # Outputs
    stat_df = pd.DataFrame()
    output = {}

    # Currents + vars
    currents = [-100, -250, -500]
    vars = ["Rebound spikes", "Rebound_ttf", "Sag ratio"]
    for c in currents:
        d = spike_data.loc[spike_data["Current"] == c]
        for i, v in enumerate(vars):
            pos = d.loc[d["GFP"] == "pos", v]
            neg = d.loc[d["GFP"] == "neg", v]
            pos = pos.dropna()
            neg = neg.dropna()
            output["STAT"] = str(v) + "_" + str(c)
            output["levene"] = levene(pos, neg)[1]
            output["ttest_ind_eqvar"] = ttest_ind(pos, neg, equal_var=True)[1]
            output["ttest_ind_ineqvar"] = ttest_ind(pos, neg, equal_var=False)[
                1
            ]
            output["mwu"] = mannwhitneyu(pos, neg)[1]
            output["kurtosis_pos"] = kurtosistest(pos)[1]
            output["kurtosis_neg"] = kurtosistest(neg)[1]
            output["skew_pos"] = skewtest(pos)[1]
            output["skew_neg"] = skewtest(neg)[1]
            output["normal_pos"] = normaltest(pos)[1]
            output["normal_neg"] = normaltest(neg)[1]
            stat_df = pd.concat(
                [stat_df, pd.DataFrame([output])], ignore_index=True
            )
    stat_df = stat_df.set_index("STAT")
    stat_df.to_excel(os.path.join(output_loc, "rebound_stats.xlsx"))
    return


if __name__ == "__main__":
    pd.options.display.float_format = "{:.5f}".format

    base_folder = Path.cwd()
    data_loc = os.path.join(base_folder, "data")
    meta_loc = os.path.join(base_folder, "metadata")
    output_loc = os.path.join(base_folder, "output")

    passive_data = pd.read_csv(os.path.join(data_loc, "passive_properties_filtered.csv"))
    spike_data = pd.read_csv(os.path.join(data_loc, "spike_properties_filtered.csv"))
    reg_spike_stats(spike_data, output_loc)
    reb_spike_stats(spike_data, output_loc)
