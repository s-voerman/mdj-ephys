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


def ttest_stats(passive_data, spike_data, output_loc):
    # Output
    output = {}
    stat_df = pd.DataFrame()
    # vars for passive data (these names are irrelevant really)
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
        output["STAT"] = str(v)
        output["levene"] = levene(pos, neg)[1]
        output["ttest_ind_eqvar"] = ttest_ind(pos, neg, equal_var=True)[1]
        output["ttest_ind_ineqvar"] = ttest_ind(pos, neg, equal_var=False)[1]
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
    stat_df.to_excel(os.path.join(output_loc, "passive_stats.xlsx"))
    return


if __name__ == "__main__":
    pd.options.display.float_format = "{:.5f}".format
    base_folder = Path.cwd()
    data_loc = os.path.join(base_folder, "data")
    meta_loc = os.path.join(base_folder, "metadata")
    output_loc = os.path.join(base_folder, "output")

    passive_data = pd.read_csv(os.path.join(data_loc, "passive_properties_filtered.csv"))
    spike_data = pd.read_csv(os.path.join(data_loc, "spike_properties_filtered.csv"))

    ttest_stats(passive_data, spike_data, output_loc)
