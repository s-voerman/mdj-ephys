import os
import pandas as pd
import numpy as np
from scipy.stats import (
    ttest_ind,
    levene,
    kurtosistest,
    skewtest,
    mannwhitneyu,
    bartlett,
    describe,
)
from pathlib import Path


def descriptives(data_loc, output_loc):
    # Set up data
    all_data = pd.DataFrame()
    files = os.listdir(os.path.join(data_loc))
    for f in files:
        data = pd.read_excel(os.path.join(data_loc, f))
        data["frequency"] = len(data) / 120
        all_data = pd.concat([all_data, data])
    all_data = all_data.groupby(["cell"]).mean(numeric_only=True).reset_index()
    all_data["tau"] = all_data["tau"] * 1000
    all_data["rise"] = all_data["rise"] * 1000

    # Actual stats
    output = {}
    stat_df = pd.DataFrame()
    vars = ["frequency", "tau", "rise", "amplitude"]
    for i, v in enumerate(vars):
        f = all_data[v]
        output["STAT"] = str(v)
        output["nobs"] = int(len(f))
        output["mean"] = np.mean(f)
        output["median"] = np.median(f)
        output["std"] = np.std(f)
        output["minmax"] = [float(np.min(f)), float(np.max(f))]
        output["skew"] = describe(f)[4]
        output["kurt"] = describe(f)[5]
        stat_df = pd.concat(
            [stat_df, pd.DataFrame([output])], ignore_index=True
        )
    stat_df = stat_df.set_index("STAT")
    stat_df.to_excel(os.path.join(output_loc, "ipsc_descriptives.xlsx"))
    return


if __name__ == "__main__":
    pd.options.display.float_format = "{:.5f}".format
    base_folder = Path.cwd()
    data_loc = os.path.join(base_folder, "data")
    meta_loc = os.path.join(base_folder, "metadata")
    output_loc = os.path.join(base_folder, "output")
    descriptives(data_loc, output_loc)
