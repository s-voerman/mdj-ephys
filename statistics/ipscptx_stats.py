import os

import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

from itertools import combinations
from pathlib import Path

from scipy.stats import (
    ttest_ind,
    levene,
    kurtosistest,
    skewtest,
    mannwhitneyu,
    bartlett,
    describe,
    ttest_rel,
)

# Issue with the related ttest as some cells dont make it past 14 minutes
# So then I can't compare all cells together
# Which is really annoying
# So maybe let's just do like an RM ANOVA or something..


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
        cell_data["tau"] = cell_data["tau"] * 1000
        cell_data["rise"] = cell_data["rise"] * 1000
        # data normalization
        # probably normalize per trace, where baseline = 1
        vars = ["amplitude", "frequency", "rise", "tau"]
        for i, v in enumerate(vars):
            cell_data.loc[:, "norm_" + v] = cell_data.loc[:, v] / np.mean(
                cell_data.loc[(cell_data.trace == 0), v]
            )
        combined_data = pd.concat([combined_data, cell_data])
    return combined_data


def descriptives(data, output_loc):
    data = (
        data.groupby(["cell", "trace"]).mean(numeric_only=True).reset_index()
    )

    output = {}
    stat_df = pd.DataFrame()
    vars = [
        "norm_frequency",
        "norm_amplitude",
        "norm_rise",
        "norm_tau",
        "frequency",
        "amplitude",
        "rise",
        "tau",
    ]
    for i, v in enumerate(vars):
        for t in data.trace.unique():
            d = data.loc[data["trace"] == t]
            f = d[v]
            output["STAT"] = str(v) + "_" + str(t)
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
    stat_df.to_excel(os.path.join(output_loc, "ipscptx_descriptives.xlsx"))
    return


def ipsc_tests(data, output_loc):
    data = (
        data.groupby(["cell", "trace"]).mean(numeric_only=True).reset_index()
    )
    vars = ["norm_frequency", "norm_amplitude", "norm_rise", "norm_tau"]
    traces = data.trace.unique()
    cc = list(combinations(traces, 2))
    cc = cc[10:16]
    stat_df = pd.DataFrame()
    for f in cc:
        output = {}
        pre = data.loc[data["trace"] == f[0]]
        post = data.loc[data["trace"] == f[1]]
        print(post["norm_amplitude"])
        for i, v in enumerate(vars):
            f_pre = pre[v]
            f_post = post[v]
            output["STAT"] = str(v) + "_" + str(f[0]) + "-" + str(f[1])
            output["levene"] = levene(f_pre, f_post)[1]
            output["bartlett"] = bartlett(f_pre, f_post)[1]
            output["ttest_rel"] = ttest_rel(f_pre, f_post)[1]
            output["mwu"] = mannwhitneyu(f_pre, f_post)[1]
            output["kurtosis_pre"] = kurtosistest(f_pre)[1]
            output["kurtosis_post"] = kurtosistest(f_post)[1]
            output["skew_pre"] = skewtest(f_pre)[1]
            output["skew_post"] = skewtest(f_post)[1]
            stat_df = pd.concat(
                [stat_df, pd.DataFrame([output])], ignore_index=True
            )
    stat_df = stat_df.set_index("STAT")
    stat_df.to_excel(os.path.join(output_loc, "ipscptx_stats.xlsx"))
    return


if __name__ == "__main__":
    pd.options.display.float_format = "{:.5f}".format
    base_folder = Path.cwd()
    data_loc = os.path.join(base_folder, "data")
    meta_loc = os.path.join(base_folder, "metadata")
    output_loc = os.path.join(base_folder, "output")

    # read data
    meta = pd.read_csv(os.path.join(meta_loc, "ipsc_ptx.csv"))
    data = read_data(meta, data_loc, output_loc)
    ipsc_tests(data, output_loc)
