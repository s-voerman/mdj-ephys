# 3 or 4 plots combined
# max dV/dT, min dV/dT, max V, min V

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

import helper_functions as hf
from pathlib import Path

from scipy.stats import (
    ttest_ind,
    levene,
    kurtosistest,
    skewtest,
    mannwhitneyu,
    bartlett,
    describe,
    iqr,
    normaltest,
)


def waveform_plots(meta_loc, data_loc, output_loc, threshold):
    meta = pd.read_csv(os.path.join(meta_loc, "main_ephys_metadata.csv"))
    plot_dict = {}
    plot_df = pd.DataFrame()
    dv_pos = []
    dv_neg = []
    v_pos = []
    v_neg = []
    for gfp in meta.GFP.unique():
        t_meta = meta.loc[meta["GFP"] == gfp]
        for index, row in t_meta.iterrows():
            if index % 1 == 0:
                # Grab stuff from meta, and read data
                rheobase = (
                    str(row.file)
                    + " "
                    + str(row.cell)
                    + "-"
                    + str(row.rheo)
                    + ".ABF"
                )
                data, sr = hf.read_abf(os.path.join(data_loc, rheobase))

                # Grab spike index
                indices = hf.edge_detection(data, threshold)
                if len(indices) < 5:
                    indices = hf.edge_detection(
                        data, (threshold - 25)
                    )  # So -15 mV
                else:
                    indices = indices

                    # get 10 spikes, make averages
                voltage = []
                dv_dt = []
                for i, spike in enumerate(indices[0:10]):
                    margin = int((2 / 1000) * sr)  # space around spike
                    spike_data = data[spike - margin : spike + margin]
                    voltage.append(spike_data)
                    spike_data_dv = np.gradient(spike_data)
                    dv_dt.append(spike_data_dv)
                voltage_avg = np.mean(voltage, axis=0)
                print(voltage_avg)
                dv_dt_avg = np.mean(dv_dt, axis=0)
                plot_dict["voltage_max"] = np.max(voltage_avg)
                plot_dict["voltage_min"] = np.min(voltage_avg)
                plot_dict["dv_dt_max"] = (np.max(dv_dt_avg)) / 20
                plot_dict["dv_dt_min"] = (np.min(dv_dt_avg)) / 20
                # Dividing by 20 converts from mV/single sample to mV/us
                # 1 sample = 20 us as sr = 50 kHz
                plot_dict["gfp"] = gfp
                plot_df = pd.concat([plot_df, pd.DataFrame([plot_dict])])
                if gfp == "pos":
                    dv_pos.append(dv_dt_avg)
                    v_pos.append(voltage_avg)
                elif gfp == "neg":
                    dv_neg.append(dv_dt_avg)
                    v_neg.append(voltage_avg)
    plot_df.to_excel(os.path.join(output_loc, "maxspikes.xlsx"))
    # Stats & plotting
    vars = ["voltage_max", "dv_dt_max", "dv_dt_min"]
    test_output = {}
    desc_output = {}
    test_df = pd.DataFrame()
    desc_df = pd.DataFrame()
    fig, ax = plt.subplots(1, 3, figsize=(3.5, 1.75))
    for t, v in enumerate(vars):
        # Setup data
        pos = plot_df.loc[plot_df["gfp"] == "pos", v]
        neg = plot_df.loc[plot_df["gfp"] == "neg", v]

        # Test stats
        test_output["STAT"] = str(v)
        test_output["levene"] = levene(pos, neg)[1]
        test_output["ttest_ind"] = ttest_ind(pos, neg)[1]
        test_output["ttest_ind_ineqvar"] = ttest_ind(
            pos, neg, equal_var=False
        )[1]
        test_output["mwu"] = mannwhitneyu(pos, neg)[1]
        test_output["kurtosis_pos"] = kurtosistest(pos)[1]
        test_output["kurtosis_neg"] = kurtosistest(neg)[1]
        test_output["skew_pos"] = skewtest(pos)[1]
        test_output["skew_neg"] = skewtest(neg)[1]
        test_output["normal_pos"] = normaltest(pos)[1]
        test_output["normal_neg"] = normaltest(neg)[1]
        test_df = pd.concat(
            [test_df, pd.DataFrame([test_output])], ignore_index=True
        )

        # Descriptive stats
        for f, g in zip([pos, neg], ["pos", "neg"]):
            desc_output["STAT"] = str(v)
            desc_output["nobs_" + g] = int(len(f))
            desc_output["mean_" + g] = np.mean(f)
            desc_output["median_" + g] = np.median(f)
            desc_output["std_" + g] = np.std(f)
            desc_output["iqr" + g] = iqr(f)
            desc_output["minmax_" + g] = [float(np.min(f)), float(np.max(f))]
            desc_output["skew_" + g] = describe(f)[4]
            desc_output["kurt_" + g] = describe(f)[5]
        desc_df = pd.concat(
            [desc_df, pd.DataFrame([desc_output])], ignore_index=True
        )
        # Plotting
        colors = ["green", "white"]
        sns.boxplot(
            data=plot_df,
            y=v,
            hue="gfp",
            legend=False,
            ax=ax[t],
            gap=0.4,
            palette=colors,
            fliersize=0,
            linewidth=0.5,
        )
        sns.stripplot(
            data=plot_df,
            y=v,
            hue="gfp",
            legend=False,
            ax=ax[t],
            color="black",
            edgecolor="black",
            linewidth=0.5,
            dodge=True,
            jitter=0,
            size=2,
            marker="d",
            palette=["white", "white"],
        )
    for a in ax.flatten():
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
        a.set_xticklabels(["", ""])
        a.set_xlim(-0.4, 0.4)
        a.set_xticks([-0.2, 0.2])
        a.set_xticklabels(["GFP+", "GFP-"])
    ax[0].set_ylabel("Voltage (max - mV)")
    ax[1].set_ylabel("dV/dT (max - mV/μs)")
    ax[2].set_ylabel("dV/dT (min - mV/μs)")
    ax[0].set_ylim(20, 80)
    ax[1].set_ylim(0, 2)
    ax[2].set_ylim(-1, 0)
    fig.align_ylabels()
    fig.tight_layout()
    fig.savefig(os.path.join(output_loc, "max_dv_v_box.pdf"))
    test_df = test_df.set_index("STAT")
    test_df.to_excel(os.path.join(stat_loc, "waveform_stats.xlsx"))
    desc_df = desc_df.set_index("STAT")
    desc_df.to_excel(os.path.join(stat_loc, "waveform_descriptives.xlsx"))
    return


if __name__ == "__main__":
    base_folder = Path.cwd()
    data_loc = os.path.join(base_folder, "data")
    meta_loc = os.path.join(base_folder, "metadata")
    output_loc = os.path.join(base_folder, "output")
    stat_loc = os.path.join(base_folder, "output")

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

    waveform_plots(meta_loc, data_loc, output_loc, 10)
