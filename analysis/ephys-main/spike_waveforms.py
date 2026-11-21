import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import helper_functions as hf
from pathlib import Path

def waveform_plots(meta_loc, data_loc, output_loc, threshold):
    meta = pd.read_csv(os.path.join(meta_loc, "main_ephys_metadata.csv"))
    dv_dt_all = []
    voltage_all = []
    fig, ax = plt.subplots(1, 3, figsize=(3.5, 1))
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
                dv_dt_avg = np.mean(dv_dt, axis=0)
                dv_dt_avg = dv_dt_avg / 20
                dv_dt_all.append(dv_dt_avg)
                voltage_all.append(voltage_avg)

                # Plotting
                if gfp == "pos":
                    color = "green"
                    t = 1
                elif gfp == "neg":
                    color = "black"
                    t = 2
                ax[t].plot(
                    voltage_avg, dv_dt_avg, lw=0.4, color=color, alpha=0.5
                )
        # Average line for all spikes, and error (std)
        voltage_all_avg = np.mean(voltage_all, axis=0)
        dv_dt_all_avg = np.mean(dv_dt_all, axis=0)

        # Plotting average lines
        if gfp == "pos":
            color = "darkgreen"
            linestyle = "-"
        elif gfp == "neg":
            color = "black"
            linestyle = "--"
        ax[0].plot(
            voltage_all_avg, dv_dt_all_avg, color=color, lw=0.75, ls=linestyle
        )
    for a in ax.flatten():
        a.set_ylim(-1, 2)
        a.set_xlim(-75, 75)
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
        a.set_yticks([-1, 0, 1, 2])
        a.set_yticklabels(["", "", "", ""])
    ax[0].set_yticklabels([-1, 0, 1, 2])
    ax[0].set_ylabel("dV/dT (mV/μs)")
    ax[1].set_ylabel(" ")
    ax[2].set_ylabel(" ")
    # ax[1].set_xlabel("Vm (mV)")
    fig.tight_layout()
    plt.show()
    fig.savefig(os.path.join(output_loc, "dv_dt_plots.pdf"))
    return


if __name__ == "__main__":
    base_folder = Path.cwd()
    data_loc = os.path.join(base_folder, "data")
    meta_loc = os.path.join(base_folder, "metadata")
    output_loc = os.path.join(base_folder, "outputs")

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
