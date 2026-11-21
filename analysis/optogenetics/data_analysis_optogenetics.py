import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import helper_functions as hf

from scipy.optimize import curve_fit
from pathlib import Path


class OptoAnalysis:
    def __init__(
        self, data_loc, cell, seal, sep, adjust, sep_cc, adjust_cc, gfp
    ):
        self.data_loc = data_loc
        self.cell = cell
        self.seal = seal
        self.sep = sep
        self.adjust = adjust
        self.sep_cc = sep_cc
        self.adjust_cc = adjust_cc
        self.gfp = gfp

    def read_abf(self, version, file):
        if version == 1:
            data, sr = hf.read_abf(os.path.join(self.data_loc, file))
        elif version == 2:
            data, sr = hf.read_abf_v2(os.path.join(self.data_loc, file))
        else:
            print("Enter a valid version (1 or 2)")
            data = None
            sr = None
        return data, sr

    def seals(self, data, sr):
        # Detection bounds, change according to formatting of seal test
        sweep_size = 1.0
        bl = [0.65, 0.69]  # Window to measure baseline, in % of sweep
        ac = [0.69, 0.71]  # Window to measure access resistance (peak)
        mb = [
            0.75,
            0.79,
        ]  # Window to measure membrane resistance (steady state)
        delta_v = -0.01  # stimulus change in V

        # Outputs
        output_df = pd.DataFrame()

        # Iterate over each sweep, do math
        sweep = sr * sweep_size
        window = np.arange(0, len(data), sweep)
        for i, x in enumerate(window):
            output = {}
            d = data[int(x) : int(x) + int(sweep)]

            baseline = np.mean(d[int(bl[0] * sweep) : int(bl[1] * sweep)])
            peak = np.min(d[int(ac[0] * sweep) : int(ac[1] * sweep)])
            sts = np.mean(d[int(mb[0] * sweep) : int(mb[1] * sweep)])

            # Math for resistances
            r_input = delta_v / ((sts - baseline) * 1e-12)  # in Ohms
            r_series = delta_v / ((peak - baseline) * 1e-12)  # in Ohms
            r_input = r_input / 1e6  # in MOhms
            r_series = r_series / 1e6  # in MOhms
            r_membrane = r_input - r_series
            r_ratio = r_series / r_input

            # Fitting for capacitance calc (double exp)
            fit_data_start = np.argmax(d)
            y = d[fit_data_start : fit_data_start + 500]
            x = np.arange(0, len(y))

            p0 = [np.max(d), 1.5, np.max(d), 0.5]
            popt, pcov = curve_fit(hf.double_exp, x, y, p0)
            a, tau, c, tau2 = popt
            if tau > tau2:
                tau_ms = tau * (1000 / sr)
            elif tau2 > tau:
                tau_ms = tau2 * (1000 / sr)
            else:
                tau_ms = tau * (1000 / sr)
                print("tau1 and tau2 are equal")
            r_parallel = 1 / ((1 / r_series) + (1 / r_membrane))
            cm = ((tau_ms * 1e-3) / (r_parallel * 1e6)) * 1e12

            # Fit quality
            squared_diffs = np.square(y - hf.double_exp(x, a, tau, c, tau2))
            squared_diffs_from_mean = np.square(y - np.mean(y))
            r_squared = 1 - np.sum(squared_diffs) / np.sum(
                squared_diffs_from_mean
            )

            # Outputting to dict and df
            output["Sweep"] = i
            output["Ri"] = r_input
            output["Rs"] = r_series
            output["Rm"] = r_membrane
            output["Rs/Ri"] = r_ratio
            output["Holding"] = baseline
            output["Capacitance"] = cm
            output["R2"] = r_squared
            output["Cell"] = str(self.cell)

            output_df = pd.concat(
                [output_df, pd.DataFrame([output])], ignore_index=True
            )
        output_df_avg = output_df.groupby(["Cell"]).mean().reset_index()
        return output_df, output_df_avg

    def opto_traces(self, data, sr, sweep_size, bounds, clamp, tag, plot=False):
        # Output
        output_df = pd.DataFrame()
        sweep = sr * sweep_size
        window = np.arange(0, len(data), sweep)
        for i, x in enumerate(window):
            output = {}
            d = data[int(x) : int(x) + int(sweep)]

            # Chronos
            chronos_bl = np.mean(
                d[int((bounds[0] * sr) - 5000) : int(bounds[0] * sr)]
            )
            if clamp == "cc":
                chronos_peak = np.max(
                    d[int(bounds[0] * sr) : int(bounds[0] * sr) + 5000]
                )
                chronos_peak_index = np.argmax(
                    d[int(bounds[0] * sr) : (int((bounds[0] * sr) + 5000))]
                )

            elif clamp == "vc":
                chronos_peak = np.min(
                    d[int(bounds[0] * sr) : int(bounds[0] * sr) + 5000]
                )
                chronos_peak_index = np.argmin(
                    d[int(bounds[0] * sr) : int((bounds[0] * sr) + 5000)]
                ) + (bounds[0] * sr)
            chronos_amp = chronos_peak - chronos_bl

            # Chrimson
            chrimson_bl = np.mean(
                d[int((bounds[1] * sr) - 5000) : int(bounds[1] * sr)]
            )
            if clamp == "cc":
                chrimson_peak = np.max(
                    d[int(bounds[1] * sr) : (int((bounds[1] * sr) + 5000))]
                )
                chrimson_peak_index = np.argmax(
                    d[int(bounds[1] * sr) : (int((bounds[1] * sr) + 5000))]
                ) + (bounds[1] * sr)
            elif clamp == "vc":
                chrimson_peak = np.min(
                    d[int(bounds[1] * sr) : (int((bounds[1] * sr) + 5000))]
                )
                chrimson_peak_index = np.argmin(
                    d[int(bounds[1] * sr) : (int((bounds[1] * sr) + 5000))]
                ) + (bounds[1] * sr)
            chrimson_amp = chrimson_peak - chrimson_bl

            # Plotting for QA
            if plot:
                if not os.path.exists(os.path.join(qa_loc, str(self.cell))):
                    os.makedirs(os.path.join(qa_loc, str(self.cell)))
                fig, ax = plt.subplots()
                ax.plot(d)
                # ax.plot(
                #     d[int(bounds[0] * sr) : int((bounds[0] * sr) + 5000)], color="green"
                # )
                # ax.plot(
                #     d[int(bounds[1] * sr) : int((bounds[1] * sr) + 5000)],
                #     color="purple",
                # )
                ax.axvline(x=chronos_peak_index, color="blue", lw=0.5, ls="--")
                ax.axvline(x=chrimson_peak_index, color="red", lw=0.5, ls="--")
                if clamp == "cc":
                    ax.set_ylim(-80, 0)
                elif clamp == "vc":
                    ax.set_ylim(0, -500)
                fig.tight_layout()
                fig.savefig(
                    os.path.join(
                        qa_loc + "/" + self.cell,
                        str(tag) + "_" + self.cell + str(i) + ".png",
                    )
                )
                # plt.show()
                plt.close()

            output["cell"] = self.cell
            output["gfp"] = self.gfp
            output["clamp"] = clamp
            output["chronos_amp"] = chronos_amp
            output["chrimson_amp"] = chrimson_amp
            output_df = pd.concat(
                [output_df, pd.DataFrame([output])], ignore_index=True
            )
        output_df = (
            output_df.groupby(["cell", "clamp", "gfp"]).mean().reset_index()
        )
        return output_df

    def run(self):
        output_df = pd.DataFrame()

        traces = [self.sep, self.sep_cc, self.adjust, self.adjust_cc]
        clamps = ["vc", "cc", "vc", "cc"]
        bounds = [[0.1, 0.6], [0.1, 1.1], [2.3, 0.3], [2.3, 0.3]]
        tag = ["sep_vc", "sep_cc", "adjust_vc", "adjust_cc"]
        seal_file = str(self.cell) + "-" + str(self.seal) + ".ABF"
        seal_data, sr = hf.read_abf(os.path.join(data_loc, seal_file))
        old, seal_output = self.seals(seal_data, sr)

        for trace, clamp, bound, tag in zip(traces, clamps, bounds, tag):
            print(str(self.cell + "-" + str(trace)))
            if trace == 0:
                continue
            filename = str(self.cell) + "-" + str(trace) + ".ABF"
            data, sr = hf.read_abf(os.path.join(data_loc, filename))
            sweep_size = (len(data) / 20) / sr
            output = self.opto_traces(
                data, sr, sweep_size, bound, clamp, tag, plot=False
            )
            output.loc[:, "tag"] = tag
            output_df = pd.concat([output_df, output])
        return output_df, seal_output


def read_meta(meta_loc, data_loc, output_loc):
    output_df = pd.DataFrame()
    seal_df = pd.DataFrame()
    meta = pd.read_excel(os.path.join(meta_loc, "MDJ_opto_V2.xlsx"))
    meta = meta.loc[meta.gfp == "pos"]  # Only GFP+ cells
    for index, row in meta.iterrows():
        cell = str(row.file) + " " + str(row.cell)
        seal = row.seal
        sep = row.sep
        adjust = row.adjust  # Pretty non-descriptive name, maybe change
        sep_cc = row.sep_cc
        adjust_cc = row.adjust_cc
        gfp = row.gfp
        p = OptoAnalysis(
            data_loc, cell, seal, sep, adjust, sep_cc, adjust_cc, gfp
        )
        output, seal_output = p.run()
        output_df = pd.concat([output_df, output])
        seal_df = pd.concat([seal_df, seal_output])
    output_df.to_excel(os.path.join(output_loc, "opto_output_gfppos.xlsx"))
    seal_df.to_excel(os.path.join(output_loc, "opto_output_seal.xlsx"))
    return output_df, seal_df


if __name__ == "__main__":
    base_folder = Path(__file__).resolve().parent.parent 
    data_loc = os.path.join(base_folder, "data")
    meta_loc = os.path.join(base_folder, "metadata")
    output_loc = os.path.join(base_folder, "output/optogenetics")
    qa_loc = os.path.join(output_loc, "output/optogenetics/qa")

    read_meta(meta_loc, data_loc, output_loc)
