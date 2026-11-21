import os
import pandas as pd
import numpy as np
import helper_functions as hf
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from pathlib import Path

BASE_OUTPUT_DIR = Path(__file__).resolve().parent.parent.parent 

"""
Primary scripts that analyses raw ephys data
EphysAnalysis is a class that takes properties of ABF files and analyzes 1 cells complete set of experiments
Data obtained from this is split into two distinct outputs:
    passive_properties (which primarily contains non-current injection results)
    spike_properties (which primarily contains current injection results)
"""

class EphysAnalysis:
    """
    Class for reading, analyzing, and plotting ephys data from .ABF formats
    """

    def __init__(
        self,
        data_loc,
        file,
        cell,
        seal,
        sepsc,
        curr,
        rheo,
        spont,
        nucleus,
    ):
        self.data_loc = data_loc
        self.file = file
        self.cell = cell
        self.cellname = str(self.file) + " " + str(self.cell)
        self.seal = seal
        self.sepsc = sepsc
        self.curr = curr
        self.rheo = rheo
        self.spont = spont
        self.nucleus = nucleus
        self.threshold = -5

    def abf_reader(self, version, file):
        """
        Reads data from abf files in either version 1 or version 2

        If you are unsure what version you have:
        ABF1 is likely made with pClamp v < 10
        ABF2 is likely made with pClamp v >= 10
        If you convert to ABF from something like .itx or HEKA .dat using the ABF utility that comes with MiniAnalysis (Synaptosoft)
        you can safely assume you have version 1.

        Parameters
        version: int
            Refers to the version of the .ABF file to be read
            You can check the file version in clampfit by opening the file and pressing Alt + Enter
        file:
            name of file to be read

        I personally use a hacky method to go from HEKA .dat files -> export as .itx -> convert to .ABF V1.
        Due to this I lose any information on the stimulus when reading the files in neo, only preserving the analog signal.
        This is not the case if you use Molecular devices, which natively spits out .ABF files.
        """
        # assert isinstance(version, int)
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
        """
        Analyzes seal tests and returns passive membrane properties:
        Membrane, series, and input resistance
        Membrane capacitance
        Holding current (at Vcommand)

        Parameters
        data: numpy array
            Array containing data of at least 1 full 1 second sweep
        sr: int
            Sampling rate
        """

        # Detection bounds, change according to formatting of seal test
        sweep_size = 1.0
        bl = [0.65, 0.69]  # Window to measure baseline, in % of sweep
        ac = [0.69, 0.71]  # Window to measure access resistance (peak)
        mb = [
            0.78,
            0.795,
        ]  # Window to measure membrane resistance (steady state) - Made this significantly smaller was [0.75,0.79]
        delta_v = -0.01  # stimulus change in V

        output_df = pd.DataFrame()
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

            output["Sweep"] = i
            output["Ri"] = r_input
            output["Rs"] = r_series
            output["Rm"] = r_membrane
            output["Rs/Ri"] = r_ratio
            output["Holding"] = baseline
            output["Capacitance"] = cm
            output["R2"] = r_squared
            output["Cell"] = str(self.cellname)

            output_df = pd.concat(
                [output_df, pd.DataFrame([output])], ignore_index=True
            )
        output_df_avg = output_df.groupby(["Cell"]).mean().reset_index()
        return output_df, output_df_avg

    def curr_inj(self, data, sr, threshold):
        """
        Analysis of current injection traces

        Parameters:
        data: numpy array
        sr: int
        threshold: int
            threshold to detect spikes, see hf.edge_detection()
        """
        # Setting up output
        output = {}
        output_df = pd.DataFrame()

        # Bounds of spiking and rebound spiking
        nm_spikes = [0.25, 0.75]
        rb_spikes = [0.75, 1.0]

        sweep_size = 1.0  # Assuming 1 trace = 1 second
        sweep = int(sr * sweep_size)
        window = np.arange(0, len(data), sweep)
        for i, x in enumerate(window):
            # Setting up data
            d = data[x : x + sweep]
            vm = int(np.mean(d[int(0.6 * sr) : int(0.7 * sr)]))
            normal_data = d[int(nm_spikes[0] * sr) : int(nm_spikes[1] * sr)]
            rebound_data = d[int(rb_spikes[0] * sr) : int(rb_spikes[1] * sr)]

            # Get spike indices
            normal_indices = hf.edge_detection(normal_data, threshold)
            rebound_indices = hf.edge_detection(rebound_data, threshold)
            normal_indices = np.array(normal_indices)
            rebound_indices = np.array(rebound_indices)
            spike_count = len(normal_indices)
            rebound_count = len(rebound_indices)

            # Sweep counts, calculation of current injected and current/spike
            sweep_count = len(data) / sr
            current = sweep_count * -5 + (i * 10) + 5

            if current != 0:
                spikes_per_pa = spike_count / current
            else:
                spikes_per_pa = 0

            if current < 0:
                try:
                    ttf_rb = rebound_indices[0] * 0.02
                except IndexError:
                    ttf_rb = np.nan
            else:
                ttf_rb = np.nan
            try:
                ttf = (normal_indices[0]) * 0.02
            except IndexError:
                ttf = np.nan

            # Sag ratio
            baseline = [0.1, 0.24]
            steady_state = [0.6, 0.65]
            if current < -1:
                bl = np.mean(d[int(baseline[0] * sr) : int(baseline[1] * sr)])
                sag = np.min(normal_data)
                sts = np.mean(
                    d[int(steady_state[0] * sr) : int(steady_state[1] * sr)]
                )
                ratio = sag / sts
            elif current >= 0:
                ratio = np.nan

            # mAHP
            if current > 0:
                mAHP_bl = np.mean(
                    d[int(baseline[0] * sr) : int(baseline[1] * sr)]
                )
                mAHP_sag = np.min(
                    d[int(rb_spikes[0] * sr) : int(rb_spikes[1] * sr)]
                )  # rb spikes index should work pretty good for this
                mAHP = np.abs(mAHP_sag - mAHP_bl)
            else:
                mAHP = np.nan

            # Average spike amplitude (per sweep)
            baseline = [0.1, 0.24]
            base = np.mean(d[int(baseline[0] * sr) : int(baseline[1] * sr)])
            if len(normal_indices) > 0:
                amp_list = []
                spike_index = hf.edge_detection(d, threshold)
                for i, x in enumerate(spike_index):
                    spike_data = d[x - 50 : x + 50]
                    peak = np.max(spike_data)
                    spike_amplitude = peak - base
                    amp_list.append(spike_amplitude)
                avg_spike_amp = np.mean(amp_list)
            else:
                avg_spike_amp = np.nan

            # Outputs
            output["Cell"] = self.cellname
            output["Spikes"] = spike_count
            output["Firing rate"] = int(
                spike_count / (sweep_size - (nm_spikes[1] - nm_spikes[0]))
            )
            output["Rebound spikes"] = rebound_count
            output["Rebound_ttf"] = ttf_rb
            output["TTF"] = ttf
            output["Sag ratio"] = ratio
            output["mAHP"] = mAHP
            output["avg_spike_amplitude"] = avg_spike_amp
            output["Current"] = current
            output["Vm"] = vm
            output["Steps"] = 10
            output_df = pd.concat(
                [output_df, pd.DataFrame([output])], ignore_index=True
            )
            output["Spikes/pA"] = spikes_per_pa
        return output_df

    def spontaneous(self, data, sr, threshold):
        """
        Analysis of spontaneous current clamp traces.
        Injected current is always 0

        Parameters:
        data: numpy array
        sr: int
        threshold: int
            threshold to detect spikes, see hf.edge_detection()
        """

        # Setting up outputs
        output = {}
        output_df = pd.DataFrame()

        # Windows to analyze
        sweep_size = 1.0
        sweep = int(sr * sweep_size)
        window = np.arange(0, len(data), sweep)
        for i, x in enumerate(window):
            d = data[x : x + sweep]
            indices = hf.edge_detection(d, threshold)
            spike_count = len(indices)
            memb_voltage = np.mean(d)

            # Output
            output["Vm"] = memb_voltage
            output["Frequency"] = spike_count
            output["Cell"] = self.cellname
            output["Sweep"] = i
            output_df = pd.concat(
                [output_df, pd.DataFrame([output])], ignore_index=True
            )
        output_df_avg = output_df.groupby(["Cell"]).mean().reset_index()
        return output_df, output_df_avg

    def coefficients(self, data, sr, threshold):
        """
        Determines CV1/CV2 based on spontaneous firing traces
        CV1 = global
        CV2 = local

        Parameters
        data: numpy array
        sr: int
        threshold: int
        """

        # Setting up outputs
        output = {}
        output_df = pd.DataFrame()

        # Iterating every spike index for CV1 and CV2
        indices = hf.edge_detection(data, threshold)

        # CV1
        cv1_list = []
        output["Cell"] = self.cellname
        if len(indices) > 1:
            for i, v in enumerate(indices):
                loc_isi = indices[i] - indices[i - 1]
                loc_isi = loc_isi * (1000 / sr)
                cv1_list.append(loc_isi)
            cv1_list = cv1_list[1 : len(indices)]
            cv1 = np.std(cv1_list) / np.mean(cv1_list)
            output["CV1"] = cv1

            # CV2
            for i in range(0, len(indices) - 1):
                isi_1 = (indices[i] - indices[i - 1]) * (1000 / sr)
                isi_2 = (indices[i + 1] - indices[i]) * (1000 / sr)
                isi_list = [isi_1, isi_2]
                isi_list_mean = np.mean(isi_list)
                isi_list_std = np.std(isi_list)

                # Output
                output["ISI"] = isi_1
                output["CV2"] = isi_list_std / isi_list_mean
                output_df = pd.concat(
                    [output_df, pd.DataFrame([output])], ignore_index=True
                )
            output_df = output_df[1:]  # Drops first values which are wrong

        elif len(indices) < 2:
            print("Spontaneous frequency is ~0")
            output["Cell"] = self.cellname
            output["CV1"] = np.nan
            output["ISI"] = np.nan
            output["CV2"] = np.nan
            output_df = pd.concat(
                [output_df, pd.DataFrame([output])], ignore_index=True
            )
            output_df.reset_index(drop=True)
        output_df_avg = output_df.groupby(["Cell"]).mean().reset_index()
        return output_df, output_df_avg

    def rheobase(self, data, sr, threshold, plot=False):
        """
        Analyzes rheobase and basic spike properties:
        - Spike threshold
        - AHP
        - Max rise slope
        - Max decay slope
        Parameters:
        data: numpy array
        sr: int
        threshold: int
        plot: bool, optional
            Whether or not to draw plots for spiking threshold and AHP minimums.
            Just for making sure things are as they should be
            Default is False
        """
        # Setting up outputs
        output = {}
        output_df = pd.DataFrame()

        # Iterable, going with first 10 spikes for now
        indices = hf.edge_detection(data, threshold)
        if len(indices) < 5:
            indices = hf.edge_detection(data, (threshold - 10))  # So -15 mV

        for i, x in enumerate(indices[0:10]):
            margin = int((10 / 1000) * sr)  # space around spike
            d = data[x - margin : x + margin]

            # Spike threshold
            dv_dt = np.gradient(d)
            threshold_index = np.where(dv_dt > 1)[0][0]
            spike_threshold = d[threshold_index]

            # AHP
            ahp_peak_index = np.argmin(d[0:600])
            ahp_peak = np.min(d[0:600])
            fahp = spike_threshold - ahp_peak

            # Amplitude
            baseline = np.median(
                d[0:200]
            )  # Not sure if ok, or if I should take baseline pre current injection
            peak = np.max(d)
            spike_amplitude = peak - baseline

            # Half width
            ld = data[x - 100 : x + 100]  # try this for now
            half_peak = peak - (0.5 * spike_amplitude)
            rise_half = np.where(ld > half_peak)[0][0]
            decay_half = np.where(ld > half_peak)[0][-1]
            half_width = (
                (decay_half - rise_half) * (1000 / sr) * (1e3)
            )  # In uS

            output["Cell"] = self.cellname
            output["fAHP"] = fahp
            output["Spike_threshold"] = spike_threshold
            output["Spike_amplitude"] = spike_amplitude
            output["Half_width"] = half_width
            output["Count"] = i

            output_df = pd.concat(
                [output_df, pd.DataFrame([output])], ignore_index=True
            )

            # Optional plots
            if plot:
                fig, ax = plt.subplots(1, 2)
                ax[0].set_title("Half width")
                ax[0].plot(ld)
                ax[0].axvline(x=rise_half, color="black")
                ax[0].axvline(x=decay_half, color="yellow")

                ax[1].set_title("Threshold/AHP")
                ax[1].plot(d)
                ax[1].axvline(x=threshold_index, color="black")
                ax[1].axvline(x=ahp_peak_index, color="yellow")
                output_folder = None
                fig.savefig(
                    os.path.join(output_folder, self.cellname + ".png")
                )
                plt.close()
        output_df_avg = output_df.groupby(["Cell"]).mean().reset_index()
        return output_df, output_df_avg

    def run(self):
        print("Current cell: " + self.cellname)

        # Setting up outputs, terrible names
        spike_data = pd.DataFrame()
        other_data = pd.DataFrame()

        traces = [self.seal, self.curr, self.rheo, self.spont]
        for i, v in enumerate(traces):
            filename = str(self.cellname) + "-" + str(v) + ".ABF"
            try:
                data, sr = self.abf_reader(
                    1, os.path.join(self.data_loc, filename)
                )
            except FileNotFoundError:
                print(str(self.file) + "- File not available for this cell")
                continue

            # Creates seperate directory for output per cell
            analysis_output = BASE_OUTPUT_DIR
            if not os.path.exists(
                os.path.join(analysis_output, self.cellname)
            ):
                os.makedirs(os.path.join(analysis_output, self.cellname))
            output_dir = os.path.join(analysis_output, self.cellname)

            if v == self.seal:
                # Seals
                per_trace, average = self.seals(data, sr)
                other_data = pd.concat([other_data, average], axis=1)
                per_trace.to_excel(
                    os.path.join(
                        output_dir, (self.cellname + " seal_alltraces.xlsx")
                    )
                )
                average.to_excel(
                    os.path.join(
                        output_dir, (self.cellname + " seal_average.xlsx")
                    )
                )

            elif v == self.curr:
                # Current injections
                per_trace = self.curr_inj(data, sr, self.threshold)
                spike_data = pd.concat([spike_data, per_trace], axis=0)
                per_trace.to_excel(
                    os.path.join(
                        output_dir, (self.cellname + " curr_alltraces.xlsx")
                    )
                )

            elif v == self.rheo:
                # Rheobase
                per_trace, average = self.rheobase(
                    data, sr, self.threshold, plot=False
                )
                other_data = pd.concat([other_data, average], axis=1)
                per_trace.to_excel(
                    os.path.join(
                        output_dir, (self.cellname + " rheo_alltraces.xlsx")
                    )
                )
                average.to_excel(
                    os.path.join(
                        output_dir, (self.cellname + " rheo_average.xlsx")
                    )
                )

            elif v == self.spont:
                # Spontaneous, normal
                per_trace, average = self.spontaneous(data, sr, self.threshold)
                other_data = pd.concat([other_data, average], axis=1)
                per_trace.to_excel(
                    os.path.join(
                        output_dir, (self.cellname + " spont_alltraces.xlsx")
                    )
                )
                average.to_excel(
                    os.path.join(
                        output_dir, (self.cellname + " spont_average.xlsx")
                    )
                )

                # Spontaneous, coefficients
                per_trace, average = self.coefficients(
                    data, sr, self.threshold
                )
                other_data = pd.concat([other_data, average], axis=1)
                per_trace.to_excel(
                    os.path.join(
                        output_dir, (self.cellname + " coeff_alltraces.xlsx")
                    )
                )
                average.to_excel(
                    os.path.join(
                        output_dir, (self.cellname + " coeff_average.xlsx")
                    )
                )

            else:
                return None
        return spike_data, other_data


