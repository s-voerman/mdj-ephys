import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import helper_functions as hf

from scipy.signal import (
    deconvolve,
)

from process_events import analyze_event


class ProcessABF:
    """
    Deconvolves epsc trace based on template and determines threshold
    Passes event index to process_events
    Returns event parameters
    """

    def __init__(self, file_loc, output_loc, plot=False):
        # parameters
        self.file_loc = file_loc
        self.output_loc = output_loc
        self.plot = plot
        self.sr = 1
        self.original_sr = 1
        self.multiplier = 1  # was 3.
        self.threshold = None
        self.file_length = 0

        # data stores
        self.original_data = None  # contains untouched data
        self.filtered_data = None  # contains filtered data
        self.detection_data = None  # contains detrended, resampled data
        self.deconvoluted_data = None  # contains results of deconvolution

        # template vars
        self.rise = 0.0025
        self.decay = 0.009
        self.amplitude = -45

        # amp -35 -> -45
        # I think for MDJ I should probably go for shorter rise and decay as Ri tends to be a bit higher than PCs
        # Adapt amplitude for ipscs, tend to be a lot bigger.

    def get_detection_vars(self):
        """
        Creates some variables for event detection
        """
        # event detection
        # ideally 20 ms
        # Maybe mess around with this a bit to see if it improves anything
        self.min_distance = 0.010 * self.sr

    def create_template(
        self, rise=None, decay=None, amplitude=None, template_plot=False
    ):
        """to properly deconvolve the data, we need a reliable
        template. note that the length of the template should be
        relatively long, as short templates induce oscillations

        Parameters
        rise: int
            rise time of event
        decay: int
            decay time of event
        amplitude: int
            amplitude of event
        plot: bool (optional)
            whether to plot the template
        """
        if rise is None:
            rise = self.rise
        if decay is None:
            decay = self.decay
        if amplitude is None:
            amplitude = self.amplitude
        rise = rise * self.sr
        decay = decay * self.sr
        amp_normalization = (decay / rise) ** (rise / (rise - decay))
        template_size = round(0.05 * self.sr)
        g_x = np.arange(0, template_size)
        normalized_amplitude = amplitude / amp_normalization
        g_y = [
            (
                normalized_amplitude
                * (-np.exp(-t / rise) + np.exp(-t / decay))
                if t > 0
                else 0
            )
            for t in g_x
        ]
        if template_plot:
            fig, ax = plt.subplots(figsize=(2, 2))
            ax.plot(g_x, g_y)
            fig.tight_layout()
            plt.show()
        return g_y, g_x

    def deconvolution(
        self, data, template, wiener_flag=True, wiener_order=255
    ):
        """
        Deconvolves data based on previously built template
        """
        if wiener_flag:
            data = hf.wiener(data, wiener_order)
        minute_in_samples = 20 * self.sr
        data_decon = np.array([])
        minute_windows = np.arange(0, len(data), minute_in_samples)
        for i, x in enumerate(minute_windows):
            print("Window", i)
            start = x
            if i == len(minute_windows) - 1:
                end = len(data)
            else:
                end = x + minute_in_samples
            dd = deconvolve(data[start:end], np.asarray(template[2:]))[0]
            data_decon = np.concatenate([data_decon, dd])
        dd = None
        data_decon[np.where(data_decon > 1.5)] = 0
        data_decon[np.where(data_decon < -1.5)] = 0
        return data_decon, data

    def get_threshold(self, data_decon):
        """
        Determines threshold above which the deconvolved signal is considered an event
        """
        thresholds = []
        for x in np.arange(0, self.file_length, 3):
            x = x * self.sr
            step = 3 * self.sr
            m = np.median(data_decon[x : x + step])
            s = np.std(data_decon[x : x + step])
            threshold = m + self.multiplier * s
            thresholds.append(np.abs(threshold))
        threshold = np.min(thresholds)
        return threshold

    def propagate_threshold(self, data_decon):
        """
        Applies threshold to deconvolved data
        """
        self.get_detection_vars()
        threshold = self.get_threshold(data_decon)
        above_thresh = []
        for i, e in enumerate(data_decon[1:-1]):
            i = i + 1
            if ((e - np.mean(data_decon[i - 50 : i])) > threshold) and (
                (e - np.mean(data_decon[i + 1 : i + 51])) > threshold
            ):
                above_thresh.append(i)

        indices = [
            i + 1
            for (x, y, i) in zip(
                above_thresh, above_thresh[1:], range(len(above_thresh))
            )
            if self.min_distance < abs(x - y)
        ]
        events = [
            above_thresh[s:e]
            for s, e in zip([0] + indices, indices + [len(above_thresh)])
        ]

        output_events = np.array([])
        if len(events) == 1:
            if len(events[0]) == 0:
                return None, None
        for e in events:
            e = np.array(e).astype(int)
            m = np.max(data_decon[e])
            event_peak = e[np.where(data_decon[e] == m)[0]]
            output_events = np.append(output_events, event_peak[0])

        delta = len(self.detection_data) - len(data_decon)
        single_delta = delta / 6
        output_events = np.array(
            [
                x + (single_delta * int(x / (20 * self.sr)))
                for x in output_events
            ]
        )
        return events, output_events, threshold

    def quick_plot(
        self, threshold, output_loc, peaks=None, plot_threshold=False
    ):
        """
        Really basic plot with event index overlayed with data and deconvolution trace
        Not really that useful
        """
        fig, ax = plt.subplots(3, 1, sharex=True)
        for i, y in enumerate(
            [self.filtered_data, self.detection_data, self.deconvoluted_data]
        ):
            ax[i].plot(
                np.linspace(0, self.file_length, len(y)), y, linewidth=0.5
            )
        ax[0].set_xlim(0, 1)
        ax[0].set_title(self.file_loc)
        if plot_threshold:
            ax[2].axhline(threshold, color="red", linestyle=":")
        if peaks.any():
            peaks = peaks / self.sr
            peaks = peaks * self.original_sr
            peaks = peaks.astype(int)
            ax[0].scatter(
                np.linspace(0, self.file_length, len(self.filtered_data))[
                    peaks
                ],
                self.filtered_data[peaks],
                color="green",
                s=5,
            )
        fig.savefig(
            os.path.join(
                self.output_loc, "/quick_plots/" + str(self.file_loc) + ".png"
            )
        )
        plt.tight_layout()
        return fig, ax

    def output_plots(self, output_dict):
        """
        More complete plot of all events fitted, outputs 1 second trace as .png
        Events are filtered here.
        """
        pdo = pd.DataFrame()
        for x in output_dict:
            t = pd.DataFrame([output_dict[x]])
            pdo = pd.concat([pdo, t.astype(pdo.dtypes)], ignore_index=True)

        print(self.file_loc)
        filename = self.file_loc.split("\\")[-1]
        print("Creating plots...")
        for x in np.arange(0, len(self.filtered_data), self.original_sr):
            start = int(x / self.original_sr)
            end = start + 1
            z = pdo.loc[pdo.time.between(start, end, inclusive="both")]

            fig, ax = plt.subplots(figsize=(15, 5))
            ax.plot(
                np.linspace(start, end, self.original_sr),
                self.filtered_data[x : x + self.original_sr],
                linewidth=0.4,
            )
            ax.set_xlim(start, end)
            ax.set_ylim(ax.get_ylim()[0] - 10, ax.get_ylim()[1] + 10)

            # select reasonable events
            z = z.loc[np.abs(z.amplitude) > 1.5 * z.rms]
            z = z.loc[z.tau > 0.0001]
            z = z.loc[z.rise > 0.00004]
            z = z.loc[z.fit_tau < 10000]  # 10000
            z = z.loc[z.fit_tau > 0]
            z = z.loc[z.tau > z.rise * 0.9]
            z = z.loc[np.abs(z.tau_pams) > 100]

            if z.size > 0:
                for i, event in z.iterrows():
                    color = "orange"
                    if np.abs(event.amplitude) > 2 * event.rms:
                        color = "green"
                    try:
                        xs = event["time"]
                        xp = event["rise"] + xs
                        xt = event["tau"] + xp
                        yb = event["baseline"]
                        yp = event["amplitude"]
                        yt = yb + (np.exp(-1) * yp)
                        fit_x = (
                            np.array((event["fit_x"]) / self.original_sr) + xp
                        )
                    except:
                        continue
                    fit_y = event["fit_y"]
                    ax.scatter(xs, yb, color=color, s=20)
                    ax.scatter(xp, yb + yp, color=color, s=20)
                    ax.scatter(xt, yt, color=color, s=20)
                    ax.plot(fit_x, fit_y, linewidth=1, color="red")
            ax.set_title(filename)
            base = self.output_loc + "/pics/"
            if not os.path.exists(base + filename):
                os.makedirs(base + filename)
            output_dir = os.path.join(base, filename)
            fig.savefig(
                os.path.join(
                    output_dir,
                    str(int(x / self.original_sr)).zfill(3) + ".png",
                )
            )
            plt.close()
        print("Plots completed")
        return

    def run_all(self):
        self.original_data, self.sr = hf.read_abf(self.file_loc)
        self.original_sr = int(self.sr)

        print("Deconvolving data...")
        self.detection_data, self.filtered_data = hf.filtering(
            self.original_data, 200, "lowpass", self.sr
        )
        self.file_length = int(len(self.filtered_data) / self.sr)
        self.detection_data, self.sr = hf.resample_data(
            self.detection_data, 2500, self.original_sr
        )

        self.detection_data = hf.detrend_data(self.detection_data)
        template, template_y = self.create_template(template_plot=False)
        self.deconvoluted_data, self.detection_data = self.deconvolution(
            self.detection_data, template, wiener_flag=True, wiener_order=13
        )
        self.deconvoluted_data = np.abs(self.deconvoluted_data)
        self.deconvoluted_data = hf.gaussian(self.deconvoluted_data, 3)
        events, outevs, threshold = self.propagate_threshold(
            self.deconvoluted_data
        )
        self.plot = False
        if not events:
            return None
        if self.plot:
            fig, ax = self.quick_plot(threshold, outevs)
        else:
            fig, ax = None, None

        output_dict_of_dicts = {}
        print("Analysing events... " + str(len(outevs)))
        for i, x in enumerate(outevs):
            analyze = analyze_event(
                self.plot,
                fig,
                ax,
                x,
                self.filtered_data,
                self.original_sr,
                self.sr,
            )
            analyze_output = analyze.process_event()
            output_dict_of_dicts[i] = analyze_output

        # self.output_plots(output_dict_of_dicts)
        # temp, this is a big issue when outputting to run_all, as it expects all these outputs in a tuple
        output_dict_of_dicts = pd.DataFrame(output_dict_of_dicts)
        output_dict_of_dicts.to_excel(os.path.join(output_loc, "test.xlsx"))
        output_dict_of_dicts = output_dict_of_dicts.transpose()

        output_dict_of_dicts.to_excel(os.path.join(output_loc, "test.xlsx"))
        return (
            output_dict_of_dicts,
            self.filtered_data,
            self.deconvoluted_data,
            self.detection_data,
            fig,
            ax,
        )

