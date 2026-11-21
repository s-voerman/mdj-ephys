from neo.rawio import AxonRawIO
import numpy as np

def read_abf(file_loc):
    """
    Reads .ABF V1 files and returns data and samping rate

    Parameters:
    file_loc: str
    """
    a = AxonRawIO(file_loc)
    a.parse_header()
    sr = int(a.header["signal_channels"][0][2])
    data = a.get_analogsignal_chunk() * a.header["signal_channels"][0][5]
    return data[:, 0], sr


def read_abf_v2(file_loc):
    """
    Reads .ABF V2+ files and returns data and samping rate, for all sweeps

    Parameters:
    file_loc: str
    """
    a = AxonRawIO(file_loc)
    a.parse_header()
    sweeps = a.header["nb_segment"][0]  # Check at some point if this is right
    sr = int(a.header["signal_channels"][0][2])
    data = np.array([])
    for i in np.arange(0, sweeps, 1):
        sweep_data = (
            a.get_analogsignal_chunk(seg_index=i, channel_indexes=[0])
            * a.header["signal_channels"][0][5]
        )
        data = np.append(data, sweep_data[:, 0])
    return data, sr


def edge_detection(data, threshold):
    """
    Simple function for detecting signals over a certain threshold.
    Returns index when a value in array 'data' exceeds 'threshold'

    Parameters
    data: numpy array
    threshold: int or float
    """
    signal = data
    thresholded_data = signal > threshold
    threshold_edges = np.convolve([1, -1], thresholded_data, mode="same")
    threshold_edge_indices = np.where(threshold_edges == 1)[0]
    return threshold_edge_indices


def exp_func(x, a, tau, c):
    return a * -np.exp(-x / tau) + c


def mono_exp(x, m, tau, b):
    return m * np.exp(-x / tau) + b


def double_exp(x, a, tau, c, tau2):
    return (a * np.exp(-x / tau) + c) + (a * np.exp(-x / tau2) + c)

