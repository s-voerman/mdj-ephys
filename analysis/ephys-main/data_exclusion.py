import os

import pandas as pd
import numpy as np
from pathlib import Path


def exclude(passive_data, spike_data):
    """
    Make a list of all cells that need to be excluded
    + reason for exclusion (injections + ephys)
    """
    # Output setup
    exclusion_list = []
    exclusion_df = pd.DataFrame()
    td = {}

    # Cells excluded based on injections
    bad_inj = passive_data.loc[passive_data["exclude"] == "Y"]
    bad_inj = list(bad_inj["Cell"])
    exclusion_list.append(bad_inj)
    exclusion_list = [x for xs in exclusion_list for x in xs]
    for x in exclusion_list:
        td["Cell"] = x
        td["Reason"] = "IO Injection"
        exclusion_df = pd.concat([exclusion_df, pd.DataFrame([td])])

    # Cells excluded based on ephys
    vars = ["Rs", "Rs/Ri", "Spike_amplitude", "Ri"]
    criteria = [30, 0.2, 80, 2000]
    conditions = [">", ">", "<", ">"]

    for i in zip(vars, criteria, conditions):
        if i[2] == ">":
            bad_cells = passive_data.loc[passive_data[i[0]] > i[1]]
            n = bad_cells.Cell.to_list()
            reason = i[0]
        elif i[2] == "<":
            bad_cells = passive_data.loc[passive_data[i[0]] < i[1]]
            n = bad_cells.Cell.to_list()
            reason = i[0]
        for r, v in enumerate(n):
            cell = n[r]
            exclusion_list.append(cell)
            td["Cell"] = cell
            td["Reason"] = reason
            exclusion_df = pd.concat([exclusion_df, pd.DataFrame([td])])

    # Apply exclusion to passive and spike data
    passive_data_f = passive_data
    spike_data_f = spike_data
    for t in exclusion_df.Cell.unique():
        passive_data_f = passive_data_f.loc[passive_data_f["Cell"] != t]
        spike_data_f = spike_data_f.loc[spike_data_f["Cell"] != t]
    passive_data_f = passive_data_f.drop(["Unnamed: 0"], axis=1)
    spike_data_f = spike_data_f.drop(["Unnamed: 0"], axis=1)

    # Exports
    passive_data_f.to_excel(os.path.join(output_loc, "passive_data_f.xlsx"))
    spike_data_f.to_excel(os.path.join(output_loc, "spike_data_f.xlsx"))
    exclusion_df.to_excel(os.path.join(output_loc, "excluded_cells.xlsx"))
    return


if __name__ == "__main__":
    base_folder = Path(__file__).resolve().parent.parent
    data_loc = os.path.join(base_folder, "data")
    meta_loc = os.path.join(base_folder, "metadata")
    output_loc = os.path.join(base_folder, "output")
    passive_data = pd.read_excel(
        os.path.join(output_loc, "passive_properties.csv")
    )
    spike_data = pd.read_excel(
        os.path.join(output_loc, "spike_properties.csv")
    )
    exclude(passive_data, spike_data)
