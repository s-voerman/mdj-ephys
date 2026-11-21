import os

import pandas as pd
import numpy as np
import pathlib as Path


def exclude(seals, opto):
    """
    Make a list of all cells that need to be excluded
    + reason for exclusion (injections + ephys)
    """
    # Output setup
    exclusion_list = []
    exclusion_df = pd.DataFrame()
    td = {}

    # Cells excluded based on ephys
    vars = ["Rs", "Rs/Ri"]
    criteria = [30, 0.2]
    conditions = [">", ">"]

    for i in zip(vars, criteria, conditions):
        if i[2] == ">":
            bad_cells = seals.loc[seals[i[0]] > i[1]]
            n = bad_cells.Cell.to_list()
            reason = i[0]
        elif i[2] == "<":
            bad_cells = seals.loc[seals[i[0]] < i[1]]
            n = bad_cells.Cell.to_list()
            reason = i[0]
        for r, v in enumerate(n):
            cell = n[r]
            exclusion_list.append(cell)
            td["Cell"] = cell
            td["Reason"] = reason
            exclusion_df = pd.concat([exclusion_df, pd.DataFrame([td])])

    # Apply exclusion to passive and spike data
    seals_f = seals
    opto_f = opto
    for t in exclusion_df.Cell.unique():
        seals_f = seals_f.loc[seals_f["Cell"] != t]
        opto_f = opto_f.loc[opto_f["cell"] != t]
    seals_f = seals_f.drop(["Unnamed: 0"], axis=1)
    opto_f = opto_f.drop(["Unnamed: 0"], axis=1)

    # Exports
    seals_f.to_excel(os.path.join(output_loc, "seals_f.xlsx"))
    opto_f.to_excel(os.path.join(output_loc, "opto_f.xlsx"))
    exclusion_df.to_excel(os.path.join(output_loc, "excluded_cells.xlsx"))
    return


if __name__ == "__main__":
    base_folder = Path(__file__).resolve().parent.parent.parent
    output_loc = os.path.join(base_folder, "output/optogenetics")

    seals = pd.read_excel(os.path.join(output_loc, "opto_output_seal.xlsx"))
    opto = pd.read_excel(os.path.join(output_loc, "opto_output_gfppos.xlsx"))
    exclude(seals, opto)
