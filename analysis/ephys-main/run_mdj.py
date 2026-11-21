import os
import pandas as pd

import data_analysis as data_analysis
from pathlib import Path

def read_meta(meta_loc, data_loc, output_loc):
    # Setting up outputs for entire dataset
    spike_data = pd.DataFrame()
    passive_data = pd.DataFrame()
    passive_data_pt = pd.DataFrame()
    meta = pd.read_csv(os.path.join(meta_loc, "main_ephys_metadata.csv"))
    for i, row in meta.iterrows():
        file = row.file
        cell = row.cell
        seal = row.seal
        sepsc = row.sepsc
        curr = row.curr
        rheo = row.rheo
        spont = row.spont
        nucleus = row.nucleus
        p = data_analysis.EphysAnalysis(
            data_loc,
            file,
            cell,
            seal,
            sepsc,
            curr,
            rheo,
            spont,
            nucleus,
        )
        spike, passive = p.run()
        passive["AP"] = row.AP
        spike["AP"] = row.AP
        passive["nucleus"] = row.nucleus
        spike["nucleus"] = row.nucleus
        passive["GFP"] = row.GFP
        spike["GFP"] = row.GFP
        passive["exclude"] = row.Exclude
        spike["exclude"] = row.Exclude

        spike_data = pd.concat([spike_data, spike])
        passive_data = pd.concat([passive_data, passive])

    # Output as excel
    spike_data.to_csv(os.path.join(output_loc, "spike_properties.csv"))
    passive_data.to_csv(os.path.join(output_loc, "passive_properties.csv"))

    return spike_data, passive_data, passive_data_pt


if __name__ == "__main__":
    base_folder = Path.cwd()
    data_loc = os.path.join(base_folder, "data")
    meta_loc = os.path.join(base_folder, "metadata")
    output_loc = os.path.join(base_folder, "output")
    read_meta(meta_loc, data_loc, output_loc)
