import os
import datetime
import time

import pandas as pd
from abf_v4 import ProcessABF
from pathlib import Path

from auxiliary_functions import convert_to_excel, convert_to_pandas


outliers = []


# This function should inherit the plot value from ABF_v4, but it works now, so lets not mess with it
def go(output_loc, outliers=outliers, files=None, plot=False):
    """
    Runs everything
    """
    for x in files:
        all_files_output = {}
        if x in outliers:
            continue

        print(datetime.datetime.now(), x)

        processing = ProcessABF(x, output_loc, plot=False)
        out = processing.run_all()

        if out == None:
            continue
        elif plot:
            all_files_output[x] = out[0]
            orig_data = out[1]
            decon_data = out[2]
            detection_data = out[3]
            fig = out[4]
            ax = out[5]
        else:
            all_files_output[x] = out

        # save space in pandas, remove fit data.
        for ev in all_files_output[x]:
            all_files_output[x][ev]["fit_x"] = None
            all_files_output[x][ev]["fit_y"] = None

        # output to an excel file.
        pdo = pd.DataFrame()
        for x in all_files_output:
            o = all_files_output[x]
            for y in o:
                y = o[y]
                y["cell"] = x
                pdo = pd.concat(
                    [pdo, pd.DataFrame([y]).astype(pdo.dtypes)],
                    ignore_index=True,
                )
        filename = x.split("\\")[-1]
        print(output_loc)
        pdo.to_excel(os.path.join(output_loc, str(filename) + ".xlsx"))

    try:
        return all_files_output, orig_data, decon_data, fig, ax, processing
    except NameError:
        return all_files_output, orig_data, decon_data


def generate_filenames(meta, data_loc):
    for index, row in meta.iterrows():
        seal = "-".join(
            [" ".join([str(row["file"]), str(row["cell"])]), str(row["seal"])]
        )
        seal = os.path.join(data_loc, seal) + ".ABF"
        epsc = "-".join(
            [" ".join([str(row["file"]), str(row["cell"])]), str(row["ipsc"])]
        )
        epsc = os.path.join(data_loc, epsc) + ".ABF"
        meta.loc[index, "seal"] = seal
        meta.loc[index, "epsc"] = epsc
    return meta


if __name__ == "__main__":
    base_folder = Path(__file__).resolve().parent.parent.parent
    data_loc = os.path.join(base_folder, "data")
    meta_loc = os.path.join(base_folder, "metadata")
    output_loc = os.path.join(base_folder, "output")
    meta = pd.read_excel(os.path.join(meta_loc, "ipsc.csv"))
    meta = generate_filenames(meta, data_loc)
    files = meta["epsc"]
    print(files)
    out = go(output_loc, files=files, plot=True)
