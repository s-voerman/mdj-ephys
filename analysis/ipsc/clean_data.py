"""
Filter finished analyzed epsc traces based on conditions

Basically just removes all the badly fitted, small, or unrealistic events

Parameters are far from perfect still. There's still quite a few very small events.
"""

import os
import pandas as pd
import numpy as np


def data_cleaning(base_folder, output_folder):
    os.chdir(base_folder)
    files = os.listdir()
    files = [x for x in files if ".xlsx" in x]
    for f in files:
        z = pd.read_excel(os.path.join(base_folder, f))
        z = z.loc[np.abs(z.amplitude) > 1.5 * z.rms]
        z = z.loc[z.tau > 0.0001]
        z = z.loc[z.rise > 0.00004]
        z = z.loc[z.fit_tau < 10000]
        z = z.loc[z.fit_tau > 0]
        z = z.loc[z.tau > z.rise * 0.9]
        z = z.loc[np.abs(z.tau_pams) > 100]

        file = f.split(".")[0]
        print(file)
        z.to_excel(os.path.join(output_folder, str(file) + "_cleaned.xlsx"))
    return

