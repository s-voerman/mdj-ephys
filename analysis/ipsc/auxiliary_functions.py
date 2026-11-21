import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd


def convert_to_excel(all_files_output):
    write = pd.ExcelWriter("./output1.xlsx")
    for x in all_files_output.keys():
        pdo = pd.DataFrame()
        for e in all_files_output[x]:
            pdo = pdo.append(e, ignore_index=True)
        pdo = pdo[["event", "time", "amplitude", "decay", "rise"]]
        pdo["date"] = x.split(" ")[0]
        pdo.to_excel(write, sheet_name=x)
    write.save()


def convert_to_pandas(input_data):
    pdo = pd.DataFrame()
    for x in input_data.keys():
        for y in input_data[x]:
            y = input_data[x][y]
            p = pd.DataFrame(y)
            meta = x.split("/")[-1].split(" ")
            date, cell = meta[0], meta[1]
            p["date"] = date
            p["cell"] = cell
            p["id"] = date + "_" + cell
            pdo = pdo.append(p, ignore_index=True)
    pdo = pdo.reset_index(drop=True)
    return pdo
