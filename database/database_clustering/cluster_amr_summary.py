import pandas as pd
import numpy as np
import re

def summarize_for_threshold(dbx0_content: str, dataframe: pd.DataFrame, row: int) -> pd.DataFrame:
    dbx0_by_line = dbx0_content.split("\n")
    cell_index = 0
    amr_classes = list()
    for line in dbx0_by_line[1:]:
        if (line == "") or (">" == line[0]):
            dataframe.iloc[row, cell_index] = "|".join(amr_classes)
            amr_classes.clear()
            cell_index += 1
        else:
            amr = line.split("|")[3]
            if amr not in amr_classes:
                amr_classes.append(amr)
    return(dataframe)

db90_file = open("db90.clstr", "r")
db90_content = db90_file.read()
db90_file.close
max_cluster = len(re.findall("\n>", db90_content)) + 1

dataframe = pd.DataFrame(np.full(shape=[6, max_cluster], fill_value=""), dtype=str, index=[90, 80, 70, 60, 50, 40])
dataframe = summarize_for_threshold(db90_content, dataframe, 0)

db80_file = open("db80.clstr", "r")
db80_content = db80_file.read()
db80_file.close
dataframe = summarize_for_threshold(db80_content, dataframe, 1)

db70_file = open("db70.clstr", "r")
db70_content = db70_file.read()
db70_file.close
dataframe = summarize_for_threshold(db70_content, dataframe, 2)

db60_file = open("db60.clstr", "r")
db60_content = db60_file.read()
db60_file.close
dataframe = summarize_for_threshold(db60_content, dataframe, 3)

db50_file = open("db50.clstr", "r")
db50_content = db50_file.read()
db50_file.close
dataframe = summarize_for_threshold(db50_content, dataframe, 4)

db40_file = open("db40.clstr", "r")
db40_content = db40_file.read()
db40_file.close
dataframe = summarize_for_threshold(db40_content, dataframe, 5)

dataframe.to_csv("cluster_amr_summary.csv")