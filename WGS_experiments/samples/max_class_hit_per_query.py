import pathlib

import numpy as np
import pandas as pd

ALIGNMENT_FILE = "X.align.daa.tsv"
PATH_TO_SS = "deeparg_results"
PATH_TO_LS = "spades/deeparg_results"
SAMPLE_LIST = "../real_samples.txt"
OUTPUT = "max_class_hit_per_query.csv"

sample_list = pathlib.Path(SAMPLE_LIST).read_text(encoding="utf-8").split("\n")
get_class = np.vectorize(lambda deeparg_hit: deeparg_hit.split("|")[-2])
for sample in sample_list:
    for al_id in [30, 50, 80]:
        for path_to_model in [PATH_TO_SS, PATH_TO_LS]:
            results_path = f"{sample}/{path_to_model}/arg_alignment_identity_{al_id}"
            results_df = pd.read_csv(
                f"{results_path}/{ALIGNMENT_FILE}", sep="\t", header=None
            )
            results_df = results_df[[0, 1, 2]].rename(
                columns={0: "query", 1: "deeparg hit"}
            )
            results_df["deeparg hit"] = get_class(results_df["deeparg hit"])
            results_df = results_df.groupby(by=["query", "deeparg hit"]).count()
            results_df = results_df.reset_index()
            same_num_of_hits = results_df.loc[
                results_df[["query", 2]].duplicated(keep=False)
            ]
            max_hit_df = results_df.loc[results_df.groupby(by=["query"]).idxmax()[2]]
            same_num_of_hits_as_max = same_num_of_hits.loc[
                same_num_of_hits.set_index(["query", 2]).index.isin(
                    max_hit_df.set_index(["query", 2]).index.to_list()
                )
            ]
            pd.concat([
                max_hit_df,
                same_num_of_hits_as_max,
            ]).drop_duplicates().set_index("query").to_csv(f"{results_path}/{OUTPUT}")
