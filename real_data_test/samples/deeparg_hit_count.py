import pandas as pd
import numpy as np

DEEPARG_HIT_FILE = "X.mapping.ARG"
SAMPLE_ID_FILE = "../real_samples.txt"
PATH_TO_RESULTS = {"SS": "deeparg_results", "LS": "spades/deeparg_results"}

OUTPUT = "deeparg_hit_count.tsv"

# Read content of real_samples.txt to find biosample ID
with open(SAMPLE_ID_FILE, "r") as sample_id_buf:
    sample_id_list = sample_id_buf.read().split('\n')

hit_counts_index = pd.MultiIndex.from_product(
    iterables=[sample_id_list, [30, 50, 80], ["LS", "SS"], [0]])

hit_counts = hit_counts_index.to_frame(
    index=False,
    name=["sample id", "alignment identity", "model", "hit count"],)

hit_counts["hit count"] = hit_counts.apply(
    lambda x: len(open(
        f"{x["sample id"]}/{PATH_TO_RESULTS[x["model"]]}/arg_alignment_identity_{x["alignment identity"]}/{DEEPARG_HIT_FILE}"
    ).readlines()) - 1, axis=1)

hit_counts.to_csv(OUTPUT, sep="\t")