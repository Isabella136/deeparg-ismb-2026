import pathlib

import pandas as pd

DIAMOND_HIT_FILE = "X.align.daa.tsv"
SAMPLE_ID_FILE = "../real_samples.txt"
PATH_TO_RESULTS = {"SS": "deeparg_results", "LS": "spades/deeparg_results"}
FEATURE_DATA = "../../database/v2_feature_data.csv"

OUTPUT = "hit_classes.tsv"

feature_classes = pd.read_csv(
    FEATURE_DATA, header=0, index_col=0, usecols=[0, 1]
)
sample_id_list = (
    pathlib.Path(SAMPLE_ID_FILE).read_text(encoding="utf-8").split("\n")
)

hit_classes = set()

for sample_id in sample_id_list:
    for model in ["SS", "LS"]:
        hit_classes.update(
            pd
            .read_csv(
                f"{sample_id}/{PATH_TO_RESULTS[model]}/arg_alignment_identity_30/{DIAMOND_HIT_FILE}",
                sep="\t",
                header=None,
                usecols=[1]
            )
            .drop_duplicates()
            .apply(lambda x: feature_classes.at[x[1], "amr class"], axis=1)
            .to_list()
        )

pathlib.Path(OUTPUT).write_text("\n".join(hit_classes), encoding="utf-8")
