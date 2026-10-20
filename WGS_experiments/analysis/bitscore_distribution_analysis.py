import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import lines

pd.options.mode.chained_assignment = None

ALIGNMENT_FILE = "X.align.daa.tsv"
MAX_CLASS_FILE = "max_class_hit_per_query.csv"
PATH_TO_SS = "deeparg_results"
PATH_TO_LS = "spades/deeparg_results"
ATTRIBUTE_DATA = "../../database/v2_feature_data.csv"

# Get list of all pair-containing queries
label_count = pd.read_csv("label_counts.tsv", sep="\t", header=0)
pair_containing_queries = label_count[
    [
        "Sample ID",
        "Alignment Identity",
        "Model",
        "Most Frequent Class",
        "DeepARG Class",
        "Query",
    ]
].drop_duplicates()
del label_count
pair_containing_queries = pair_containing_queries.sort_values(
    by=["Query", "Alignment Identity"]
)
pair_containing_queries = pair_containing_queries.set_index(
    keys=["Query", "Alignment Identity"]
)

# Keep class switch queries only
class_switch_queries = pair_containing_queries.loc[
    pair_containing_queries["Most Frequent Class"]
    != pair_containing_queries["DeepARG Class"]
]
del pair_containing_queries

# Reformat query names
class_switch_queries = class_switch_queries.reset_index()
class_switch_queries["Query"] = [
    query[13:] for query in class_switch_queries["Query"].to_list()
]
class_switch_queries = class_switch_queries.set_index(
    keys=["Sample ID", "Model", "Alignment Identity", "Query"]
)

attribute_data = pd.read_csv(ATTRIBUTE_DATA, index_col=0, header=0)
attribute_data = attribute_data[
    ["amr class", "v2-only cluster index", "superfamily(ies) id(s)"]
].rename(
    columns={
        "amr class": "class",
        "v2-only cluster index": "clstr",
        "superfamily(ies) id(s)": "super",
    }
)


def analyze_model_at_alignment_identity(model: str, al_id: int) -> pd.DataFrame:  # noqa: PLR0914
    """Analyze class switches in either SS or LS model at a specified alignment identity."""  # noqa: E501
    other_al_ids = [50, 80] if al_id == 30 else [80]  # noqa: PLR2004
    df = class_switch_queries.loc[
        class_switch_queries.index.isin(values=[model], level=1)
        & class_switch_queries.index.isin(values=[al_id], level=2)
    ]
    sample_list = df.index.get_level_values(0).drop_duplicates()
    query_sample_df: pd.DataFrame = (
        df.index
        .droplevel(2)
        .droplevel(1)
        .drop_duplicates()
        .to_frame()
        .set_index("Sample ID")
    )
    previous_align_file = ""
    alignments_df = pd.DataFrame()
    previous_first_max_class_file = ""
    first_max_class_df = pd.DataFrame()
    previous_second_max_class_file = ""
    second_max_class_df = pd.DataFrame()
    output_df = pd.DataFrame()
    for sample in sample_list:
        path_to_file = (
            f"../samples/{sample}/{PATH_TO_SS if model == 'SS' else PATH_TO_LS}"
        )
        file_to_open = (
            f"{path_to_file}/arg_alignment_identity_{al_id}/{ALIGNMENT_FILE}"
        )
        if file_to_open != previous_align_file:
            previous_align_file = file_to_open
            alignments_df = pd.read_csv(
                file_to_open,
                sep="\t",
                index_col=[0],
                usecols=[0, 1, 2, 11],
                header=None,
            ).rename(columns={1: "sequence", 2: "alignment identity", 11: "bit-score"})
        queries = query_sample_df.loc[sample]
        bitscores = alignments_df.loc[
            alignments_df.index.isin(
                queries["Query"] if len(queries) > 1 else [queries["Query"]]
            )
        ]
        bitscores["class"] = bitscores["sequence"].apply(
            lambda seq: seq.split("|")[-2]
        )
        bitscores["clstr"] = bitscores["sequence"].apply(
            lambda seq: attribute_data.at[seq, "clstr"]
        )
        bitscores["super"] = bitscores["sequence"].apply(
            lambda seq: attribute_data.at[seq, "super"]
        )
        bitscores["Most frequent"] = bitscores.index.map(
            lambda query: df.at[(sample, model, al_id, query), "Most Frequent Class"]
        )
        bitscores["DeepARG"] = bitscores.index.map(
            lambda query: df.at[(sample, model, al_id, query), "DeepARG Class"]
        )
        file_to_open = (
            f"{path_to_file}/arg_alignment_identity_{other_al_ids[0]}/{MAX_CLASS_FILE}"
        )
        if previous_first_max_class_file != file_to_open:
            previous_first_max_class_file = file_to_open
            first_max_class_df = pd.read_csv(file_to_open, header=0, index_col=0)
        bitscores[f"Have alignments at {other_al_ids[0]}"] = bitscores.index.map(
            lambda query: query in first_max_class_df.index.to_list()
        )
        bitscores[f"Most frequent {other_al_ids[0]}"] = bitscores.index.map(
            lambda query: (
                (
                    frozenset(first_max_class_df.loc[query]["deeparg hit"].to_list())
                    if np.sum(first_max_class_df.index.isin([query])) > 1
                    else frozenset({first_max_class_df.at[query, "deeparg hit"]})
                )
                if np.sum(first_max_class_df.index.isin([query])) > 0
                else frozenset()
            )
        )
        bitscores[f"Which is equal to {other_al_ids[0]}"] = bitscores.apply(
            lambda row: (
                "Both"
                if (
                    (row["Most frequent"] in row[f"Most frequent {other_al_ids[0]}"])
                    and (row["DeepARG"] in row[f"Most frequent {other_al_ids[0]}"])
                )
                else "Alignment-based"
                if (row["Most frequent"] in row[f"Most frequent {other_al_ids[0]}"])
                else "DeepARG"
                if (row["DeepARG"] in row[f"Most frequent {other_al_ids[0]}"])
                else "Neither"
            ),
            axis=1,
        )
        if len(other_al_ids) == 2:  # noqa: PLR2004
            file_to_open = (
                f"{path_to_file}/arg_alignment_identity_{other_al_ids[1]}/{MAX_CLASS_FILE}")
            if previous_second_max_class_file != file_to_open:
                previous_second_max_class_file = file_to_open
                second_max_class_df = pd.read_csv(file_to_open, header=0, index_col=0)
            bitscores[f"Have alignments at {other_al_ids[1]}"] = bitscores.index.map(
                lambda query: query in second_max_class_df.index.to_list()
            )
            bitscores[f"Most frequent {other_al_ids[1]}"] = bitscores.index.map(
                lambda query: (
                    (
                        frozenset(
                            second_max_class_df.loc[query]["deeparg hit"].to_list()
                        )
                        if np.sum(second_max_class_df.index.isin([query])) > 1
                        else frozenset({second_max_class_df.at[query, "deeparg hit"]})
                    )
                    if np.sum(second_max_class_df.index.isin([query])) > 0
                    else frozenset()
                )
            )
            bitscores[f"Which is equal to {other_al_ids[1]}"] = bitscores.apply(
                lambda row: (
                    "Both"
                    if (
                        (
                            row["Most frequent"]
                            in row[f"Most frequent {other_al_ids[1]}"]
                        )
                        and (row["DeepARG"] in row[f"Most frequent {other_al_ids[1]}"])
                    )
                    else "Alignment-based"
                    if (
                        row["Most frequent"] in row[f"Most frequent {other_al_ids[1]}"]
                    )
                    else "DeepARG"
                    if (row["DeepARG"] in row[f"Most frequent {other_al_ids[1]}"])
                    else "Neither"
                ),
                axis=1,
            )
        output_df = pd.concat([output_df, bitscores])
    hit_count_per_query = (
        output_df.reset_index(names="query").groupby("query").count()["sequence"]
    )
    output_df["hits count for query"] = output_df.index.map(
        lambda query: hit_count_per_query[query]
    )
    output_df["weight"] = output_df.apply(
        lambda row: 1 / row["hits count for query"], axis=1
    )
    min_score_per_query = (
        output_df.reset_index(names="query").groupby("query").min()["bit-score"]
    )
    output_df["min score for query"] = output_df.index.map(
        lambda query: min_score_per_query[query]
    )
    output_df["x"] = output_df.apply(
        lambda row: row["bit-score"] / row["min score for query"], axis=1
    )
    output_df[f"hue {other_al_ids[0]}"] = output_df.apply(
        lambda row: (
            "alignment-based"
            if row["class"] == row["Most frequent"]
            else "DeepARG"
            if row["class"] == row["DeepARG"]
            else f"{other_al_ids[0]}% alignment-based"
            if row["class"] in row[f"Most frequent {other_al_ids[0]}"]
            else "other"
        ),
        axis=1,
    )
    if len(other_al_ids) == 2:  # noqa: PLR2004
        output_df[f"hue {other_al_ids[1]}"] = output_df.apply(
            lambda row: (
                "alignment-based"
                if row["class"] == row["Most frequent"]
                else "DeepARG"
                if row["class"] == row["DeepARG"]
                else f"{other_al_ids[1]}% alignment-based"
                if row["class"] in row[f"Most frequent {other_al_ids[1]}"]
                else "other"
            ),
            axis=1,
        )
    return output_df


def plot_kernel_density(
    df: pd.DataFrame, name: str, hue_col: str, other_al_id: int
) -> None:
    """Plot kernel density for specified dataframe."""
    hues = [
        "alignment-based",
        "DeepARG",
        f"{other_al_id}% alignment-based",
        "other",
    ]
    palette = dict(zip(hues, sns.color_palette(n_colors=4)))
    fig, ax = plt.subplots()
    sns.kdeplot(
        x=df["x"],
        hue=df[hue_col],
        weights=df["weight"],
        palette=palette,
        hue_order=hues,
        cut=0,
        ax=ax,
    )

    ax.set_xlabel("Bit-score ratio (normalized to query's minimum hit bit-score)")
    # kde_plot.legend.set_title("Class")
    # for text, label in zip(kde_plot.legend_.texts, [
    #         "Most frequent",
    #         "DeepARG",
    #         f"{other_al_id}% most frequent",
    #         "Other",
    #     ]):
    #     text.set_text(label)
    plt.legend(
        title="Class",
        handles=[lines.Line2D([0], [0], c=color) for color in palette.values()],
        labels=[
            "Most frequent",
            "DeepARG",
            f"{other_al_id}% most frequent",
            "Other",
        ],
    )

    # Check to see if one of the elements has 0 zero variance
    zero_var = df.loc[~df[hue_col].duplicated(False)]
    if not zero_var.empty:
        # Place a dot representation of that element just above the x-axis
        ymax = ax.get_ylim()[1]
        dot_y = ymax * 0.02
        ax.scatter(
            x=zero_var["x"],
            y=dot_y,
            color=zero_var[hue_col].apply(lambda hue: palette[hue]),
            linewidth=0,
        )

    plt.savefig(fname=name)
    plt.close(fig)


for model in ["LS", "SS"]:
    for al_id in [30, 50]:
        if not pathlib.Path(f"bitscores_{model}_{al_id}.csv").exists():
            bitscores = analyze_model_at_alignment_identity(model, al_id)
            bitscores.to_csv(f"bitscores_{model}_{al_id}.csv")
        else:
            bitscores = pd.read_csv(
                f"bitscores_{model}_{al_id}.csv", header=0, index_col=0
            )
        if al_id == 30:  # noqa: PLR2004
            better_cross_diff_summary = bitscores.reset_index(names="Query")[
                [
                    "Query",
                    "Most frequent",
                    "DeepARG",
                    "Have alignments at 50",
                    "Most frequent 50",
                    "Which is equal to 50",
                    "Have alignments at 80",
                    "Most frequent 80",
                    "Which is equal to 80",
                ]
            ].drop_duplicates()
        else:
            better_cross_diff_summary = bitscores.reset_index(names="Query")[
                [
                    "Query",
                    "Most frequent",
                    "DeepARG",
                    "Have alignments at 80",
                    "Most frequent 80",
                    "Which is equal to 80",
                ]
            ].drop_duplicates()
        better_cross_diff_summary.to_csv(f"better_cross_diff_summary_{model}_{al_id}.csv")

        for other_al_id in [50, 80] if al_id == 30 else [80]:  # noqa: PLR2004
            # Plot for no other_al_id alignments
            plot_kernel_density(
                df=bitscores.loc[~bitscores[f"Have alignments at {other_al_id}"]],
                name=f"{model} {al_id}% no {other_al_id}% alignments.png",
                hue_col=f"hue {other_al_id}",
                other_al_id=other_al_id,
            )
            # Plot for neither == other_al_id alignments
            plot_kernel_density(
                df=bitscores.loc[
                    bitscores[f"Have alignments at {other_al_id}"]
                    & (bitscores[f"Which is equal to {other_al_id}"].isin(["Neither"]))
                ],
                name=f"{model} {al_id}% neither is {other_al_id}% alignments.png",
                hue_col=f"hue {other_al_id}",
                other_al_id=other_al_id,
            )
            # Plot for alignment-based == other_al_id alignments
            plot_kernel_density(
                df=bitscores.loc[
                    bitscores[f"Which is equal to {other_al_id}"].isin([
                        "Alignment-based"
                    ])
                ],
                name=f"{model} {al_id}% alignment-based is {other_al_id}% alignments.png",  # noqa: E501
                hue_col=f"hue {other_al_id}",
                other_al_id=other_al_id,
            )
            # Plot for DeepARG == other_al_id alignments
            plot_kernel_density(
                df=bitscores.loc[
                    bitscores[f"Which is equal to {other_al_id}"].isin(["DeepARG"])
                ],
                name=f"{model} {al_id}% DeepARG is {other_al_id}% alignments.png",
                hue_col=f"hue {other_al_id}",
                other_al_id=other_al_id,
            )
            # Plot for DeepARG == other_al_id alignments
            plot_kernel_density(
                df=bitscores.loc[
                    bitscores[f"Which is equal to {other_al_id}"].isin(["Both"])
                ],
                name=f"{model} {al_id}% both are {other_al_id}% alignments.png",
                hue_col=f"hue {other_al_id}",
                other_al_id=other_al_id,
            )
