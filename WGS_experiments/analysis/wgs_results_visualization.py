import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Circle

sys.path.append("../../")

from utils.curved_text import CurvedText

HIT_COUNT = "../samples/deeparg_hit_count.tsv"
INPUT_SEQ_COUNT = "../samples/sequence_count.tsv"
FEATURE_DATA = "../../database/v2_feature_data.csv"

mpl.rcParams["mathtext.fontset"] = "dejavusans"
mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["svg.fonttype"] = "none"

# Get alignment label count
label_df = pd.read_csv("label_counts.tsv", sep="\t", header=0).rename(
    columns={"bitscore": "count"}
)

# Populate hit_count_amr_df for bar graphs and sort it by run
hit_count_amr_df = pd.read_csv(HIT_COUNT, sep="\t", header=0, index_col=0)
hit_count_amr_df = hit_count_amr_df.sort_values(
    by=["model", "sample id", "alignment identity"]
)

# Get percentages of sequences that had an alignment hit
input_seq_count_df = pd.read_csv(
    INPUT_SEQ_COUNT, sep="\t", header=0, index_col=0
)
hit_count_amr_df["alignment hit percentage"] = hit_count_amr_df.apply(
    lambda x: (
        float(x["diamond hit count"])
        / float(
            input_seq_count_df.at[
                x["sample id"], ("Read" if x["model"] == "SS" else "ORF")
            ]
        )
    ),
    axis=1
)

# Get percentages of alignment hits that had a DeepARG prediction
hit_count_amr_df["deeparg prediction percentage"] = hit_count_amr_df.apply(
    lambda x: float(x["deeparg hit count"]) / float(x["diamond hit count"]),
    axis=1,
)

# Get percentage of DeepARG predicitions that differs from most frequent hit
diff_prediction_count = (
    label_df
    .loc[(label_df["Is DeepARG Class"] != label_df["Is Most Frequent Class"])]
    .drop_duplicates(
        subset=["Sample ID", "Alignment Identity", "Model", "Query"]
    )
    .value_counts(subset=["Sample ID", "Alignment Identity", "Model"])
)
hit_count_amr_df["diff freq percentage"] = hit_count_amr_df.apply(
    lambda x: (
        float(
            diff_prediction_count.at[
                (x["sample id"], x["alignment identity"], x["model"])
            ]
            if (x["sample id"], x["alignment identity"], x["model"])
            in diff_prediction_count.index
            else 0
        )
        / float(x["deeparg hit count"])
    ),
    axis=1,
)

cb_palette = sns.color_palette("colorblind")

# Make figure 3
plt.figure(figsize=(15, 30))
perc_freq_diff_ax = plt.axes((0.15, 0.07, 0.8, 0.28))
perc_pred_ax = plt.axes((0.15, 0.385, 0.8, 0.28), sharex=perc_freq_diff_ax)
perc_hit_ax = plt.axes((0.15, 0.7, 0.8, 0.28), sharex=perc_freq_diff_ax)


# Draw figure 3A: percentages of sequences that had an alignment hit, averaged
# for each model–alignment identity run condition
sns.barplot(
    data=hit_count_amr_df,
    x="alignment identity",
    y="alignment hit percentage",
    hue="model",
    errorbar=("pi", 75),
    err_kws={"linewidth": 5, "color": "black"},
    ax=perc_hit_ax,
    palette=["#52C4FF", "#9B0009"]
)

perc_hit_ax.legend(
    handles=perc_hit_ax.get_legend_handles_labels()[0],
    labels=["Long Sequence Model", "Short Sequence Model"],
    fontsize=30
)

# Since perc_hit_ax plot is vertically aligned to perc_freq_diff_ax plot, we
# don't need to label perc_hit_ax x-axis, nor do we need tick marks
perc_hit_ax.set_xlabel(xlabel="")
perc_hit_ax.tick_params(axis="x", labelbottom=False)

# y-axis is average percentages from 0 to 5. Axis title is set later in the code
perc_hit_ax.set_ylabel("")
perc_hit_ax.set_yticks(
    ticks=[0, 0.01, 0.02, 0.03, 0.04, 0.05],
    labels=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
    fontsize=30,
)
perc_hit_ax.set_ylim(top=0.05)

# Set title
plt.gcf().text(0.02, 0.98, "A", fontsize=40, va="center", weight="bold")

# Styling the graph a bit more
perc_hit_ax.tick_params(length=0, pad=8)
perc_hit_ax.grid(visible=True, which="major", axis="y", color="0.75")

# Draw figure 3B: percentages of alignment hits with a deeparg prediction,
# averaged for each model–alignment identity run condition
sns.barplot(
    data=hit_count_amr_df,
    x="alignment identity",
    y="deeparg prediction percentage",
    hue="model",
    errorbar=("pi", 75),
    err_kws={"linewidth": 5, "color": "black"},
    ax=perc_pred_ax,
    palette=["#52C4FF", "#9B0009"],
    legend=False
)

# Since perc_pred_ax plot is vertically aligned to perc_freq_diff_ax plot, we
# don't need to label perc_pred_ax x-axis, nor do we need tick marks
perc_pred_ax.set_xlabel(xlabel="")
perc_pred_ax.tick_params(axis="x", labelbottom=False)

# y-axis is average percentages from 0 to 100. Axis title set later in the code
perc_pred_ax.set_ylabel("")
perc_pred_ax.set_yticks(
    ticks=[0, 0.2, 0.4, 0.6, 0.8, 1],
    labels=["0.0", "20.0", "40.0", "60.0", "80.0", "100.0"],
    fontsize=30,
)

# Set title
plt.gcf().text(0.02, 0.665, "B", fontsize=40, va="center", weight="bold")

# Styling the graph a bit more
perc_pred_ax.tick_params(length=0, pad=8)
perc_pred_ax.grid(visible=True, which="major", axis="y", color="0.75")

# Draw figure 3C: percentages of deeparg prediction that differs from most
# frequent alignment hit, averaged for each model–alignment identity run
# condition
sns.barplot(
    data=hit_count_amr_df,
    x="alignment identity",
    y="diff freq percentage",
    hue="model",
    errorbar=("pi", 75),
    err_kws={"linewidth": 5, "color": "black"},
    ax=perc_freq_diff_ax,
    palette=["#52C4FF", "#9B0009"],
    legend=False
)

# x-axis is the list of model–alignment identity run condition
perc_freq_diff_ax.set_xlabel(
    xlabel="Alignment Identity Threshold",
    fontsize=35,
    loc="center",
    labelpad=10,
)
perc_freq_diff_ax.set_xticks(
    ticks=range(3),
    labels=["30%", "50%", "80%"],
    fontsize=30,
)

# y-axis is average percentages from 0 to 10. Axis title is set later in the code  # noqa: E501
perc_freq_diff_ax.set_ylabel("")
perc_freq_diff_ax.set_ylim(top=0.08)
perc_freq_diff_ax.set_yticks(
    ticks=[0, 0.02, 0.04, 0.06, 0.08],
    labels=[0.0, 2.0, 4.0, 6.0, 8.0], fontsize=30)

# Set title
plt.gcf().text(0.02, 0.35, "C", fontsize=40, va="center", weight="bold")

# Styling the graph a bit more
perc_freq_diff_ax.tick_params(length=0, pad=8)
perc_freq_diff_ax.grid(visible=True, which="major", axis="y", color="0.75")

# We want y-axis labels for all graphs to be vertically aligned
plt.gcf().text(
    0.02, 0.21, "Average Percentage", rotation=90, fontsize=35, va="center"
)
plt.gcf().text(
    0.02, 0.525, "Average Percentage", rotation=90, fontsize=35, va="center"
)
plt.gcf().text(
    0.02, 0.84, "Average Percentage", rotation=90, fontsize=35, va="center"
)

# Save plot and be done with it
plt.savefig("barplot.svg")

# Only mapping a specific identity and model
model = sys.argv[1]
ident = int(sys.argv[2])
label_df = label_df.loc[
    (label_df["Model"] == model) & (label_df["Alignment Identity"] == ident)
].reset_index(drop=True)

# Rename AMR columns based on abbrev
amr_abbrev = pd.read_csv("amr_abbrev.csv", header=None, index_col=0)
amr_abbrev.at["aminoglycoside:aminocoumarin", 1] = "AG/\nAC"
label_df["amr"] = label_df["amr"].apply(lambda x: amr_abbrev.at[x, 1])
label_df["Most Frequent Class"] = label_df["Most Frequent Class"].apply(
    lambda x: amr_abbrev.at[x, 1]
)
label_df["DeepARG Class"] = label_df["DeepARG Class"].apply(
    lambda x: amr_abbrev.at[x, 1]
)

# Retrieve feature data, and replace amr class with abbreviation
amr_abbrev = pd.read_csv(
    "../../database/amr_abbrev.csv", header=None, index_col=0
)
amr_abbrev.at["aminoglycoside:aminocoumarin", 1] = "AG/\nAC"
feature_data = pd.read_csv(FEATURE_DATA, header=0, index_col=0)
feature_data["amr class"] = feature_data.apply(
    lambda x: amr_abbrev.at[x["amr class"], 1], axis=1
)

# Get individual clstr|class for DeepARG database from feature data
clstr_class_df = (
    feature_data
    .reset_index(drop=True)
    [["amr class", "v2-only cluster index"]]
    .drop_duplicates()
)

# Using code from good old stackoverflow:
# https://stackoverflow.com/questions/71688904/dealing-with-multiple-values-in-pandas-dataframe-cell
# Get super|class for DeepARG database from feature data
# (Counting each super in multi-super features)

indiv_super_class_df = feature_data.reset_index()[
    ["index", "amr class", "superfamily(ies) id(s)"]
].melt("index")
indiv_super_class_df["value"] = indiv_super_class_df[
    "value"
].str.split("$")
indiv_super_class_df = indiv_super_class_df.explode("value")
corresponding_class_df = indiv_super_class_df.loc[
    indiv_super_class_df["variable"] == "amr class"
]
corresponding_class_df = corresponding_class_df[["index", "value"]].set_index(
    "index"
)
indiv_super_class_df = indiv_super_class_df.loc[
    indiv_super_class_df["variable"] == "superfamily(ies) id(s)"
]
indiv_super_class_df = pd.DataFrame(
    indiv_super_class_df[["index", "value"]]
    .apply(
        lambda x: pd.Series({
            "amr class": corresponding_class_df.at[x["index"], "value"],
            "superfamily": x["value"],
        }),
        axis=1,
    )
    .drop_duplicates()
)

# Find classes that share clusters with another class
classes_per_clstr_count_df = (
    clstr_class_df["v2-only cluster index"]
    .value_counts(dropna=True, sort=False)
    .reset_index()
    .set_axis(["clstr", "count"], axis=1)
)
multi_class_clstr = classes_per_clstr_count_df.loc[
    classes_per_clstr_count_df["count"] > 1
]["clstr"].to_list()
shared_clstr_class = clstr_class_df.loc[
    clstr_class_df["v2-only cluster index"].isin(multi_class_clstr)
]["amr class"].drop_duplicates()

# Find classes that share superfamilies with another class
classes_per_super_count_df = (
    indiv_super_class_df["superfamily"]
    .value_counts(dropna=True, sort=False)
    .reset_index()
    .set_axis(["super", "count"], axis=1)
)
multi_class_super = classes_per_super_count_df.loc[
    classes_per_super_count_df["count"] > 1
]["super"].to_list()
shared_super_class = indiv_super_class_df.loc[
    indiv_super_class_df["superfamily"].isin(multi_class_super)
]["amr class"].drop_duplicates()

# Those are the classes that share superfamilies or clusters; will be colored
# differently than other classes
sharing_classes = set(shared_clstr_class).union(shared_super_class)

# Must get the queries of interest:
pair_label_df = label_df.loc[
    label_df["Is Most Frequent Class"] != label_df["Is DeepARG Class"]
][
    ["Sample ID", "Query", "Most Frequent Class", "DeepARG Class"]
].drop_duplicates()

# For the heatmaps, we will make X-axis be Diamond and y-axis be DeepARG.
heatmap_y_axis = pd.Index(
    pair_label_df["Most Frequent Class"]
    .drop_duplicates()
    .sort_values()
    .to_list()
)
heatmap_x_axis = pd.Index(
    pair_label_df["DeepARG Class"].drop_duplicates().sort_values().to_list()
)

# Aggregate pair switch
pair_df = (
    pair_label_df
    .groupby(by=["Most Frequent Class", "DeepARG Class"])[["Query"]]
    .count()
)

# Create dataframe to count how many times switch from x to y happens.
# Only needs to be filled up once (can be done at alignment)
switch_df = pd.DataFrame(
        data=np.full(
            shape=(heatmap_x_axis.shape[0], heatmap_y_axis.shape[0]),
            fill_value=np.nan,
        ),
        index=heatmap_x_axis,
        columns=heatmap_y_axis,
    )

for graph_type in ["share", "imbalance", "frequency"]:
    # See if switches coincide with label sharing
    if graph_type == "share":
        # Create dataframe to keep track of switches that coincide with label
        # sharing
        share_df = pd.DataFrame(
            data=np.full(
                shape=(heatmap_x_axis.shape[0], heatmap_y_axis.shape[0]),
                fill_value=0.0,
            ),
            index=heatmap_x_axis,
            columns=heatmap_y_axis,
        )

        # Get class y labels for queries where class y is DeepARG prediction
        y_alignments = (
            label_df
            .loc[
                (
                    label_df["Is DeepARG Class"]
                    & ~label_df["Is Most Frequent Class"]
                )
            ][
                [
                    "Sample ID",
                    "Most Frequent Class",
                    "DeepARG Class",
                    "Query",
                    "clstr",
                    "super",
                ]
            ]
            .drop_duplicates()
            .set_index(
                ["Sample ID", "Most Frequent Class", "DeepARG Class", "Query"],
                append=True,
            )
        )

        # Get class x labels for queries where class x is most freq class hit
        x_alignments = (
            label_df
            .loc[
                (
                    ~label_df["Is DeepARG Class"]
                    & label_df["Is Most Frequent Class"]
                )
            ][
                [
                    "Sample ID",
                    "Most Frequent Class",
                    "DeepARG Class",
                    "Query",
                    "clstr",
                    "super",
                ]
            ]
            .drop_duplicates()
            .set_index(
                ["Sample ID", "Most Frequent Class", "DeepARG Class", "Query"],
                append=True,
            )
        )

        for row in pair_label_df.iterrows():
            # For x- and y-axis
            deeparg = row[1]["DeepARG Class"]
            diamond = row[1]["Most Frequent Class"]

            # Check for annotation sharing
            sample = row[1]["Sample ID"]
            query = row[1]["Query"]
            x_query_alignment = x_alignments.loc[
                (slice(None), sample, diamond, deeparg, query)
            ][["clstr", "super"]]
            y_query_alignment = y_alignments.loc[
                (slice(None), sample, diamond, deeparg, query)
            ][["clstr", "super"]]
            if (
                True
                in (
                    x_query_alignment["clstr"].isin(y_query_alignment["clstr"])
                    | x_query_alignment["super"].isin(
                        y_query_alignment["super"]
                    )
                ).to_list()
            ):
                share_df.at[deeparg, diamond] += 1.0

        # Now let's get the heatmap values
        heatmap_df = share_df.div(switch_df, fill_value=np.nan)

    # See if database skews toward Diamond or DeepARG
    elif graph_type == "imbalance":
        # Create dataframe to keep track of class with biggest superfamily for
        # each class switch
        biggest_df = pd.DataFrame(
            data=np.full(
                shape=(heatmap_x_axis.shape[0], heatmap_y_axis.shape[0]),
                fill_value=0,
            ),
            index=heatmap_x_axis,
            columns=heatmap_y_axis,
        )

        # Dataframe of whether class y has biggest superfamily
        y_alignments = (
            label_df
            .loc[
                (
                    label_df["Is DeepARG Class"]
                    & ~label_df["Is Most Frequent Class"]
                )
            ][
                [
                    "Sample ID",
                    "Most Frequent Class",
                    "DeepARG Class",
                    "Query",
                    "Is Class with Biggest Superfam",
                ]
            ]
            .drop_duplicates()
            .set_index(
                ["Sample ID", "Most Frequent Class", "DeepARG Class", "Query"]
            )
        )

        # Dataframe of whether class x has biggest superfamily
        x_alignments = (
            label_df
            .loc[
                (
                    ~label_df["Is DeepARG Class"]
                    & label_df["Is Most Frequent Class"]
                )
            ][
                [
                    "Sample ID",
                    "Most Frequent Class",
                    "DeepARG Class",
                    "Query",
                    "Is Class with Biggest Superfam",
                ]
            ]
            .drop_duplicates()
            .set_index(
                ["Sample ID", "Most Frequent Class", "DeepARG Class", "Query"]
            )
        )

        for row in pair_label_df.iterrows():
            # For x- and y-axis
            deeparg = row[1]["DeepARG Class"]
            diamond = row[1]["Most Frequent Class"]

            # Check for biggest superfamily
            sample = row[1]["Sample ID"]
            query = row[1]["Query"]
            biggest = 0
            y_query_alignment = y_alignments.at[
                (sample, diamond, deeparg, query),
                "Is Class with Biggest Superfam",
            ]
            x_query_alignment = x_alignments.at[
                (sample, diamond, deeparg, query),
                "Is Class with Biggest Superfam",
            ]
            if y_query_alignment:
                biggest = 1
            elif x_query_alignment:
                biggest = -1
            biggest_df.at[deeparg, diamond] += float(biggest)

        # Now let's get the heatmap values
        heatmap_df = biggest_df.div(switch_df, fill_value=np.nan)

    # Check frequency of queries where switch from x to y happens.
    else:
        # Calculate the number of hits to amr class for a given query
        amr_count_df = (
            label_df[["Sample ID", "Query", "amr", "count"]]
            .groupby(["Sample ID", "Query", "amr"])
            .sum()
        )

        # Get all queries where x = most frequent class hits and y = another hit
        x_is_best_y_in_alignment = label_df.loc[
            pd.merge(  # noqa: PD015
                left=label_df[["Most Frequent Class", "amr"]].rename(
                    columns={"amr": "DeepARG Class"}
                ),
                right=pair_label_df[
                    ["Most Frequent Class", "DeepARG Class"]
                ].drop_duplicates(),
                how="left",
                indicator="exist",
            )["exist"]
            == "both"
        ]

        # Calculate the number of queries where x = most frequent class hits and
        # y = another hit
        query_counts = x_is_best_y_in_alignment[
            ["Most Frequent Class", "Query", "amr"]
        ].drop_duplicates()
        query_counts = query_counts.groupby(
            ["Most Frequent Class", "amr"]).count()

        for row in pair_df.iterrows():
            # For x- and y-axis
            deeparg_class = row[0][1]
            diamond_class = row[0][0]

            # For switch_df
            switch = row[1]["Query"]

            # For switch_freq_df
            queries = query_counts.at[(diamond_class, deeparg_class), "Query"]

            # Insert in df
            switch_df.at[deeparg_class, diamond_class] = float(
                switch / queries)

        # Now let's get the heatmap values
        heatmap_df = switch_df

    fig = plt.figure(figsize=[18, 18])
    ax = fig.add_axes((0.08, 0.08, 0.84, 0.84))

    # Create graph objects along with nodes and edges
    G = nx.MultiDiGraph()
    G.add_nodes_from(set(heatmap_x_axis.to_list()).union(set(heatmap_y_axis)))
    G.add_weighted_edges_from([
        (y, x, round(heatmap_df.at[x, y], 3))
        for x in heatmap_x_axis.to_list()
        for y in heatmap_y_axis.to_list()
        if not np.isnan(heatmap_df.at[x, y])
    ])

    # Specify a specific position for each node. Based on class specificity.
    partition_super = 18
    slice_super = 360.0 / partition_super
    pos = {
        "AG/\nAC": (
            np.sin(np.deg2rad(5 * slice_super)) * 1.25,
            np.cos(np.deg2rad(5 * slice_super)) * 1.25,
        ),
        "AG": (
            np.sin(np.deg2rad(5 * slice_super)) * 2.00,
            np.cos(np.deg2rad(5 * slice_super)) * 2.00,
        ),
        "BCM": (
            np.sin(np.deg2rad(7 * slice_super)) * 2.00,
            np.cos(np.deg2rad(7 * slice_super)) * 2.00,
        ),
        "BL": (
            np.sin(np.deg2rad(18 * slice_super)) * 1.25,
            np.cos(np.deg2rad(18 * slice_super)) * 1.25,
        ),
        "DAP": (
            np.sin(np.deg2rad(6 * slice_super)) * 2.00,
            np.cos(np.deg2rad(6 * slice_super)) * 2.00,
        ),
        "FFA": (
            np.sin(np.deg2rad(8 * slice_super)) * 2.00,
            np.cos(np.deg2rad(8 * slice_super)) * 2.00,
        ),
        "FQ": (
            np.sin(np.deg2rad(3 * slice_super)) * 2.00,
            np.cos(np.deg2rad(3 * slice_super)) * 2.00,
        ),
        "GlyP": (
            np.sin(np.deg2rad(15 * slice_super)) * 2.00,
            np.cos(np.deg2rad(15 * slice_super)) * 2.00,
        ),
        "MDR": (
            np.sin(np.deg2rad(75)) * 0.50,
            np.cos(np.deg2rad(75)) * 0.50
        ),
        "MLS": (
            np.sin(np.deg2rad(17 * slice_super)) * 1.25,
            np.cos(np.deg2rad(17 * slice_super)) * 1.25,
        ),
        "NI": (
            np.sin(np.deg2rad(11 * slice_super)) * 2.00,
            np.cos(np.deg2rad(11 * slice_super)) * 2.00,
        ),
        "NUC": (
            np.sin(np.deg2rad(10 * slice_super)) * 2.00,
            np.cos(np.deg2rad(10 * slice_super)) * 2.00,
        ),
        "OXA": (
            np.sin(np.deg2rad(14 * slice_super)) * 2.00,
            np.cos(np.deg2rad(14 * slice_super)) * 2.00,
        ),
        "PEP": (
            np.sin(np.deg2rad(1 * slice_super)) * 2.00,
            np.cos(np.deg2rad(1 * slice_super)) * 2.00,
        ),
        "PHE": (
            np.sin(np.deg2rad(13 * slice_super)) * 2.00,
            np.cos(np.deg2rad(13 * slice_super)) * 2.00,
        ),
        "PLM": (
            np.sin(np.deg2rad(16 * slice_super)) * 2.00,
            np.cos(np.deg2rad(16 * slice_super)) * 2.00,
        ),
        "PMX": (
            np.sin(np.deg2rad(1 * slice_super)) * 2.75,
            np.cos(np.deg2rad(1 * slice_super)) * 2.75,
        ),
        "SUL": (
            np.sin(np.deg2rad(4 * slice_super)) * 2.00,
            np.cos(np.deg2rad(4 * slice_super)) * 2.00,
        ),
        "TET": (
            np.sin(np.deg2rad(9 * slice_super)) * 2.00,
            np.cos(np.deg2rad(9 * slice_super)) * 2.00,
        ),
        "TCM": (
            np.sin(np.deg2rad(12 * slice_super)) * 2.00,
            np.cos(np.deg2rad(12 * slice_super)) * 2.00,
        ),
        "TRI": (
            np.sin(np.deg2rad(2 * slice_super)) * 2.75,
            np.cos(np.deg2rad(2 * slice_super)) * 2.75,
        ),
        "UNC": (
            np.sin(np.deg2rad(-105)) * 0.50,
            np.cos(np.deg2rad(-105)) * 0.50
        ),
    }

    # Specify class specificity with circles
    ax.add_collection(
        PatchCollection(
            [
                Circle(xy=(0, 0), radius=3.125),
                Circle(xy=(0, 0), radius=2.375),
                Circle(xy=(0, 0), radius=1.625),
                Circle(xy=(0, 0), radius=0.875),
            ],
            facecolors=["darkgrey", "white", "white", "white"],
            alpha=[0.75, 0.4, 0.4, 1],
        )
    )

    # Separate unrelated classes with lines
    ax.add_collection(
        LineCollection(
            [
                [
                    (
                        np.sin(np.deg2rad(slice_super / 2 + (i * slice_super)))
                        * 0.875,
                        np.cos(np.deg2rad(slice_super / 2 + (i * slice_super)))
                        * 0.875,
                    ),
                    (
                        np.sin(np.deg2rad(slice_super / 2 + (i * slice_super)))
                        * 3.125,
                        np.cos(np.deg2rad(slice_super / 2 + (i * slice_super)))
                        * 3.125,
                    ),
                ]
                for i in range(partition_super)
            ],
            colors="white",
        )
    )

    # Calculate edge weights and find edge color and style
    mid_style = (0, (4, 3))
    pos_style = (0, ())
    neg_style = (0, (1, 2))
    weights = nx.get_edge_attributes(G, "weight")
    weights.pop(("MDR", "AG", 0))
    normal_edges = list(weights.keys())
    weights = np.array(list(weights.values()))
    if graph_type == "frequency":
        edge_color_f = np.vectorize(
            lambda x: (
                "#9B0009"
                if x < 0.01  # noqa: PLR2004
                else "#FF9896"
                if x < 0.1  # noqa: PLR2004
                else "#0D4A70"
                if x >= 0.9  # noqa: PLR2004
                else "#52C4FF"
                if x >= 0.5  # noqa: PLR2004
                else "black"
            )
        )
        style_f = lambda x: (  # noqa: E731
            neg_style
            if x < 0.1  # noqa: PLR2004
            else pos_style
            if x >= 0.5  # noqa: PLR2004
            else mid_style
        )
    else:
        edge_color_f = np.vectorize(
            lambda x: (
                "black"
                if np.abs(x) < 0.05  # noqa: PLR2004
                else "#0D4A70"
                if x >= 0.9  # noqa: PLR2004
                else "#52C4FF"
                if x >= 0.05  # noqa: PLR2004
                else "#FF9896"
                if x > -0.9  # noqa: PLR2004
                else "#9B0009"
            )
        )
        style_f = lambda x: (  # noqa: E731
            mid_style
            if np.abs(x) < 0.05  # noqa: PLR2004
            else pos_style
            if x >= 0.05  # noqa: PLR2004
            else neg_style
        )
    edge_color = edge_color_f(weights)
    edge_style = [style_f(x) for x in weights]

    # Label partition's class (thanking the gods at stackoverflow for this:
    # https://stackoverflow.com/questions/19353576/curved-text-rendering-in-matplotlib)
    inner_partition_labels = [
        "beta-lactams",
        "peptides",
        "antiseptics",
        "quinolones",
        "sulfonamides",
        "amino- ",
        "diamino-",
        "bicyclomycins",
        "free fatty acids",
        "tetracyclines",
        "nucleosides",
        "nitroimidazoles",
        "tetra-",
        "phenicols",
        "oxazolidinones",
        "glycopeptides",
        "pleuromutilins",
        "MLS drugs"]

    inner_curves = [
        [
            np.sin(
                np.linspace(
                    np.deg2rad((partition - 0.5) * slice_super),
                    np.deg2rad((partition + 0.5) * slice_super),
                    100
                )
            )
            * 3.125,
            np.cos(
                np.linspace(
                    np.deg2rad((partition - 0.5) * slice_super),
                    np.deg2rad((partition + 0.5) * slice_super),
                    100
                )
            )
            * 3.125,
        ]
        for partition in range(partition_super)
    ]

    for curve, partition_label in zip(inner_curves, inner_partition_labels):
        # adding the text
        text = CurvedText(
            x=curve[0],
            y=curve[1],
            text=partition_label,
            va="bottom",
            fontsize=24,
            ax=ax
        )

    outer_partition_labels = [
        " ",
        " ",
        " ",
        "fluoro-",
        " ",
        "glycosides",
        "pyrimidines",
        " ",
        " ",
        " ",
        " ",
        " ",
        "cenomycins",
        " ",
        " ",
        " ",
        " ",
        " "]

    outer_curves = [
        [
            np.sin(
                np.linspace(
                    np.deg2rad((partition - 0.5) * slice_super),
                    np.deg2rad((partition + 0.5) * slice_super),
                    100
                )
            )
            * 3.25,
            np.cos(
                np.linspace(
                    np.deg2rad((partition - 0.5) * slice_super),
                    np.deg2rad((partition + 0.5) * slice_super),
                    100
                )
            )
            * 3.25,
        ]
        for partition in range(partition_super)
    ]

    for curve, partition_label in zip(outer_curves, outer_partition_labels):
        # adding the text
        text = CurvedText(
            x=curve[0],
            y=curve[1],
            text=partition_label,
            va="bottom",
            fontsize=24,
            ax=ax
        )

    # Draw graph, starting with nodes
    nx.draw_networkx_nodes(
        G=G,
        ax=ax,
        node_color=[
            "white" if node in sharing_classes else "black" for node in list(G)
        ],
        node_size=4000,
        edgecolors="grey",
        pos=pos,
    )

    # Add labels
    nx.draw_networkx_labels(
        G=G,
        ax=ax,
        font_color={
            node: "black" if node in sharing_classes else "white"
            for node in list(G)
        },
        font_size=25,
        pos=pos,
    )

    # Then add all edges except (MDR, AG)
    fancy_arrow_patches = nx.draw_networkx_edges(
        G=G,
        edgelist=normal_edges,
        arrows=True,
        arrowsize=35,
        arrowstyle="-|>",
        ax=ax,
        connectionstyle="arc3,rad=0.2",
        edge_color=edge_color,
        style=edge_style,
        node_size=4000,
        pos=pos,
        width=4,
    )

    # Change cap style of arrow edge
    for patch in fancy_arrow_patches:
        patch.set(capstyle="butt")

    # Make special edge for (MDR, AG)
    special_weight = np.array(
        [nx.get_edge_attributes(G, "weight")[("MDR", "AG", 0)]]
    )
    edge_color = edge_color_f(special_weight)
    edge_style = style_f(special_weight)
    fancy_arrow_patches = nx.draw_networkx_edges(
        G=G,
        edgelist=[("MDR", "AG", 0)],
        arrows=True,
        arrowsize=35,
        arrowstyle="-|>",
        ax=ax,
        connectionstyle="arc3,rad=0.5",
        edge_color=edge_color,
        style=edge_style,
        node_size=4000,
        pos=pos,
        width=4,
    )

    # Change cap style of arrow edge
    for patch in fancy_arrow_patches:
        patch.set(capstyle="butt")

    # Add node legend
    node_legend = fig.legend(
        handles=[
            Line2D(
                [0], [0], marker="o", c="white", ms=20, mfc="white", mec="grey"
            ),
            Line2D([0], [0], marker="o", c="white", ms=20, mfc="black"),
        ],
        labels=[
            "Shares a label with another class in database",
            "Doesn't share a label with another class in database",
        ],
        loc="upper right",
        fontsize=24
    )
    fig.add_artist(node_legend)

    # Add edge legend then save
    legend_artists = [
        Line2D([0], [0], c="#9B0009", ls=neg_style, lw=4),
        Line2D([0], [0], c="#FF9896", ls=neg_style, lw=4),
        Line2D([0], [0], c="#0D4A70", ls=pos_style, lw=4),
        Line2D([0], [0], c="#52C4FF", ls=pos_style, lw=4),
        Line2D([0], [0], c="black", ls=mid_style, lw=4),
    ]
    if graph_type == "share":
        fig.legend(
            handles=legend_artists[2:],
            labels=[
                "At least 90% of queries",
                "Less than 90% of queries",
                "Less than 5% of queries",
            ],
            ncols=3,
            loc="lower left",
            fontsize=24
        )

        fig.savefig(
            f"most_freq_share_label_with_deeparg_graph_{model}_{ident}.svg"
        )
    elif graph_type == "imbalance":
        fig.legend(
            handles=legend_artists,
            labels=[
                "Heavy skew toward DIAMOND hit",
                "Skew toward DIAMOND hit",
                "Heavy skew toward DeepARG prediction",
                "Skew toward DeepARG prediction",
                "No skew",
            ],
            ncols=3,
            loc="lower left",
            fontsize=23
        )

        fig.savefig(f"most_freq_largest_super_graph_{model}_{ident}.svg")
    elif graph_type == "frequency":
        fig.legend(
            handles=legend_artists,
            labels=[
                "Switches < 1% of pairs",
                "Switches < 10% of pairs",
                "Switches ≥ 90% of pairs",
                "Switches ≥ 50% of pairs",
                "Switches ≥ 10% and < 50% of pairs",
            ],
            ncols=3,
            loc="lower left",
            fontsize=23
        )

        fig.savefig(f"switch_frequecy_graph_{model}_{ident}.svg")
