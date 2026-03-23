import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Circle

sys.path.append("../../")

from utils.curved_text import CurvedText

HIT_COUNT = "../samples/deeparg_hit_count.tsv"
INPUT_SEQ_COUNT = "../samples/sequence_count.tsv"
FEATURE_DATA = "../../database/v2_feature_data.csv"

mpl.rcParams["mathtext.fontset"] = "dejavusans"
mpl.rcParams["font.family"] = "DejaVu Sans"

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
hit_count_amr_df["diff percentage"] = hit_count_amr_df.apply(
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
perc_diff_ax = plt.axes((0.15, 0.07, 0.8, 0.28))
perc_pred_ax = plt.axes((0.15, 0.385, 0.8, 0.28), sharex=perc_diff_ax)
perc_hit_ax = plt.axes((0.15, 0.7, 0.8, 0.28), sharex=perc_diff_ax)

# Draw figure 3A: percentages of sequences that had an alignment hit, averaged
# for each model–alignment identity run condition
sns.barplot(
    x=range(6),
    y=hit_count_amr_df.groupby(["model", "alignment identity"])[
        "alignment hit percentage"
    ].mean(),
    ax=perc_hit_ax,
)

# Since top perc_hit_ax plot is vertically aligned to bottom perc_diff_ax
# plot, we don't need to label perc_hit_ax x-axis, nor do we need tick marks
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
    x=range(6),
    y=hit_count_amr_df.groupby(["model", "alignment identity"])[
        "deeparg prediction percentage"
    ].mean(),
    ax=perc_pred_ax,
)

# Since middle perc_pred_ax plot is vertically aligned to bottom perc_diff_ax
# plot, we don't need to label perc_pred_ax x-axis, nor do we need tick marks
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
    x=range(6),
    y=hit_count_amr_df.groupby(["model", "alignment identity"])[
        "diff percentage"
    ].mean(),
    color=cb_palette[0],
    ax=perc_diff_ax,
)

# x-axis is the list of model–alignment identity run condition
perc_diff_ax.set_xlabel(
    xlabel="Model–Alignment Identity Condition",
    fontsize=35,
    loc="center",
    labelpad=0,
)
perc_diff_ax.set_xticks(
    ticks=range(6),
    labels=["LS–30%", "LS–50%", "LS–80%", "SS–30%", "SS–50%", "SS–80%"],
    rotation_mode="anchor",
    ha="right",
    va="top",
    rotation=45,
    fontsize=30,
)

# y-axis is average percentages from 0 to 5. Axis title is set later in the code
perc_diff_ax.set_ylabel("")
perc_diff_ax.set_ylim(top=0.06)
perc_diff_ax.set_yticks(
    ticks=[0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06],
    labels=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0], fontsize=30)

# Set title
plt.gcf().text(0.02, 0.35, "C", fontsize=40, va="center", weight="bold")

# Styling the graph a bit more
perc_diff_ax.tick_params(length=0, pad=8)
perc_diff_ax.grid(visible=True, which="major", axis="y", color="0.75")

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
plt.savefig("barplot.png")

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

# Get individual superfamily|class for DeepARG database from feature data
# (Consider combos in multi-superfamily features as their own domain)
combo_super_class_df = (
    feature_data
    .reset_index(drop=True)
    [["amr class", "superfamily(ies) id(s)"]]
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
    combo_super_class_df["superfamily(ies) id(s)"]
    .value_counts(dropna=True, sort=False)
    .reset_index()
    .set_axis(["super", "count"], axis=1)
)
multi_class_super = classes_per_super_count_df.loc[
    classes_per_super_count_df["count"] > 1
]["super"].to_list()
shared_super_class = combo_super_class_df.loc[
    combo_super_class_df["superfamily(ies) id(s)"].isin(multi_class_super)
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

for graph_type in ["alignment", "share", "imbalance"]:
    # Compare observed switches vs expected switches
    if graph_type == "alignment":
        # Calculate the number of hits to amr class for a given query
        amr_count_df = (
            label_df[["Sample ID", "Query", "amr", "count"]]
            .groupby(["Sample ID", "Query", "amr"])
            .sum()
        )

        exp_df = pd.DataFrame(
            data=np.full(
                shape=(heatmap_x_axis.shape[0], heatmap_y_axis.shape[0]),
                fill_value=np.nan,
            ),
            index=heatmap_x_axis,
            columns=heatmap_y_axis,
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

        # Calculate the number of hits to y for queries where y is not the most
        # frequent class hit
        y_alignments_counts = x_is_best_y_in_alignment[
            ["Sample ID", "Most Frequent Class", "Query", "amr"]
        ].drop_duplicates()
        y_alignments_counts["count"] = y_alignments_counts.apply(
            lambda x: amr_count_df.at[
                (x["Sample ID"], x["Query"], x["amr"]), "count"
            ],
            axis=1,
        )
        y_alignments_counts = y_alignments_counts.set_index(
            ["Sample ID", "Most Frequent Class", "Query", "amr"]
        )

        # Calculate the number of total hits for queries where y is not the most
        # frequent class hit
        all_alignment_counts = amr_count_df.groupby(
            level=["Sample ID", "Query"]
        )[["count"]].sum()
        all_alignment_counts = (
            y_alignments_counts
            .reset_index()
            .apply(
                lambda x: pd.Series(
                    data={
                        "Sample ID": x["Sample ID"],
                        "Most Frequent Class": x["Most Frequent Class"],
                        "Query": x["Query"],
                        "amr": x["amr"],
                        "count": all_alignment_counts.loc[
                            (x["Sample ID"], x["Query"])
                        ].at["count"],
                    }
                ),
                axis=1,
            )
            .set_index(["Sample ID", "Most Frequent Class", "Query", "amr"])
        )

        # Caluclate the probability of choosing y at random
        probabilities = y_alignments_counts.div(all_alignment_counts)
        probabilities = probabilities.groupby(
            by=["Most Frequent Class", "amr"]
        )["count"].sum()

        for row in pair_df.iterrows():
            # For x- and y-axis
            deeparg_class = row[0][1]
            diamond_class = row[0][0]

            # For switch_df
            switch = row[1]["Query"]

            # For exp_db
            expected = probabilities.at[(diamond_class, deeparg_class)]

            # Insert in dfs
            switch_df.at[deeparg_class, diamond_class] = float(switch)
            exp_df.at[deeparg_class, diamond_class] = float(expected)

        # Now let's get the heatmap values
        heatmap_df = (
            switch_df
            .add(1, fill_value=np.nan)
            .div(exp_df.add(1), fill_value=np.nan)
            .map(np.log, na_action="ignore")
            .sort_index()
        )

    # See if switches coincide with label sharing
    elif graph_type == "share":
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

    else:
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

    fig = plt.figure(figsize=[15, 15])
    ax = fig.add_axes((0.01, 0.01, 0.98, 0.98))

    # Create graph objects along with nodes and edges
    G = nx.MultiDiGraph()
    G.add_nodes_from(set(heatmap_x_axis.to_list()).union(set(heatmap_y_axis)))
    G.add_weighted_edges_from([
        (y, x, round(heatmap_df.at[x, y], 2))
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
        "TET-C": (
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
            facecolors=["slategray", "white", "white", "white"],
            alpha=[0.45, 0.2, 0.2, 0.2],
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
    weights = nx.get_edge_attributes(G, "weight")
    weights.pop(("MDR", "AG", 0))
    normal_edges = list(weights.keys())
    weights = np.array(list(weights.values()))
    if graph_type == "alignment":
        width = np.abs(weights) + 1
        edge_color_f = np.vectorize(
            lambda x: (
                "black"
                if np.abs(x) <= np.log(1.5)
                else "royalblue"
                if x >= np.log(10)
                else "cornflowerblue"
                if x > np.log(1.5)
                else "lightcoral"
                if x > np.log(0.1)
                else "red"
            )
        )
        style_f = np.vectorize(
            lambda x: (
                "-"
                if np.abs(x) <= np.log(1.5)
                else "--"
                if x > np.log(1.5)
                else "-."
            )
        )
    else:
        width = (np.abs(weights) + 1) * 2
        edge_color_f = np.vectorize(
            lambda x: (
                "black"
                if np.abs(x) <= 0.05  # noqa: PLR2004
                else "royalblue"
                if x >= 0.9  # noqa: PLR2004
                else "cornflowerblue"
                if x > 0.05  # noqa: PLR2004
                else "lightcoral"
                if x > -0.9  # noqa: PLR2004
                else "red"
            )
        )
        style_f = np.vectorize(
            lambda x: (
                "-"
                if np.abs(x) <= 0.05  # noqa: PLR2004
                else "--"
                if x > 0.05  # noqa: PLR2004
                else "-."
            )
        )
    edge_color = edge_color_f(weights)
    edge_style = style_f(weights)

    # The tough part: label partition's class (thanking the gods at stackoverflow
    # for this: https://stackoverflow.com/questions/19353576/curved-text-rendering-in-matplotlib)
    partition_labels = [
        "beta-lactams",
        "peptides",
        "antiseptics",
        "fluoroquinolones",
        "sulfonamides",
        "aminoglycosides",
        "diaminopyrimidines",
        "bicyclomycins",
        "free fatty acids",
        "tetracyclines",
        "nucleosides",
        "nitroimidazoles",
        "tetracenomycin C",
        "phenicols",
        "oxazolidinones",
        "glycopeptides",
        "pleuromutilins",
        "MLS drugs"]

    curves = [
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

    for curve, partition_label in zip(curves, partition_labels):
        # adding the text
        text = CurvedText(
            x=curve[0],
            y=curve[1],
            text=partition_label,
            va="bottom",
            fontsize=20,
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
    nx.draw_networkx_edges(
        G=G,
        edgelist=normal_edges,
        arrows=True,
        arrowsize=35,
        ax=ax,
        connectionstyle="arc3,rad=0.2",
        edge_color=edge_color,
        style=edge_style,
        node_size=4000,
        pos=pos,
        width=width,
    )

    # Make special edge for (MDR, AG)
    special_weight = np.array(
        [nx.get_edge_attributes(G, "weight")[("MDR", "AG", 0)]]
    )
    if graph_type == "alignment":
        width = np.abs(special_weight) + 1
    else:
        width = (np.abs(special_weight) + 1) * 2
    edge_color = edge_color_f(special_weight)
    edge_style = style_f(special_weight)
    nx.draw_networkx_edges(
        G=G,
        edgelist=[("MDR", "AG", 0)],
        arrows=True,
        arrowsize=35,
        ax=ax,
        connectionstyle="arc3,rad=0.5",
        edge_color=edge_color,
        style=edge_style,
        node_size=4000,
        pos=pos,
        width=width,
    )

    # Save
    if graph_type == "alignment":
        fig.savefig(
            f"most_freq_amr_switch_relative_to_alignment_graph_{model}_{ident}.png"
        )
    elif graph_type == "share":
        fig.savefig(
            f"most_freq_share_label_with_deeparg_graph_{model}_{ident}.png"
        )
    elif graph_type == "imbalance":
        fig.savefig(f"most_freq_largest_super_graph_{model}_{ident}.png")

    continue
    ax_left = fig.add_axes((0.12, 0.12, 0.75, 0.84))
    cbar = fig.add_axes((0.9, 0.24, 0.02, 0.6))

    heatmap_df.transpose().to_csv(f"{graph_type}.csv")

    sns.heatmap(
        data=heatmap_df.transpose(),
        # mask=mask_df.transpose(),
        # cmap=custom_cmap,
        center=0,
        vmin=0
        if graph_type == "share"
        else -1
        if graph_type == "imbalance"
        else -6,
        vmax=1 if graph_type != "alignment" else 3,
        ax=ax_left,
        cbar_ax=cbar,
    )

    cbar.tick_params(labelsize=30)

    # Label x-axis
    ax_left.set_xlabel("DeepARG prediction (Y)", fontsize=35)

    # Label y-axis
    ax_left.set_ylabel(ylabel="alignment-based prediction (X)", fontsize=35)

    # Label x-ticks
    ax_left.set_xticks(
        ticks=ax_left.get_xticks(),
        labels=heatmap_df.index.drop_duplicates(),
        rotation_mode="anchor",
        rotation=45,
        ha="right",
        va="top",
        fontsize=30,
    )

    # Label y-ticks
    ax_left.set_yticks(
        ticks=ax_left.get_yticks(),
        labels=heatmap_df.columns.drop_duplicates(),
        rotation_mode="anchor",
        rotation=0,
        ha="right",
        va="center",
        fontsize=30,
    )

    y_boundaries = range(heatmap_df.columns.size)
    x_boundaries = range(heatmap_df.index.size)

    for y_loc in y_boundaries[1:]:
        ax_left.hlines(
            y_loc,
            xmin=ax_left.get_xlim()[0],
            xmax=ax_left.get_xlim()[1],
            colors="white",
            linewidth=5,
        )
    for x_loc in x_boundaries[1:]:
        ax_left.vlines(
            x_loc,
            ymin=ax_left.get_ylim()[0],
            ymax=ax_left.get_ylim()[1],
            colors="white",
            linewidth=5,
        )

    if graph_type == "alignment":
        fig.savefig(
            f"most_freq_amr_switch_relative_to_alignment_{model}_{ident}.png"
        )
    elif graph_type == "share":
        fig.savefig(f"most_freq_share_label_with_deeparg_{model}_{ident}.png")
    elif graph_type == "imbalance":
        fig.savefig(f"most_freq_largest_super_{model}_{ident}.png")
