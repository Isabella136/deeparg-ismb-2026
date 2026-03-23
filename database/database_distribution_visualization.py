import sys
import warnings
from itertools import product

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Circle

sys.path.append("../")

from utils.curved_text import CurvedText

VERSION = sys.argv[1]
FEATURE_DATA = f"v{VERSION}_feature_data.csv"
HIT_CLASSES = "../WGS_experiments/samples/hit_classes.tsv"


def get_shared_label_counts(
    multi_class_label_rows: pd.DataFrame,
    label_type: str,
    classes: set,
) -> pd.DataFrame:
    """Get counts of labels shared between two classes."""
    # Get list of labels shared by two or more classes
    labels_series = multi_class_label_rows[label_type]
    labels_list_dedup = labels_series.drop_duplicates().to_list()

    # Construct dataframe containg counts of labels shared between two classes
    shared_label_counts = pd.DataFrame(
        data=np.full(shape=len(classes) ** 2, fill_value=0),
        index=pd.MultiIndex.from_product([classes, classes]),
        columns=["count"],
    )
    for label in labels_list_dedup:
        # List of classes containing current label
        class_list = multi_class_label_rows[labels_series == label][
            "amr class"
        ].to_list()

        # Make index coords from x-product of classes containing current label
        coords = list(product(class_list, class_list))
        for x, y in coords:
            # Only add count if both classes are relevant
            if len({x, y}.intersection(classes)) == 2:  # noqa: PLR2004
                shared_label_counts.at[(x, y), "count"] += 1

    # Return dataframe
    return shared_label_counts


# Retrieve feature data, and replace amr class with abbreviation
amr_abbrev = pd.read_csv("amr_abbrev.csv", header=None, index_col=0)
feature_data = pd.read_csv(FEATURE_DATA, header=0, index_col=0)
feature_data["amr class"] = feature_data.apply(
    lambda x: amr_abbrev.at[x["amr class"], 1], axis=1)

# Get class count for DeepARG database from feature data
class_count_df = (
    feature_data
    .reset_index(drop=True)[["amr class", f"amr class v{VERSION} count"]]
    .drop_duplicates()
    .reset_index(drop=True)
)

# Get clstr|class count for DeepARG database from feature data
clstr_class_count_df = pd.DataFrame(
    feature_data.reset_index()[
        ["amr class", f"v{VERSION}-only cluster index"]
    ].value_counts(dropna=False, sort=False)
).reset_index(names=["amr class", "clstr"])

# Get superfamily|class count for DeepARG database from feature data
# (Counting combos in multi-superfamily features as their own domain)
combo_super_class_count_df = pd.DataFrame(
    feature_data.reset_index()[
        ["amr class", "superfamily(ies) id(s)"]
    ].value_counts(dropna=True, sort=False)
).reset_index(names=["amr class", "superfamily"])

# Create figure 1
plt.figure(figsize=(40, 20))
label_ax = plt.axes((0.05, 0.1, 0.92, 0.42))
feature_ax = plt.axes((0.05, 0.55, 0.92, 0.42), sharex=label_ax)
cb_palette = sns.color_palette("colorblind")

# Find the size of the biggest cluster per class
clstr_max = pd.DataFrame(
    data={"max": clstr_class_count_df.groupby("amr class")["count"].max()},
    index=clstr_class_count_df["amr class"].drop_duplicates().sort_values(),
)
clstr_max = clstr_max.reset_index(names="amr class")

# Find the size of the biggest superfamily per class
combo_super_max = pd.DataFrame(
    data={
        "max": combo_super_class_count_df.groupby(by="amr class")["count"].max()
    },
    index=combo_super_class_count_df["amr class"]
    .drop_duplicates()
    .sort_values(),
)
combo_super_max = combo_super_max.reset_index(names="amr class")

# Create vectorize function that sorts classes by gene sequence count
sorted_class = class_count_df.sort_values(
    by=f"amr class v{VERSION} count", ascending=False, ignore_index=True
)["amr class"]
sort_class_by_size = np.vectorize(
    lambda x: sorted_class[sorted_class == x].index[0]
)

# We need to combine clstr_max and combo_super_max s.t. we know which label is
# a cluster and which is superfamily. Additionally, we need to sort classes by
# gene sequence count
clstr_max.insert(loc=0, column="label", value="clstr|amr")
combo_super_max.insert(loc=0, column="label", value="super|amr")
label_max = pd.DataFrame(
    data=pd.concat([clstr_max, combo_super_max], join="inner"),
    columns=["label", "amr class", "max"],
).sort_values(by="amr class", key=sort_class_by_size)

# We also need to combine clstr_class_count and combo_superfamily_class_count
# s.t. we know which label is a cluster and which is superfamily. Additionally,
# we need to sort classes by gene sequence count
clstr_class_count_df.insert(loc=0, column="label", value="clstr|amr")
combo_super_class_count_df.insert(loc=0, column="label", value="super|amr")
label_count_df = pd.DataFrame(
    data=pd.concat(
        [clstr_class_count_df, combo_super_class_count_df], join="inner"
    ),
    columns=["label", "amr class", "count"],
).sort_values(by="amr class", key=sort_class_by_size)

# We actually need to plot label counts as ratios of the total class size
label_max["max ratio"] = label_max.apply(
    lambda x: (
        float(x["max"])
        / float(
            class_count_df.loc[class_count_df["amr class"] == x["amr class"]][  # noqa: PD009
                f"amr class v{VERSION} count"
            ].iat[0]
        )
    ),
    axis=1,
)
label_count_df["ratio"] = label_count_df.apply(
    lambda x: (
        float(x["count"])
        / float(
            class_count_df.loc[class_count_df["amr class"] == x["amr class"]][  # noqa: PD009
                f"amr class v{VERSION} count"
            ].iat[0]
        )
    ),
    axis=1,
)

# Make bar plot of median ratio and strip plot of max ratio
sns.barplot(
    data=label_count_df,
    x="amr class",
    y="ratio",
    hue="label",
    hue_order=["clstr|amr", "super|amr"],
    estimator="median",
    ax=label_ax,
    palette=cb_palette[:2],
    errorbar=None,
    alpha=0.5,
    legend=True,
)
sns.stripplot(
    data=label_max,
    x="amr class",
    y="max ratio",
    hue="label",
    hue_order=["clstr|amr", "super|amr"],
    ax=label_ax,
    dodge=True,
    size=15,
    palette=cb_palette[:2],
    legend=True,
)

# x-axis is list of classes sorted by size
label_ax.set_xlabel(xlabel="AMR Class", fontsize=35, loc="center", labelpad=0)
label_ax.set_xticks(
    ticks=label_ax.get_xticks(),
    labels=label_ax.get_xticklabels(),
    fontsize=30,
    rotation_mode="anchor",
    rotation=45,
    ha="right",
    va="center",
)

# y-axis is sequence count ratios, with upper limit set to 1.1 so that y = 1
# doesn't intersect with top of plot. Axis title is set later in the code
label_ax.set_ylabel(ylabel="")
label_ax.set_yticks(
    ticks=label_ax.get_yticks(),
    labels=label_ax.get_yticklabels(),
    fontsize=30,
    va="center",
)
label_ax.set_ylim(bottom=0, top=1.1)

# Set title and legend
label_ax.set_title("B", loc="left", fontsize=35, weight="bold")
label_ax.legend(
    handles=label_ax.legend_.legend_handles,
    labels=[
        "cluster median",
        "superfamily median",
        "cluster max",
        "superfamily max",
    ],
    fontsize=30,
    loc="upper left",
)

# Make histogram of sequence counts, and label actual count on top of bar
hist = sns.histplot(
    data=feature_data,
    x="amr class",
    ax=feature_ax,
    legend=True,
    fill=True,
    linewidth=3,
    element="bars",
    common_norm=False,
    shrink=0.8,
    color="black",
)
hist.bar_label(
    hist.containers[0],
    labels=hist.containers[0].datavalues,
    padding=2,
    fontsize=27,
)

# y-axis is sequence count. Axis title is set later in the code
feature_ax.set_ylabel(ylabel="", fontsize=35)
feature_ax.set_yticks(
    ticks=feature_ax.get_yticks(),
    labels=feature_ax.get_yticklabels(),
    fontsize=30,
)

# Since top feature_ax plot is vertically aligned to bottom label_ax plot, we
# don't need to label feature_ax x-axis, nor do we need tick marks
feature_ax.set_xlabel(xlabel="")
feature_ax.tick_params(axis="x", labelbottom=False)

# Set title
feature_ax.set_title("A", loc="left", fontsize=35, weight="bold")

# We don't want a gap between y-axis line and MDR bar
label_ax.margins(x=0)

# We want y-axis labels for both graphs to be vertically aligned
plt.gcf().text(
    0.008, 0.31, "Sequence Count Ratio", rotation=90, fontsize=35, va="center"
)
plt.gcf().text(
    0.008, 0.77, "Sequence Count", rotation=90, fontsize=35, va="center"
)

# Let's save and be done with it
plt.savefig(f"db_distr_v{VERSION}.png")

# Create figure 2
plt.figure(figsize=(15, 30))
super_ax = plt.axes((0.01, 0.0, 0.98, 0.49))
clstr_ax = plt.axes((0.01, 0.49, 0.98, 0.49))

amr_abbrev.at["polyamine:peptide", 1] = "PA/\nPEP"
clstr_class_count_df["amr class"] = clstr_class_count_df["amr class"].apply(
    lambda x: x if x != "PA/PEP" else "PA/\nPEP"
)
combo_super_class_count_df["amr class"] = combo_super_class_count_df[
    "amr class"
].apply(lambda x: x if x != "PA/PEP" else "PA/\nPEP")

# We will only include classes that will have alignment hits in WGS experiments
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    interesting_nodes = pd.read_csv(HIT_CLASSES).apply(
        lambda x: amr_abbrev.at[x[0], 1], axis=1
    )

# Find clusters shared across two or more classes
class_per_clstr_count_df = pd.DataFrame(
    clstr_class_count_df["clstr"].value_counts(dropna=True, sort=False)
).reset_index(names="clstr")
multi_class_clstr = class_per_clstr_count_df.loc[
    class_per_clstr_count_df["count"] > 1
]["clstr"].to_list()
multi_class_clstr_rows = clstr_class_count_df.loc[
    clstr_class_count_df["clstr"].isin(multi_class_clstr)
]

# Get counts of clusters shared between two classes
clstr_share_matrix = get_shared_label_counts(
    multi_class_clstr_rows, "clstr", interesting_nodes
)

# Create node graph from set of classes that actually do share clusters
G_clstr = nx.DiGraph()
G_clstr_nodes = {
    row[0][0]
    for row in clstr_share_matrix.iterrows()
    if (row[1].iat[0] > 0)  # noqa: PD009
}
G_clstr.add_nodes_from(G_clstr_nodes)

# Specify a specific direction for edges. Only matters becaus edges have a
# counterclockwise curve, and certain edge orientations heavily overlap with
# nodes
G_clstr_edges = [
    ("FFA", "FQ"),
    ("MDR", "PHE"),
    ("MLS", "PHE"),
    ("MDR", "TRI"),
    ("PMX", "PEP"),
    ("MDR", "TET"),
    ("MLS", "TET"),
    ("TET", "FFA"),
    ("UNC", "MDR"),
    ("UNC", "MLS"),
    ("AG", "UNC"),
    ("UNC", "GlyP"),
    ("UNC", "PHE"),
    ("MDR", "BL"),
    ("MDR", "AG"),
    ("MDR", "AC"),
    ("MDR", "FQ"),
    ("MDR", "FFA"),
    ("MDR", "MLS"),
]

# Add edges and weights (shared cluster count) to node graph
G_clstr.add_weighted_edges_from([
    (row[0], row[1], int(clstr_share_matrix.at[row, "count"]))
    for row in G_clstr_edges
])

# Specify a specific position for each node. Based on class specificity.
partition_cluster = 11
slice_cluster = 360.0 / partition_cluster
pos = {
    "AC": (
        np.sin(np.deg2rad(3 * slice_cluster)) * 2.00,
        np.cos(np.deg2rad(3 * slice_cluster)) * 2.00,
    ),
    "AG": (
        np.sin(np.deg2rad(4 * slice_cluster)) * 2.00,
        np.cos(np.deg2rad(4 * slice_cluster)) * 2.00,
    ),
    "BL": (
        np.sin(np.deg2rad(0 * slice_cluster)) * 1.25,
        np.cos(np.deg2rad(0 * slice_cluster)) * 1.25,
    ),
    "FFA": (
        np.sin(np.deg2rad(6 * slice_cluster)) * 2.00,
        np.cos(np.deg2rad(6 * slice_cluster)) * 2.00,
    ),
    "FQ": (
        np.sin(np.deg2rad(5 * slice_cluster)) * 2.00,
        np.cos(np.deg2rad(5 * slice_cluster)) * 2.00,
    ),
    "GlyP": (
        np.sin(np.deg2rad(9 * slice_cluster)) * 2.00,
        np.cos(np.deg2rad(9 * slice_cluster)) * 2.00,
    ),
    "MDR": (
        np.sin(np.deg2rad(135)) * 0.50,
        np.cos(np.deg2rad(135)) * 0.50
    ),
    "MLS": (
        np.sin(np.deg2rad(10 * slice_cluster)) * 1.25,
        np.cos(np.deg2rad(10 * slice_cluster)) * 1.25,
    ),
    "PEP": (
        np.sin(np.deg2rad(2 * slice_cluster)) * 2.00,
        np.cos(np.deg2rad(2 * slice_cluster)) * 2.00,
    ),
    "PHE": (
        np.sin(np.deg2rad(8 * slice_cluster)) * 2.00,
        np.cos(np.deg2rad(8 * slice_cluster)) * 2.00,
    ),
    "PMX": (
        np.sin(np.deg2rad(2 * slice_cluster)) * 2.75,
        np.cos(np.deg2rad(2 * slice_cluster)) * 2.75,
    ),
    "TET": (
        np.sin(np.deg2rad(7 * slice_cluster)) * 2.00,
        np.cos(np.deg2rad(7 * slice_cluster)) * 2.00,
    ),
    "TRI": (
        np.sin(np.deg2rad(1 * slice_cluster)) * 2.75,
        np.cos(np.deg2rad(1 * slice_cluster)) * 2.75,
    ),
    "UNC": (
        np.sin(np.deg2rad(-45)) * 0.50,
        np.cos(np.deg2rad(-45)) * 0.50
    ),
}

# Specify class specificity with circles
clstr_ax.add_collection(
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
clstr_ax.add_collection(
    LineCollection(
        [
            [
                (
                    np.sin(np.deg2rad(slice_cluster / 2 + (i * slice_cluster)))
                    * 0.875,
                    np.cos(np.deg2rad(slice_cluster / 2 + (i * slice_cluster)))
                    * 0.875,
                ),
                (
                    np.sin(np.deg2rad(slice_cluster / 2 + (i * slice_cluster)))
                    * 3.125,
                    np.cos(np.deg2rad(slice_cluster / 2 + (i * slice_cluster)))
                    * 3.125,
                ),
            ]
            for i in range(partition_cluster)
        ],
        colors="white",
    )
)

# The tough part: label partition's class (thanking the gods at stackoverflow
# for this: https://stackoverflow.com/questions/19353576/curved-text-rendering-in-matplotlib)
partition_labels = [
    "beta-lactams",
    "antiseptics",
    "peptides",
    "aminocoumarins",
    "aminoglycosides",
    "fluoroquinolones",
    "free fatty acids",
    "tetracyclines",
    "phenicols",
    "glycopeptides",
    "MLS drugs"]

curves = [
    [
        np.sin(
            np.linspace(
                np.deg2rad((partition - 0.5) * slice_cluster),
                np.deg2rad((partition + 0.5) * slice_cluster),
                100
            )
        )
        * 3.125,
        np.cos(
            np.linspace(
                np.deg2rad((partition - 0.5) * slice_cluster),
                np.deg2rad((partition + 0.5) * slice_cluster),
                100
            )
        )
        * 3.125,
    ]
    for partition in range(partition_cluster)
]

for curve, partition_label in zip(curves, partition_labels):
    # adding the text
    text = CurvedText(
        x=curve[0],
        y=curve[1],
        text=partition_label,
        va="bottom",
        fontsize=25,
        ax=clstr_ax
    )

# Draw the actual nodes and edges
nx.draw(
    G_clstr,
    pos=pos,
    with_labels=True,
    ax=clstr_ax,
    width=list(nx.get_edge_attributes(G_clstr, "weight").values()),
    arrows=True,
    arrowstyle="-",
    connectionstyle="arc3,rad=0.3",
    node_color="white",
    font_color="black",
    font_size=25,
    node_size=5000,
)

# Set title
clstr_ax.set_title("A", loc="left", va="top", fontsize=35, weight="bold")

# Find superfamilies shared across two or more classes
class_per_super_count_df = pd.DataFrame(
    combo_super_class_count_df["superfamily"].value_counts(
        dropna=True, sort=False
    )
).reset_index(names="superfamily")
multi_class_supers = class_per_super_count_df.loc[
    class_per_super_count_df["count"] > 1
]["superfamily"].to_list()
multi_class_super_rows = combo_super_class_count_df.loc[
    combo_super_class_count_df["superfamily"].isin(multi_class_supers)
]

# Get counts of superfamilies shared between two classes
super_share_matrix = get_shared_label_counts(
    multi_class_super_rows, "superfamily", interesting_nodes
)

# Create node graph from set of classes that actually do share superfamilies
G_super = nx.DiGraph()
G_super_nodes = {
    row[0][0]
    for row in super_share_matrix.iterrows()
    if (row[1].iat[0] > 0)  # noqa: PD009
}
G_super.add_nodes_from(G_super_nodes)

# Specify a specific direction for edges. Only matters becaus edges have a
# counterclockwise curve, and certain edge orientations heavily overlap with
# nodes
G_super_edges = [
    ("FQ", "MDR"),
    ("FQ", "MLS"),
    ("FQ", "TET"),
    ("FQ", "UNC"),
    ("FQ", "AG"),
    ("FQ", "FFA"),
    ("BCM", "FQ"),
    ("PHE", "FQ"),
    ("TET-C", "FQ"),
    ("NUC", "FQ"),
    ("AG", "MDR"),
    ("AC", "MDR"),
    ("MLS", "AG"),
    ("MLS", "AC"),
    ("AG", "UNC"),
    ("PHE", "AG"),
    ("AG", "PA/\nPEP"),
    ("NUC", "AG"),
    ("MLS", "MDR"),
    ("BCM", "MLS"),
    ("MLS", "TET"),
    ("MLS", "UNC"),
    ("MLS", "PHE"),
    ("MLS", "TET-C"),
    ("MLS", "NUC"),
    ("TRI", "MLS"),
    ("MLS", "OXA"),
    ("MDR", "TRI"),
    ("OXA", "MDR"),
    ("UNC", "MDR"),
    ("GlyP", "UNC"),
    ("UNC", "PHE"),
    ("PEP", "PMX"),
    ("BAC", "PEP"),
    ("BAC", "UNC"),
    ("BAC", "MLS"),
    ("FOF", "GlyP"),
    ("MLS", "FOM"),
    ("FOM", "TET-C"),
    ("FOM", "PHE"),
    ("FOM", "NUC"),
    ("MDR", "FOM"),
    ("FOM", "TET"),
    ("FOM", "BCM"),
    ("FQ", "FOM"),
    ("TET-C", "MDR"),
    ("TET", "TET-C"),
    ("TET-C", "PHE"),
    ("NUC", "TET-C"),
    ("TET-C", "BCM"),
    ("MDR", "BL"),
    ("PEP", "BL"),
    ("TET", "PHE"),
    ("PHE", "MDR"),
    ("PHE", "NUC"),
    ("BCM", "PHE"),
    ("MDR", "PEP"),
    ("MDR", "FFA"),
    ("MDR", "BCM"),
    ("MDR", "TET"),
    ("MDR", "NUC"),
    ("NUC", "TET"),
    ("BCM", "NUC"),
    ("TET", "FFA"),
    ("BCM", "TET"),
]

# Add edges and weights (shared superfamily count) to node graph
G_super.add_weighted_edges_from([
    (row[0], row[1], int(super_share_matrix.at[row, "count"]))
    for row in G_super_edges
])

# Specify a specific position for each node. Based on class specificity.
partition_super = 16
slice_super = 360.0 / partition_super
pos = {
    "AC": (
        np.sin(np.deg2rad(3 * slice_super)) * 2.00,
        np.cos(np.deg2rad(3 * slice_super)) * 2.00,
    ),
    "AG": (
        np.sin(np.deg2rad(4 * slice_super)) * 2.00,
        np.cos(np.deg2rad(4 * slice_super)) * 2.00,
    ),
    "BAC": (
        np.sin(np.deg2rad(1.75 * slice_super)) * 2.75,
        np.cos(np.deg2rad(1.75 * slice_super)) * 2.75,
    ),
    "BCM": (
        np.sin(np.deg2rad(6 * slice_super)) * 2.00,
        np.cos(np.deg2rad(6 * slice_super)) * 2.00,
    ),
    "BL": (
        np.sin(np.deg2rad(16 * slice_super)) * 1.25,
        np.cos(np.deg2rad(16 * slice_super)) * 1.25,
    ),
    "FFA": (
        np.sin(np.deg2rad(7 * slice_super)) * 2.00,
        np.cos(np.deg2rad(7 * slice_super)) * 2.00,
    ),
    "FOF": (
        np.sin(np.deg2rad(9.25 * slice_super)) * 2.75,
        np.cos(np.deg2rad(9.25 * slice_super)) * 2.75,
    ),
    "FOM": (
        np.sin(np.deg2rad(8.75 * slice_super)) * 2.75,
        np.cos(np.deg2rad(8.75 * slice_super)) * 2.75,
    ),
    "FQ": (
        np.sin(np.deg2rad(5 * slice_super)) * 2.00,
        np.cos(np.deg2rad(5 * slice_super)) * 2.00,
    ),
    "GlyP": (
        np.sin(np.deg2rad(14 * slice_super)) * 2.00,
        np.cos(np.deg2rad(14 * slice_super)) * 2.00,
    ),
    "MDR": (
        np.sin(np.deg2rad(160)) * 0.50,
        np.cos(np.deg2rad(160)) * 0.50
    ),
    "MLS": (
        np.sin(np.deg2rad(15 * slice_super)) * 1.25,
        np.cos(np.deg2rad(15 * slice_super)) * 1.25,
    ),
    "NUC": (
        np.sin(np.deg2rad(10 * slice_super)) * 2.00,
        np.cos(np.deg2rad(10 * slice_super)) * 2.00,
    ),
    "OXA": (
        np.sin(np.deg2rad(13 * slice_super)) * 2.00,
        np.cos(np.deg2rad(13 * slice_super)) * 2.00,
    ),
    "PA/\nPEP": (
        np.sin(np.deg2rad(2 * slice_super)) * 1.25,
        np.cos(np.deg2rad(2 * slice_super)) * 1.25,
    ),
    "PEP": (
        np.sin(np.deg2rad(2 * slice_super)) * 2.00,
        np.cos(np.deg2rad(2 * slice_super)) * 2.00,
    ),
    "PHE": (
        np.sin(np.deg2rad(11 * slice_super)) * 2.00,
        np.cos(np.deg2rad(11 * slice_super)) * 2.00,
    ),
    "PMX": (
        np.sin(np.deg2rad(2.25 * slice_super)) * 2.75,
        np.cos(np.deg2rad(2.25 * slice_super)) * 2.75,
    ),
    "TET": (
        np.sin(np.deg2rad(8 * slice_super)) * 2.00,
        np.cos(np.deg2rad(8 * slice_super)) * 2.00,
    ),
    "TET-C": (
        np.sin(np.deg2rad(12 * slice_super)) * 2.00,
        np.cos(np.deg2rad(12 * slice_super)) * 2.00,
    ),
    "TRI": (
        np.sin(np.deg2rad(1 * slice_super)) * 2.75,
        np.cos(np.deg2rad(1 * slice_super)) * 2.75,
    ),
    "UNC": (
        np.sin(np.deg2rad(-20)) * 0.50,
        np.cos(np.deg2rad(-20)) * 0.50
    ),
}

# Specify class specificity with circles
super_ax.add_collection(
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
super_ax.add_collection(
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

# The tough part: label partition's class (thanking the gods at stackoverflow
# for this: https://stackoverflow.com/questions/19353576/curved-text-rendering-in-matplotlib)
partition_labels = [
    "beta-lactams",
    "antiseptics",
    "peptides",
    "aminocoumarins",
    "aminoglycosides",
    "fluoroquinolones",
    "bicyclomycins",
    "free fatty acids",
    "tetracyclines",
    "phosphonates",
    "nucleoside",
    "phenicols",
    "tetracenomycin C",
    "oxazolidinones",
    "glycopeptides",
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
        ax=super_ax
    )

# Draw the actual nodes and edges
nx.draw(
    G_super,
    pos=pos,
    with_labels=True,
    ax=super_ax,
    width=list(nx.get_edge_attributes(G_super, "weight").values()),
    arrows=True,
    arrowstyle="-",
    connectionstyle="arc3,rad=0.3",
    node_color="white",
    font_color="black",
    font_size=25,
    node_size=5000,
)

# Set title
super_ax.set_title("B", loc="left", va="top", fontsize=35, weight="bold")

# Let's save and be done with it
plt.savefig(f"inter_class_sim_graph_v{VERSION}.png")
