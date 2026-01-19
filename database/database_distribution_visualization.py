import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from itertools import product

import networkx as nx
import seaborn as sn
import pandas as pd
import numpy as np
import sys

def make_label_share_matrix(df: pd.DataFrame, classes: pd.Series, label: str) -> pd.DataFrame:
    labels = df[label].drop_duplicates().to_list()
    

    label_share = pd.DataFrame(
        data=np.full(shape=len(classes)**2,fill_value=0),
        index=pd.MultiIndex.from_product([classes, classes]),
        columns=["count"])
    
    for l in labels:
        class_list = df[df[label] == l]["amr class"].to_list()
        coords = list(product(class_list, class_list))
        for x, y in coords:
            label_share.at[(x, y), "count"] += 1
    
    return label_share

def make_label_adj_matrix(df: pd.DataFrame, classes: pd.Series, label: str) -> pd.DataFrame:
    labels = df[label].drop_duplicates().to_list()
    

    label_share = pd.DataFrame(
        data=np.full(shape=(len(classes), len(classes)),fill_value=0),
        index=classes,
        columns=classes)
    
    for l in labels:
        class_list = df[df[label] == l]["amr class"].to_list()
        coords = list(product(class_list, class_list))
        for x, y in coords:
            label_share.at[x, y] += 1
    
    return label_share

VERSION = sys.argv[1]
FEATURE_DATA = f"v{VERSION}_feature_data.csv"

abbrev_file = open("amr_abbrev.csv", "r")
amr_abbrev = dict()
for line in abbrev_file.readlines():
    row = line.strip().split(",")
    amr_abbrev.update({row[0]: row[1]})
abbrev_file.close()

feature = pd.read_csv(FEATURE_DATA, header=0, index_col=0)
feature["amr class"] = feature.apply(
    lambda x: amr_abbrev[x["amr class"]], axis=1)

# Get class count for DeepARG database from feature data
class_count_df = (feature
    .reset_index(drop=True)[["amr class", f"amr class v{VERSION} count"]]
    .drop_duplicates()
    .reset_index(drop=True))

# Get clstr|class count for DeepARG database from feature data
clstr_class_count_df = (pd.DataFrame(feature
        .reset_index()[["amr class", f"v{VERSION}-only cluster index"]]
        .value_counts(dropna=False, sort=False))
    .reset_index(names=["amr class", "clstr"]))

# Get superfamily|class count for DeepARG database from feature data
# (Counting combos in multi-superfamily features as their own domain)
combo_superfamily_class_count_df = (pd.DataFrame(feature
        .reset_index()[["amr class", "superfamily(ies) id(s)"]]
        .value_counts(dropna=True, sort=False))
    .reset_index(names=["amr class", "superfamily"]))

clstr_stats = pd.DataFrame(
    data={
        "max": clstr_class_count_df.groupby("amr class")["count"].max()},
    index=clstr_class_count_df["amr class"].drop_duplicates().sort_values())
clstr_stats.reset_index(names="amr class", inplace=True)
clstr_stats.insert(
    loc=0,
    column="label",
    value="clstr|amr")
combo_super_stats = pd.DataFrame(
    data={
        "max": combo_superfamily_class_count_df.groupby(by="amr class")["count"].max()},
    index=combo_superfamily_class_count_df["amr class"].drop_duplicates().sort_values())
combo_super_stats.reset_index(names="amr class", inplace=True)
combo_super_stats.insert(
    loc=0,
    column="label",
    value="super|amr")
sorted_class = (class_count_df
    .sort_values(
        by=f"amr class v{VERSION} count", 
        ascending=False, 
        ignore_index=True)["amr class"])
sort_class_by_size = np.vectorize(lambda x: sorted_class[sorted_class==x].index[0])
label_stats = pd.DataFrame(
    data= pd.concat([
        clstr_stats,
        combo_super_stats], join='inner'),
    columns=["label", "amr class", "max"]).sort_values(
        by="amr class",
        key=sort_class_by_size)
label_stats["max ratio"] = label_stats.apply(
    lambda x: float(x["max"])/float(
        class_count_df.loc[class_count_df["amr class"]==x["amr class"]][f"amr class v{VERSION} count"].iat[0]),
    axis=1)

# We need to combine clstr_class_count and combo_superfamily_class_count
# s.t. we know which label is clstr and which is super:
clstr_class_count_df.insert(
    loc=0,
    column="label",
    value="clstr|amr")
combo_superfamily_class_count_df.insert(
    loc=0,
    column="label",
    value="super|amr")
label_distr_df = pd.DataFrame(
    data= pd.concat([
        clstr_class_count_df,
        combo_superfamily_class_count_df], join='inner'),
    columns=["label", "amr class", "count"]).sort_values(
        by="amr class",
        key=sort_class_by_size)
label_distr_df["ratio"] = label_distr_df.apply(
    lambda x: float(x["count"])/float(
        class_count_df.loc[class_count_df["amr class"]==x["amr class"]][f"amr class v{VERSION} count"].iat[0]),
    axis=1)

# Create figure 1
plt.figure(figsize=(40, 20))
label_ax = plt.axes((0.05,0.1,0.92,0.42))
feature_ax = plt.axes((0.05,0.55,0.92,0.42), sharex=label_ax)
cb_palette = sn.color_palette("colorblind")

# Make bar plot of median ratio and strip plot of max ratio
sn.barplot(
    data=label_distr_df,
    x="amr class",
    y="ratio",
    hue="label",
    hue_order=["clstr|amr", "super|amr"],
    estimator="median",
    ax=label_ax,
    palette=cb_palette[:2],
    errorbar=None,
    alpha=0.5,
    legend=True)
sn.stripplot(
    data=label_stats,
    x="amr class",
    y="max ratio",
    hue="label",
    hue_order=["clstr|amr", "super|amr"],
    ax=label_ax,
    dodge=True,
    size=15,
    palette=cb_palette[:2],
    legend=True)
label_ax.set_xticks(
    ticks=label_ax.get_xticks(), 
    labels=label_ax.get_xticklabels(), 
    fontsize=30, 
    rotation_mode="anchor", 
    rotation=45,
    ha='right',
    va='center')
label_ax.set_yticks(
    ticks=label_ax.get_yticks(),
    labels=label_ax.get_yticklabels(),
    fontsize=30,
    va='center')
label_ax.set_ylim(
    bottom=0, top=1.1)
label_ax.set_xlabel(
    xlabel="AMR Class",
    fontsize=35,
    loc='center',
    labelpad=0)
label_ax.set_ylabel(ylabel="")
plt.gcf().text(
    0.008, 0.31, f"Sequence Count Ratio", rotation=90, fontsize=35, va='center')
label_ax.set_title(
    "B",
    loc="left",
    fontsize=35,
    weight='bold')
label_ax.legend(
    handles=label_ax.legend_.legend_handles,
    labels=["cluster median", "superfamily median", "cluster max", "superfamily max"],
    fontsize=30,
    loc="upper left")

# Make histogram of sequence counts
hist = sn.histplot(
    data=feature,
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
    fontsize=27)
feature_ax.set_yticks(
    ticks=feature_ax.get_yticks(),
    labels= feature_ax.get_yticklabels(),
    fontsize=30)
feature_ax.set_ylabel(
    ylabel="",
    fontsize=35)
plt.gcf().text(
    0.008, 0.77, "Sequence Count", rotation=90, fontsize=35, va='center')
feature_ax.tick_params(
    axis='x', 
    labelbottom=False)
feature_ax.set_xlabel(
    xlabel="")
feature_ax.set_title(
    "A",
    loc="left",
    fontsize=35,
    weight='bold')
label_ax.margins(x=0)
plt.savefig(f"db_distr_v{VERSION}.png")

# Now create adjacent matrix and edge tuples for shared clusters and superfamily
class_per_clstr_count_df = (pd.DataFrame(
        clstr_class_count_df["clstr"].value_counts(dropna=True, sort=False))
    .reset_index(names="clstr"))
multi_class_clstr = class_per_clstr_count_df.loc[class_per_clstr_count_df["count"] > 1]["clstr"]
multi_class_clstr_rows = (clstr_class_count_df.loc[clstr_class_count_df["clstr"]
    .apply(lambda x: x in multi_class_clstr.to_list())])
clstr_share_matrix = make_label_share_matrix(
    multi_class_clstr_rows, sorted_class, "clstr")
clstr_adj_matrix = make_label_adj_matrix(
    multi_class_clstr_rows, sorted_class, "clstr")

class_per_super_count_df = (pd.DataFrame(
        combo_superfamily_class_count_df["superfamily"].value_counts(dropna=True, sort=False))
    .reset_index(names="superfamily"))
multi_class_supers = class_per_super_count_df.loc[class_per_super_count_df["count"] > 1]["superfamily"]
multi_class_super_rows = (combo_superfamily_class_count_df.loc[combo_superfamily_class_count_df["superfamily"]
        .apply(lambda x: x in multi_class_supers.to_list())])
super_share_matrix = make_label_share_matrix(
    multi_class_super_rows, sorted_class, "superfamily")
super_adj_matrix = make_label_adj_matrix(
    multi_class_super_rows, sorted_class, "superfamily")

# Start with edge and node graph
plt.figure(figsize=(15, 30))
super_ax = plt.axes((0.02,0.01,0.96,0.47))
clstr_ax = plt.axes((0.02,0.50,0.96,0.47))

G_important_nodes = ["PEP", "UNC", "MDR"]
G_clstr_outer_nodes = [
    "GlyP", "MLS", "FQ", "TET", "AC", "AD", "FFA", "TRI", "BL", "PMX", "AG", "PHE"]
G_super_outer_nodes = [
    "PHE", "GlyP", "MLS", "FQ", "TET", "OXA", "AC", "AD", "BCM", "FFA", "FOM", "NUC", "TRI", "BL", "PMX", "BAC", "TET-C", "AG"]

G_clstr_nodes = set([row[0][0] for row in clstr_share_matrix.iterrows() if (
    (row[1].iat[0] > 0) and ((row[0][0] in G_important_nodes) or (row[0][1] in G_important_nodes)))])
G_clstr = nx.Graph()
G_clstr.add_nodes_from(G_clstr_nodes)
G_clstr.add_weighted_edges_from([
    row for row in clstr_share_matrix.reset_index(names=["class 1", "class 2"]).itertuples(index=False,name=None)
    if (row[2] > 0) and (row[0] > row[1]) and ((row[0] in G_important_nodes) or (row[1] in G_important_nodes))])
nx.draw(
    G_clstr, 
    nodelist=G_important_nodes + G_clstr_outer_nodes,
    pos=nx.shell_layout(G_clstr, [
        G_important_nodes,
        G_clstr_outer_nodes]),
    with_labels=True, 
    ax=clstr_ax, 
    font_color="white", 
    font_size=30, 
    node_size=7000,
    node_color=[
        "#D20A2E", "#D20A2E", "#D20A2E",
        "#0F52BA", "#0F52BA", "#0F52BA", "#0F52BA", "#0F52BA", "#0F52BA",
        "#0F52BA", "#0F52BA", "#0F52BA", "#0F52BA", "#0F52BA", "#0F52BA"])
nx.draw_networkx_edge_labels(
    G_clstr, 
    pos=nx.shell_layout(G_clstr, [
        G_important_nodes,
        G_clstr_outer_nodes]),
    edge_labels=nx.get_edge_attributes(G_clstr, 'weight'),
    font_size=30,
    label_pos=0.5,
    ax=clstr_ax)
clstr_ax.set_title("A",
    loc="left",
    fontsize=35,
    weight='bold')

G_super_nodes = set([row[0][0] for row in super_share_matrix.iterrows() if (
    (row[1].iat[0] > 0) and ((row[0][0] in G_important_nodes) or (row[0][1] in G_important_nodes)))])
G_super = nx.Graph()
G_super.add_nodes_from(G_super_nodes)
G_super.add_weighted_edges_from([
    row for row in super_share_matrix.reset_index(names=["class 1", "class 2"]).itertuples(index=False,name=None)
    if (row[2] > 0) and (row[0] > row[1]) and ((row[0] in G_important_nodes) or (row[1] in G_important_nodes))])
nx.draw(
    G_super, 
    nodelist=G_important_nodes + G_super_outer_nodes,
    pos=nx.shell_layout(G_super, [
        G_important_nodes,
        G_super_outer_nodes]),
    with_labels=True, 
    ax=super_ax, 
    font_color="white", 
    font_size=30, 
    node_size=7000,
    node_color=[
        "#D20A2E", "#D20A2E", "#D20A2E",
        "#0F52BA", "#0F52BA", "#0F52BA", "#0F52BA", "#0F52BA", "#0F52BA",
        "#0F52BA", "#0F52BA", "#0F52BA", "#0F52BA", "#0F52BA", "#0F52BA",
        "#0F52BA", "#0F52BA", "#0F52BA", "#0F52BA", "#0F52BA", "#0F52BA"])
nx.draw_networkx_edge_labels(
    G_super, 
    pos=nx.shell_layout(G_super, [
        G_important_nodes,
        G_super_outer_nodes]),
    edge_labels=nx.get_edge_attributes(G_super, 'weight'),
    font_size=30,
    label_pos=0.5,
    ax=super_ax)
super_ax.set_title(
    "B",
    loc="left",
    fontsize=35,
    weight='bold')
plt.savefig(f"inter_class_sim_graph_v{VERSION}.png")

# Create a custom color palette
custom_colors = [
'#D20A2E', # first color
'#d9d9d9', # middle color
'#0F52BA' # last color
]
custom_cmap = LinearSegmentedColormap.from_list("custom_gradient", custom_colors)
custom_cmap.set_bad(color="white")

# Now do heatmap
plt.figure(figsize=(40, 20))
clstr_ax = plt.axes((0.05,0.07,0.4,0.86))
super_ax = plt.axes((0.55,0.07,0.4,0.86))
cbar = plt.axes((0.95, 0.15, 0.02, 0.7))

sn.heatmap(
    data=clstr_adj_matrix,
    #mask=np.triu(np.ones_like(clstr_adj_matrix)),
    vmax=9,
    annot=True,
    center=0,
    cmap=custom_cmap,
    ax=clstr_ax,
    cbar_ax=cbar)
clstr_ax.set_title(
    label="Shared clusters count",
    fontsize=40)
clstr_ax.set_xticks(
    ticks=clstr_ax.get_xticks()[:-1],
    labels=clstr_ax.get_xticklabels()[:-1],
    rotation_mode='anchor',
    rotation=45,
    ha='right',
    va='top',
    fontsize=30)
clstr_ax.set_yticks(
    ticks=clstr_ax.get_yticks()[1:],
    labels=clstr_ax.get_yticklabels()[1:],
    rotation_mode='anchor',
    rotation=0,
    ha='right',
    va='center',
    fontsize=30)
clstr_ax.set_xlabel("")
clstr_ax.set_ylabel("")

sn.heatmap(
    data=super_adj_matrix,
    #mask=np.triu(np.ones_like(super_adj_matrix)),
    vmax=9,
    annot=True,
    center=0,
    cmap=custom_cmap,
    ax=super_ax,
    cbar_ax=cbar)
super_ax.set_title(
    label="Shared superfamilies count",
    fontsize=40)
super_ax.set_xticks(
    ticks=super_ax.get_xticks()[:-1],
    labels=super_ax.get_xticklabels()[:-1],
    rotation_mode='anchor',
    rotation=45,
    ha='right',
    va='top',
    fontsize=30)
super_ax.set_yticks(
    ticks=super_ax.get_yticks()[1:],
    labels=super_ax.get_yticklabels()[1:],
    rotation_mode='anchor',
    rotation=0,
    ha='right',
    va='center',
    fontsize=30)
super_ax.set_xlabel("")
super_ax.set_ylabel("")

cbar.tick_params(labelsize=30)

plt.savefig(f"inter_class_sim_adj_v{VERSION}.png")