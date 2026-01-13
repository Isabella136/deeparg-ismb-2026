import matplotlib
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import plotly.graph_objects as go

import seaborn as sn
import pandas as pd
import numpy as np

import sys

FEATURE_DATA = "v2_feature_data.csv"

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
    .reset_index(drop=True)[["amr class", "amr class v2 count"]]
    .drop_duplicates()
    .reset_index(drop=True))

# Using code from good old stackoverflow: 
    # https://stackoverflow.com/questions/71688904/dealing-with-multiple-values-in-pandas-dataframe-cell
# Get domain count for DeepARG database from feature data 
# (Counting each domain in multi-dom features)
indiv_domain_count_df = (feature
    .reset_index()[["index", "conserved domain(s) id(s)"]]
    .melt("index")[["variable", "value"]])
indiv_domain_count_df["value"] = indiv_domain_count_df["value"].str.split('$')
indiv_domain_count_df = (
    pd.DataFrame(
        indiv_domain_count_df.explode("value")["value"].value_counts(dropna=False, sort=False))
    .reset_index(names="domain"))

# Get domain count for DeepARG database from feature data 
# (Counting combos in multi-dom features as their own domain)
combo_domain_count_df = (
    pd.DataFrame(
        feature.reset_index()[["conserved domain(s) id(s)"]].value_counts(dropna=False, sort=False))
    .reset_index(names="domain"))

# Get arg count for DeepARG database from feature data
get_arg = np.vectorize(lambda a: a.split('|')[-1].lower())
feature_args: np.ndarray = get_arg(feature.index.to_numpy())
arg_count_df = np.unique_counts(feature_args)
arg_count_df = pd.DataFrame({"arg": arg_count_df.values, "count": arg_count_df.counts})

# Get clstr|class count for DeepARG database from feature data
clstr_class_count_df = (pd.DataFrame(feature
        .reset_index()[["amr class", "v2-only cluster index"]]
        .value_counts(dropna=False, sort=False))
    .reset_index(names=["amr class", "clstr"]))

# Using code from good old stackoverflow: 
    # https://stackoverflow.com/questions/71688904/dealing-with-multiple-values-in-pandas-dataframe-cell
# Get domain|class count for DeepARG database from feature data
# (Counting each domain in multi-dom features)

indiv_domain_class_count_df = (feature
    .reset_index()[["index", "amr class", "conserved domain(s) id(s)"]]
    .melt("index"))
indiv_domain_class_count_df["value"] = indiv_domain_class_count_df["value"].str.split('$')
indiv_domain_class_count_df = indiv_domain_class_count_df.explode("value")
corresponding_class_df = indiv_domain_class_count_df.loc[
    indiv_domain_class_count_df["variable"] == "amr class"]
corresponding_class_df = corresponding_class_df[["index", "value"]].set_index("index")
indiv_domain_class_count_df = (
    pd.DataFrame(indiv_domain_class_count_df[["index", "value"]]
        .apply(lambda x: pd.Series({
                "domain": x["value"],
                "amr class": corresponding_class_df.at[x["index"],"value"]
            }), axis=1)
        .value_counts(dropna=False, sort=False)))

# Get domain|class count for DeepARG database from feature data
# (Counting combos in multi-dom features as their own domain)
combo_domain_class_count_df = (pd.DataFrame(feature
        .reset_index()[["amr class", "conserved domain(s) id(s)"]]
        .value_counts(dropna=True, sort=False))
    .reset_index(names=["amr class", "domain"]))

# Get arg|class count for DeepARG database from feature data
arg_class_count_df = feature[["amr class"]]
arg_class_count_df.insert(0, "arg", feature_args)
arg_class_count_df = (
    pd.DataFrame(
        arg_class_count_df.value_counts(sort=False))
    .reset_index(names=["arg", "amr class"]))

print(combo_domain_class_count_df.columns)

clstr_stats = pd.DataFrame(
    data={
        "average": clstr_class_count_df.groupby("amr class")["count"].mean(),
        "max": clstr_class_count_df.groupby("amr class")["count"].max(),
        "min": clstr_class_count_df.groupby("amr class")["count"].min()},
    index=clstr_class_count_df["amr class"].drop_duplicates().sort_values())
clstr_stats.reset_index(names="amr class", inplace=True)
clstr_stats.insert(
    loc=0,
    column="label",
    value="clstr|amr")
combo_domain_stats = pd.DataFrame(
    data={
        "average": combo_domain_class_count_df.groupby(by="amr class")["count"].mean(),
        "max": combo_domain_class_count_df.groupby(by="amr class")["count"].max(),
        "min": combo_domain_class_count_df.groupby(by="amr class")["count"].min()},
    index=combo_domain_class_count_df["amr class"].drop_duplicates().sort_values())
combo_domain_stats.reset_index(names="amr class", inplace=True)
combo_domain_stats.insert(
    loc=0,
    column="label",
    value="domain|amr")
sorted_class = (class_count_df
    .sort_values(
        by="amr class v2 count", 
        ascending=False, 
        ignore_index=True)["amr class"])
sort_class_by_size = np.vectorize(lambda x: sorted_class[sorted_class==x].index[0])
label_stats = pd.DataFrame(
    data= pd.concat([
        clstr_stats,
        combo_domain_stats], join='inner'),
    columns=["label", "amr class", "average", "max", "min"]).sort_values(
        by="amr class",
        key=sort_class_by_size)
label_stats["max ratio"] = label_stats.apply(
    lambda x: float(x["max"])/float(
        class_count_df.loc[class_count_df["amr class"]==x["amr class"]]["amr class v2 count"].iat[0]),
    axis=1)
label_stats["min ratio"] = label_stats.apply(
    lambda x: float(x["min"])/float(
        class_count_df.loc[class_count_df["amr class"]==x["amr class"]]["amr class v2 count"].iat[0]),
    axis=1)
label_stats.to_csv("label_stats.csv")

# Create a grid to show label distributions in database
cb_palette = sn.color_palette("colorblind")
custom_colors = [
    cb_palette[0], # first color
    cb_palette[1], # second color
    '#000000' # last color
]

plt.figure(figsize=(40, 20))
label = plt.axes((0.05,0.1,0.92,0.42))
feature_ax = plt.axes((0.05,0.55,0.92,0.42), sharex=label)

# We need to combine arg_class_count and combo_domain_class_count
# s.t. we know which label is arg and which is domain:
clstr_class_count_df.insert(
    loc=0,
    column="label",
    value="clstr|amr")
combo_domain_class_count_df.insert(
    loc=0,
    column="label",
    value="domain|amr")
label_distr_df = pd.DataFrame(
    data= pd.concat([
        clstr_class_count_df,
        combo_domain_class_count_df], join='inner'),
    columns=["label", "amr class", "count"]).sort_values(
        by="amr class",
        key=sort_class_by_size)
label_distr_df["ratio"] = label_distr_df.apply(
    lambda x: float(x["count"])/float(
        class_count_df.loc[class_count_df["amr class"]==x["amr class"]]["amr class v2 count"].iat[0]),
    axis=1)
label_distr_df.to_csv("label_distr_df.csv")

sn.barplot(
    data=label_distr_df,
    x="amr class",
    y="ratio",
    hue="label",
    hue_order=["clstr|amr", "domain|amr"],
    estimator="median",
    ax=label,
    palette=cb_palette[:2],
    errorbar=None,
    alpha=0.5,
    legend=True)
    #log_scale=10)
# sn.stripplot(
#     data=label_stats,
#     x="amr class",
#     y="min ratio",
#     hue="label",
#     hue_order=["clstr|amr", "domain|amr"],
#     ax=label,
#     dodge=True,
#     size=15,
#     palette=cb_palette[:2],
#     legend=True,
#     marker="X")
sn.stripplot(
    data=label_stats,
    x="amr class",
    y="max ratio",
    hue="label",
    hue_order=["clstr|amr", "domain|amr"],
    ax=label,
    dodge=True,
    size=15,
    palette=cb_palette[:2],
    legend=True)
label.set_xticks(
    ticks=label.get_xticks(), 
    labels=label.get_xticklabels(), 
    fontsize=25, 
    rotation_mode="anchor", 
    rotation=45,
    ha='right',
    va='center')
label.set_yticks(
    ticks=label.get_yticks(),
    labels=label.get_yticklabels(),
    fontsize=25,
    va='center')
# label.set_ylim(bottom=1, top=1000)
label.set_ylim(bottom=0, top=1.1)
label.set_xlabel(
    xlabel="AMR Class",
    fontsize=35,
    loc='center',
    labelpad=0)
label.set_ylabel(
    ylabel="Average Feature Count Ratio",
    fontsize=35,
    loc='center')
label.set_title(
    "B",
    loc="left",
    fontsize=35)
label.legend(
    handles=label.legend_.legend_handles,
    #labels=["clstr|amr", "domain|amr", "clstr|amr min", "domain|amr min", "clstr|amr max", "domain|amr max"],
    labels=["clstr|amr", "domain|amr", "clstr|amr max", "domain|amr max"],
    fontsize=25)

sn.histplot(
    data=feature,
    x="amr class",
    ax=feature_ax,
    legend=True,
    fill=True,
    linewidth=3,
    element="bars",
    common_norm=False,
    shrink=0.8,
    color="black")
feature_ax.set_yticks(
    ticks=feature_ax.get_yticks(),
    labels= feature_ax.get_yticklabels(),
    fontsize=25)
feature_ax.set_ylabel(
    ylabel="Feature Count",
    fontsize=35)
feature_ax.tick_params(
    axis='x', 
    labelbottom=False)
feature_ax.set_xlabel(
    xlabel="")
feature_ax.set_title(
    "A",
    loc="left",
    fontsize=35)

label.margins(x=0)

plt.savefig("db_distr.png")

sys.exit()

# Should restrict sankey to divergent labels
class_per_dom_count_df = (
    pd.DataFrame(
        combo_domain_class_count_df["domain"].value_counts(dropna=False, sort=False))
    .reset_index(names="domain"))
multi_class_doms = class_per_dom_count_df.loc[
    class_per_dom_count_df["count"] > 1]["domain"].fillna("no domain")
multi_class_dom_rows = (combo_domain_class_count_df
    .fillna("no domain")
    .loc[combo_domain_class_count_df["domain"]
        .fillna("no domain")
        .apply(lambda x: x in multi_class_doms.to_list())])

sankey_label = pd.concat([
    multi_class_dom_rows["amr class"],
    multi_class_dom_rows["domain"]], ignore_index=True)

sankey_source = multi_class_dom_rows["amr class"].apply(
    lambda x: sankey_label[sankey_label == x].index[0])

sankey_target = multi_class_dom_rows["domain"].apply(
    lambda x: sankey_label[sankey_label == x].index[0])

sankey_value = multi_class_dom_rows["count"]

sankey = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color="black", width=0.5),
        label=sankey_label,
        color="blue"),
    link=dict(
        source=sankey_source,
        target=sankey_target,
        value=sankey_value)
    )])

sankey.write_html("test.html")

sys.exit()

sankey_label = pd.concat([
    combo_domain_count_df["domain"],
    class_count_df["amr class"],
    arg_count_df["arg"]], ignore_index=True).fillna("no domain")

sankey_source = pd.concat([
    combo_domain_class_count_df["domain"],
    arg_class_count_df["amr class"]]).fillna("no domain").apply(
        lambda x: sankey_label[sankey_label == x].index[0])

sankey_target = pd.concat([
    combo_domain_class_count_df["amr class"],
    arg_class_count_df["arg"]]).apply(
        lambda x: sankey_label[sankey_label == x].index[0])

sankey_value = pd.concat([
    combo_domain_class_count_df["count"],
    arg_class_count_df["count"]])

sankey = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color="black", width=0.5),
        label=sankey_label,
        color="blue"),
    link=dict(
        source=sankey_source,
        target=sankey_target,
        value=sankey_value)
    )])

sankey.write_html("test.html")