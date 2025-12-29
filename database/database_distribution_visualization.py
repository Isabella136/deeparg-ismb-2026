from matplotlib import axes
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
        .value_counts(dropna=False, sort=False))
    .reset_index(names=["amr class", "domain"]))

# Get arg|class count for DeepARG database from feature data
arg_class_count_df = feature[["amr class"]]
arg_class_count_df.insert(0, "arg", feature_args)
arg_class_count_df = (
    pd.DataFrame(
        arg_class_count_df.value_counts(sort=False))
    .reset_index(names=["arg", "amr class"]))

# Create a grid to show label distributions in database
cb_palette = sn.color_palette("colorblind")
grid = sn.JointGrid()
grid.figure.set_figwidth(40)
grid.figure.set_figheight(15)

# For the main joint axes, make a strip plot of arg and domain size distribution per class
# This means that we need to combine arg_class_count and combo_domain_class_count
# s.t. we know which label is arg and which is domain:
arg_class_count_df.insert(
    loc=0,
    column="label",
    value="arg")
combo_domain_class_count_df.insert(
    loc=0,
    column="label",
    value="domain")
sorted_class = (class_count_df
    .sort_values(
        by="amr class v2 count", 
        ascending=False, 
        ignore_index=True)["amr class"])
sort_class_by_size = np.vectorize(lambda x: sorted_class[sorted_class==x].index[0])
label_distr_df = pd.DataFrame(
    data= pd.concat([
        arg_class_count_df,
        combo_domain_class_count_df], join='inner'),
    columns=["label", "amr class", "count"]).sort_values(
        by="amr class",
        key=sort_class_by_size)
sn.stripplot(
    data=label_distr_df,
    x="amr class",
    y="count",
    hue="label",
    ax=grid.ax_joint,
    size=10.0,
    dodge=True,
    alpha=.3,
    palette=cb_palette[:2],
    log_scale=10)
legend_handles = grid.ax_joint.get_legend_handles_labels()[0]
for handle in legend_handles:
    handle.set_alpha(1)
grid.ax_joint.legend(
    handles=legend_handles,
    labels=["arg|amr", "domain|amr"],
    fontsize=20)
grid.ax_joint.set_xticks(
    ticks=grid.ax_joint.get_xticks(), 
    labels=grid.ax_joint.get_xticklabels(), 
    fontsize=20, 
    rotation_mode="anchor", 
    rotation=45,
    ha='right',
    va='center')
grid.ax_joint.set_yticks(
    ticks=grid.ax_joint.get_yticks(),
    labels=grid.ax_joint.get_yticklabels(),
    fontsize=20,
    va='center')
grid.ax_joint.set_ylim(bottom=1, top=1000)
grid.ax_joint.set_xlabel(
    xlabel="AMR Class",
    fontsize=30,
    loc='center')
grid.ax_joint.set_ylabel(
    ylabel="Label Count (log)",
    fontsize=30,
    loc='center')

# For the y marginal axes, make a histogram for domain|amr and arg|amr
sn.histplot(
    data=label_distr_df,
    y="count",
    hue="label",
    ax=grid.ax_marg_y,
    legend=False,
    alpha=0.5,
    stat="density",
    common_norm=False,
    palette=cb_palette[:2])

# For the x marginal axes, make a histogram for domain per amr and arg per amr
sn.histplot(
    data=label_distr_df,
    x="amr class",
    hue="label",
    ax=grid.ax_marg_x,
    legend=False,
    alpha=0.5,
    stat="density",
    common_norm=False,
    multiple="dodge",
    shrink=.8,
    palette=cb_palette[:2])
grid.ax_marg_x.set_yticks(
    ticks=grid.ax_marg_x.get_yticks(),
    labels=[])
grid.ax_marg_x.tick_params(
    axis='y',
    length=0)

# Additionally, add cumulative ecdf plot for amr in x marginal axes
# This means twining the axes
ecdf_ax = grid.ax_marg_x.twinx()
ecdf_ax.set_axis_off()
sn.ecdfplot(
    x=(feature
        .reset_index(drop=True)[["amr class", "amr class v2 count"]]
        .sort_values(by="amr class v2 count", ascending=False)
        .reset_index(drop=True)["amr class"]), 
    complementary=True, 
    ax=ecdf_ax,
    color=cb_palette[2],
    linewidth=2,
    legend=True
)

grid.savefig("db_distr.png")

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