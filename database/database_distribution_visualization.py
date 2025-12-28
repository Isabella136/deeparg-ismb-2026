import matplotlib
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
get_arg = np.vectorize(lambda a: a.split('|')[-1])
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


sankey_label = pd.concat([
    combo_domain_count_df["domain"],
    class_count_df["amr class"],
    arg_count_df["arg"]], ignore_index=True).fillna("no domain")
sankey_label.to_csv("test.csv")

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