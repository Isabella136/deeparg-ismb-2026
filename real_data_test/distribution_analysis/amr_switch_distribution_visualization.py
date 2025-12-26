import matplotlib.axes
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sn
import pandas as pd
import numpy as np
from io import StringIO
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
import sys

TABLE_LOCATION = "differential_distribution_output"
FEATURE_DATA = "../../database/v2_feature_data.csv"
HIT_COUNT = "../samples/deeparg_hit_count.tsv"

matplotlib.rcParams["mathtext.fontset"] = "dejavusans"
matplotlib.rcParams["font.family"] = "DejaVu Sans"

def get_hit_count(diff_count: pd.DataFrame, x: pd.Series) -> int:
    try:
        return diff_count.loc[(x["sample id"], x["alignment identity"], x["model"]), "switch pair count"]
    except:
        return 0

abbrev_file = open("amr_abbrev.csv", "r")
amr_abbrev = dict()
for line in abbrev_file.readlines():
    row = line.strip().split(",")
    amr_abbrev.update({row[0]: row[1]})
abbrev_file.close()

# Get class count for DeepARG database by retrieving feature data
feature = pd.read_csv(FEATURE_DATA, header=0)
class_count_df = feature.drop_duplicates(subset=["amr class", "amr class v2 count"])
class_count_dict = dict(zip(class_count_df["amr class"], class_count_df["amr class v2 count"]))

df = pd.read_csv(f"{TABLE_LOCATION}/switch_pair_abundance_amr.tsv", sep='\t', header=0)
diff_count = df[[
    "sample", "alignment identity", "model", "diamond best-hit label",
    "deeparg hit label", "switch pair count"]].drop_duplicates()
diff_count = diff_count.groupby(by=["sample", "alignment identity", "model"])[["switch pair count"]].sum()

hit_count_df = pd.read_csv(HIT_COUNT, sep="\t", header=0, index_col=0)
hit_count_df.insert(4, "diff count", 0)
hit_count_df["diff count"] = hit_count_df.apply(
    lambda x: get_hit_count(diff_count, x), axis=1
)
hit_count_df = hit_count_df.sort_values(by=["model", "sample id", "alignment identity"])
hit_count_df.insert(5, "hit percentage", 0)
hit_count_df["hit percentage"] = hit_count_df.apply(
    lambda x: float(x["hit count"]) / float(hit_count_df.loc[
        (hit_count_df["sample id"] == x["sample id"]) & 
        (hit_count_df["alignment identity"] == 30) & 
        (hit_count_df["model"] == x["model"])].iat[0, 3]), axis=1)
hit_count_df.insert(6, "diff percentage", 0)
hit_count_df["diff percentage"] = hit_count_df.apply(
    lambda x: float(x["diff count"]) / float(hit_count_df.loc[
        (hit_count_df["sample id"] == x["sample id"]) & 
        (hit_count_df["alignment identity"] == 30) & 
        (hit_count_df["model"] == x["model"])].iat[0, 3]), axis=1)

cb_palette = sn.color_palette("colorblind")
plt.figure(figsize=(40, 12))
axes: list[matplotlib.axes.Axes] = [plt.axes((0.05,0.12,0.46,0.76)), plt.axes((0.53,0.12,0.46,0.76))]

for model, ax in zip(["LS", "SS"], axes):

    filtered_hit_count_df = hit_count_df.loc[hit_count_df["model"] == model]

    same = sn.barplot(
        x=range(39),
        y=np.insert(
            arr=filtered_hit_count_df["hit percentage"].values,
            obj=slice(3,30,3),
            values=np.full(9, 0)), 
        color=cb_palette[0],
        ax=ax)
    different = sn.barplot(
        x=range(39),
        y=np.insert(
            arr=filtered_hit_count_df["diff percentage"].values,
            obj=slice(3,30,3),
            values=np.full(9, 0)), 
        color=cb_palette[1],
        ax=ax)

    if model == "SS":
        top_bar = mpatches.Patch(color=cb_palette[0], label='same')
        bottom_bar = mpatches.Patch(color=cb_palette[1], label='different')
        ax.legend(handles=[top_bar, bottom_bar], fontsize=24)
    ax.set_xticks(
        ticks=list(range(39)), 
        labels=np.insert(
            arr=filtered_hit_count_df["alignment identity"].values.astype('<U2'),
            obj=slice(3,30,3),
            values=np.full(9, "")),
        rotation_mode="anchor", ha='center', fontsize=24)
    for loc, sample in enumerate(filtered_hit_count_df["sample id"].drop_duplicates().values):
        print(sample)
        ax.text(
            loc*4+1, -0.075, f"#{loc+1}", ha="center", va="center", fontsize=28)
    ax.set_ylim(top=1.0)
    if model == "LS":
        ax.set_yticks(
            ticks=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1], 
            labels=list(range(0, 110, 10)), fontsize=24)
        ax.set_ylabel(f"Percentage of total DeepARG 30% hits", fontsize=30)
    else:
        ax.set_yticks(ticks=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1], labels="")
    ax.grid(visible=True, which="major", axis='y', color="0.75")
    ax.tick_params(length=0, pad=8)
    ax.set_title(f"DeepARG-{model} Similarity to Diamond Best Hit", fontsize=36, loc="left")

plt.gcf().text(0.52, 0.015, "Run (upper value is alignment identity, lower value is sample)", ha="center", fontsize=30)
plt.savefig("barplot.png")
sys.exit()

# Rename AMR columns based on abbrev
for row in range(df.shape[0]):
    df.at[row, "diamond best-hit label"] = amr_abbrev[df.at[row, "diamond best-hit label"]]
    df.at[row, "deeparg hit label"] = amr_abbrev[df.at[row, "deeparg hit label"]]

# Create first figure where we group by Diamond AMR, then second figure where we group by DeepARG AMR
for (group, indiv) in [("Diamond AMR", "DeepARG AMR"), ("DeepARG AMR", "Diamond AMR")]:
    # Sort DataFrame 
    df.sort_values(
        by=["Model", "Alignment Identity", "Sample", group, indiv], 
        inplace=True,
        ignore_index=True)

    # Make heatmap multilevel x-axis indices
    indices = pd.MultiIndex.from_arrays(
        [df.get(["Model", "Alignment Identity", "Sample"]).drop_duplicates()["Model"].values,
        df.get(["Model", "Alignment Identity", "Sample"]).drop_duplicates()["Alignment Identity"].values,
        df.get(["Model", "Alignment Identity", "Sample"]).drop_duplicates()["Sample"].values],
        names=["Model", "Alignment Identity", "Sample"])

    # Make heatmap multilevel y-axis headers
    headers = pd.MultiIndex.from_arrays(
        [df.get(["Diamond AMR", "DeepARG AMR"]).drop_duplicates().sort_values(
            by=[group, indiv], key=lambda x: x.str.lower())[group].values,
        df.get(["Diamond AMR", "DeepARG AMR"]).drop_duplicates().sort_values(
            by=[group, indiv], key=lambda x: x.str.lower())[indiv].values],
        names=[group, indiv])

    heatmap_df = pd.DataFrame(
        np.empty(
            shape=[45, 97]), 
        index=indices, columns=headers)

    mask_df = pd.DataFrame(
        np.full(
            shape=[45, 97], fill_value=True), 
        index=indices, columns=headers)

    max_val = 0
    min_val = 0

    for row in df.iterrows() :
        model = row[1]["Model"]
        ident = row[1]["Alignment Identity"]
        sampl = row[1]["Sample"]
        diamo = row[1]["Diamond AMR"]
        deepa = row[1]["DeepARG AMR"]
        group_curr = row[1][group]
        indiv_curr = row[1][indiv]

        heatmap_df.at[(model, ident, sampl), (group_curr, indiv_curr)] = np.log(
            (1 + row[1]["Count"])/(
                1 + (class_count_dict[deepa] / (sum(class_count_dict.values()) - class_count_dict[diamo]) * sum(
                    df.loc[df["Model"] == model].loc[df["Alignment Identity"] == ident].loc[df["Sample"] == sampl].loc[df["Diamond AMR"] == diamo]["Count"].values))))
        mask_df.at[(model, ident, sampl), (group_curr, indiv_curr)] = False
        max_val = max(max_val, heatmap_df.at[(model, ident, sampl), (group_curr, indiv_curr)])
        min_val = min(min_val, heatmap_df.at[(model, ident, sampl), (group_curr, indiv_curr)])

    # Create a custom color palette
    custom_colors = [
    '#D20A2E', # first color
    '#d9d9d9', # middle color
    '#0F52BA' # last color
    ]
    custom_cmap = LinearSegmentedColormap.from_list("custom_gradient", custom_colors)
    custom_cmap.set_bad()

    # Create heatmap and colormap key
    fig = plt.figure.Figure(figsize=[32, 15])
    ax = fig.add_axes((0.07, 0.05, 0.85, 0.9))
    cbar = fig.add_axes((0.94, 0.2, 0.02, 0.6))
    heatmap = sn.heatmap(
        heatmap_df, mask=mask_df, cmap=custom_cmap, center=0, vmax=max_val, vmin=min_val, 
        square=True, ax=ax, cbar_ax=cbar,
        xticklabels=heatmap_df.columns.get_level_values(indiv),
        yticklabels=heatmap_df.index.get_level_values("Sample"))
    cbar.tick_params(labelsize=16)
    
    # Set x-ticks direction and axis label depending on current grouped category
    if group == "Diamond AMR":
        ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation_mode="anchor", rotation=45, ha='right', fontsize=14)
        ax.set_xlabel(indiv, fontsize=20, color="#ffaa5b", weight="bold", style="italic")
    else:
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top') 
        ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation_mode="anchor", rotation=45, ha='left', fontsize=14)
        ax.set_xlabel(r'$\underline{\textbf{\textsf{\Huge{Diamond\ AMR}}}}$', color="#3dac87", usetex=True)

    # Remove y-axis ticks
    ax.yaxis.set_ticks([])
    ax.set_ylabel("")

    # Seperate by group
    y_groups = heatmap_df.index.to_frame(index=False).groupby(["Model", "Alignment Identity"]).size()
    x_groups = heatmap_df.columns.to_frame(index=False).groupby([group], sort=False).size()

    y_boundaries = np.cumsum(y_groups)
    x_boundaries = np.cumsum(x_groups)

    for y_loc in y_boundaries[:-1]:
        ax.hlines(y_loc, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], colors="white", linewidth=5)
    for x_loc in x_boundaries[:-1]:
        ax.vlines(x_loc, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], colors="white", linewidth=5)

    # Add group tick labels
    prev_y_loc = 0
    for y_loc, y_group in zip(y_boundaries, y_groups.index):
        ax.text(-0.5, (y_loc+prev_y_loc)/2.0, y_group, ha="right", fontsize=16)
        prev_y_loc = y_loc

    prev_x_loc = 0
    for x_loc, x_group in zip(x_boundaries, x_groups.index):
        if group == "Diamond AMR":
            ax.text((x_loc+prev_x_loc)/2.0, -0.5, x_group, rotation=45, rotation_mode="anchor", va="center", fontsize=14)
        else:
            ax.text((x_loc+prev_x_loc)/2.0, y_boundaries.iloc[-1] + 0.5, x_group, rotation=45, rotation_mode="anchor", va="center", ha="right", fontsize=14)
        prev_x_loc = x_loc

    # Add group axis labels
    ax.text(-7, y_boundaries.iloc[-1]/2.0, "(Model, Alignment Identity)", rotation=90, va="center", fontsize=20)
    if group == "Diamond AMR":
        ax.text(
            x_boundaries.iloc[-1]/2.0, -2.5, r'$\underline{\textbf{\textsf{\Huge{Diamond\ AMR}}}}$', 
            rotation=0, ha="center", color="#3dac87", usetex=True)
    else:
        ax.text(
            x_boundaries.iloc[-1]/2.0, y_boundaries.iloc[-1] + 3.5, group, 
            rotation=0, ha="center", fontsize=20, color="#ffaa5b", weight="bold", style="italic")

    fig.savefig("{a}_amr_switch_distribution.png".format(a = group.split(' ')[0].lower()))


    # Create heatmap and colormap key for cluster amr switches only
    if group == "Diamond AMR":
        heatmap_df_cluster = heatmap_df.loc[:, [
            ("AG", "MDR"), ("AG", "UNC"), ("BL", "MDR"), ("FFA", "MDR"),
            ("FQ", "MDR"), ("GlyP", "UNC"), ("MDR", "AG"), ("MDR", "BL"),
            ("MDR", "FFA"), ("MDR", "FQ"), ("MDR", "MLS"), ("MDR", "PHE"),
            ("MDR", "TRI"), ("MDR", "UNC"), ("MLS", "MDR"), ("MLS", "PHE"),
            ("MLS", "TET"), ("PEP", "PMX"), ("PHE", "MLS"), ("PMX", "PEP"),
            ("TET", "MDR"), ("TET", "MLS"), ("TRI", "MDR"), ("UNC", "GlyP"),
            ("UNC", "MDR"), ("UNC", "PHE")]]
        mask_df_cluster = mask_df.loc[:, [
            ("AG", "MDR"), ("AG", "UNC"), ("BL", "MDR"), ("FFA", "MDR"),
            ("FQ", "MDR"), ("GlyP", "UNC"), ("MDR", "AG"), ("MDR", "BL"),
            ("MDR", "FFA"), ("MDR", "FQ"), ("MDR", "MLS"), ("MDR", "PHE"),
            ("MDR", "TRI"), ("MDR", "UNC"), ("MLS", "MDR"), ("MLS", "PHE"),
            ("MLS", "TET"), ("PEP", "PMX"), ("PHE", "MLS"), ("PMX", "PEP"),
            ("TET", "MDR"), ("TET", "MLS"), ("TRI", "MDR"), ("UNC", "GlyP"),
            ("UNC", "MDR"), ("UNC", "PHE")]]
    else:
        heatmap_df_cluster = heatmap_df.loc[:, [
            ("AG", "MDR"), ("BL", "MDR"), ("FFA", "MDR"), ("FQ", "MDR"),
            ("GlyP", "UNC"), ("MDR", "AG"),  ("MDR", "BL"), ("MDR", "FFA"), 
            ("MDR", "FQ"), ("MDR", "MLS"), ("MDR", "TET"), ("MDR", "TRI"), 
            ("MDR", "UNC"), ("MLS", "MDR"), ("MLS", "PHE"), ("MLS", "TET"),
            ("PEP", "PMX"), ("PHE", "MDR"), ("PHE", "MLS"), ("PHE", "UNC"), 
            ("PMX", "PEP"), ("TET", "MLS"), ("TRI", "MDR"), ("UNC", "AG"), 
            ("UNC", "MDR"), ("UNC", "GlyP"),]]
        mask_df_cluster = mask_df.loc[:, [
            ("AG", "MDR"), ("BL", "MDR"), ("FFA", "MDR"), ("FQ", "MDR"),
            ("GlyP", "UNC"), ("MDR", "AG"),  ("MDR", "BL"), ("MDR", "FFA"), 
            ("MDR", "FQ"), ("MDR", "MLS"), ("MDR", "TET"), ("MDR", "TRI"), 
            ("MDR", "UNC"), ("MLS", "MDR"), ("MLS", "PHE"), ("MLS", "TET"),
            ("PEP", "PMX"), ("PHE", "MDR"), ("PHE", "MLS"), ("PHE", "UNC"), 
            ("PMX", "PEP"), ("TET", "MLS"), ("TRI", "MDR"), ("UNC", "AG"), 
            ("UNC", "MDR"), ("UNC", "GlyP"),]]
        
    fig = plt.figure.Figure(figsize=[12, 15])
    ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))
    cbar = fig.add_axes((0.9, 0.2, 0.05, 0.6))
    heatmap = sn.heatmap(
        heatmap_df_cluster, mask=mask_df_cluster, cmap=custom_cmap, center=0, vmax=max_val, vmin=min_val, 
        square=True, ax=ax, cbar_ax=cbar,
        xticklabels=heatmap_df_cluster.columns.get_level_values(indiv),
        yticklabels=heatmap_df_cluster.index.get_level_values("Sample"))
    cbar.tick_params(labelsize=16)
    
    # Set x-ticks direction and axis label depending on current grouped category
    if group == "Diamond AMR":
        ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation_mode="anchor", rotation=45, ha='right', fontsize=14)
        ax.set_xlabel(indiv, fontsize=20, color="#ffaa5b", weight="bold", style="italic")
    else:
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top') 
        ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation_mode="anchor", rotation=45, ha='left', fontsize=14)
        ax.set_xlabel(r'$\underline{\textbf{\textsf{\Huge{Diamond\ AMR}}}}$', color="#3dac87", usetex=True)

    # Remove y-axis ticks
    ax.yaxis.set_ticks([])
    ax.set_ylabel("")

    # Seperate by group
    y_groups = heatmap_df_cluster.index.to_frame(index=False).groupby(["Model", "Alignment Identity"]).size()
    x_groups = heatmap_df_cluster.columns.to_frame(index=False).groupby([group], sort=False).size()

    y_boundaries = np.cumsum(y_groups)
    x_boundaries = np.cumsum(x_groups)

    for y_loc in y_boundaries[:-1]:
        ax.hlines(y_loc, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], colors="white", linewidth=5)
    for x_loc in x_boundaries[:-1]:
        ax.vlines(x_loc, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], colors="white", linewidth=5)

    # Add group tick labels
    prev_y_loc = 0
    for y_loc, y_group in zip(y_boundaries, y_groups.index):
        ax.text(-0.5, (y_loc+prev_y_loc)/2.0, y_group, ha="right", fontsize=16)
        prev_y_loc = y_loc

    prev_x_loc = 0
    for x_loc, x_group in zip(x_boundaries, x_groups.index):
        if group == "Diamond AMR":
            ax.text((x_loc+prev_x_loc)/2.0, -0.5, x_group, rotation=45, rotation_mode="anchor", va="center", fontsize=14)
        else:
            ax.text((x_loc+prev_x_loc)/2.0, y_boundaries.iloc[-1] + 0.5, x_group, rotation=45, rotation_mode="anchor", va="center", ha="right", fontsize=14)
        prev_x_loc = x_loc

    # Add group axis labels
    ax.text(-7, y_boundaries.iloc[-1]/2.0, "(Model, Alignment Identity)", rotation=90, va="center", fontsize=20)
    if group == "Diamond AMR":
        ax.text(
            x_boundaries.iloc[-1]/2.0, -2, r'$\underline{\textbf{\textsf{\Huge{Diamond\ AMR}}}}$', 
            rotation=0, ha="center", color="#3dac87", usetex=True)
    else:
        ax.text(
            x_boundaries.iloc[-1]/2.0, y_boundaries.iloc[-1] + 3, group, 
            rotation=0, ha="center", fontsize=20, color="#ffaa5b", weight="bold", style="italic")

    fig.savefig("{a}_clustered_amr_switch_distribution.png".format(a = group.split(' ')[0].lower()))

    # Create heatmap and colormap key for non-clustered amr switches only
    if group == "Diamond AMR":
        heatmap_df_noncluster = heatmap_df.drop([
            ("AG", "MDR"), ("AG", "UNC"), ("BL", "MDR"), ("FFA", "MDR"),
            ("FQ", "MDR"), ("GlyP", "UNC"), ("MDR", "AG"), ("MDR", "BL"),
            ("MDR", "FFA"), ("MDR", "FQ"), ("MDR", "MLS"), ("MDR", "PHE"),
            ("MDR", "TRI"), ("MDR", "UNC"), ("MLS", "MDR"), ("MLS", "PHE"),
            ("MLS", "TET"), ("PEP", "PMX"), ("PHE", "MLS"), ("PMX", "PEP"),
            ("TET", "MDR"), ("TET", "MLS"), ("TRI", "MDR"), ("UNC", "GlyP"),
            ("UNC", "MDR"), ("UNC", "PHE")], axis=1)
        mask_df_noncluster = mask_df.drop([
            ("AG", "MDR"), ("AG", "UNC"), ("BL", "MDR"), ("FFA", "MDR"),
            ("FQ", "MDR"), ("GlyP", "UNC"), ("MDR", "AG"), ("MDR", "BL"),
            ("MDR", "FFA"), ("MDR", "FQ"), ("MDR", "MLS"), ("MDR", "PHE"),
            ("MDR", "TRI"), ("MDR", "UNC"), ("MLS", "MDR"), ("MLS", "PHE"),
            ("MLS", "TET"), ("PEP", "PMX"), ("PHE", "MLS"), ("PMX", "PEP"),
            ("TET", "MDR"), ("TET", "MLS"), ("TRI", "MDR"), ("UNC", "GlyP"),
            ("UNC", "MDR"), ("UNC", "PHE")], axis=1)
    else:
        heatmap_df_noncluster = heatmap_df.drop([
            ("AG", "MDR"), ("BL", "MDR"), ("FFA", "MDR"), ("FQ", "MDR"),
            ("GlyP", "UNC"), ("MDR", "AG"),  ("MDR", "BL"), ("MDR", "FFA"), 
            ("MDR", "FQ"), ("MDR", "MLS"), ("MDR", "TET"), ("MDR", "TRI"), 
            ("MDR", "UNC"), ("MLS", "MDR"), ("MLS", "PHE"), ("MLS", "TET"),
            ("PEP", "PMX"), ("PHE", "MDR"), ("PHE", "MLS"), ("PHE", "UNC"), 
            ("PMX", "PEP"), ("TET", "MLS"), ("TRI", "MDR"), ("UNC", "AG"), 
            ("UNC", "MDR"), ("UNC", "GlyP")], axis=1)
        mask_df_noncluster = mask_df.drop([
            ("AG", "MDR"), ("BL", "MDR"), ("FFA", "MDR"), ("FQ", "MDR"),
            ("GlyP", "UNC"), ("MDR", "AG"),  ("MDR", "BL"), ("MDR", "FFA"), 
            ("MDR", "FQ"), ("MDR", "MLS"), ("MDR", "TET"), ("MDR", "TRI"), 
            ("MDR", "UNC"), ("MLS", "MDR"), ("MLS", "PHE"), ("MLS", "TET"),
            ("PEP", "PMX"), ("PHE", "MDR"), ("PHE", "MLS"), ("PHE", "UNC"), 
            ("PMX", "PEP"), ("TET", "MLS"), ("TRI", "MDR"), ("UNC", "AG"), 
            ("UNC", "MDR"), ("UNC", "GlyP")], axis=1)
        
    fig = plt.figure.Figure(figsize=[24, 15])
    ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))
    cbar = fig.add_axes((0.94, 0.2, 0.02, 0.6))
    heatmap = sn.heatmap(
        heatmap_df_noncluster, mask=mask_df_noncluster, cmap=custom_cmap, center=0, vmax=max_val, vmin=min_val, 
        square=True, ax=ax, cbar_ax=cbar,
        xticklabels=heatmap_df_noncluster.columns.get_level_values(indiv),
        yticklabels=heatmap_df_noncluster.index.get_level_values("Sample"))
    cbar.tick_params(labelsize=16)
    
    # Set x-ticks direction and axis label depending on current grouped category
    if group == "Diamond AMR":
        ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation_mode="anchor", rotation=45, ha='right', fontsize=14)
        ax.set_xlabel(indiv, fontsize=20, color="#ffaa5b", weight="bold", style="italic")
    else:
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top') 
        ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation_mode="anchor", rotation=45, ha='left', fontsize=14)
        ax.set_xlabel(r'$\underline{\textbf{\textsf{\Huge{Diamond\ AMR}}}}$', color="#3dac87", usetex=True)

    # Remove y-axis ticks
    ax.yaxis.set_ticks([])
    ax.set_ylabel("")

    # Seperate by group
    y_groups = heatmap_df_noncluster.index.to_frame(index=False).groupby(["Model", "Alignment Identity"]).size()
    x_groups = heatmap_df_noncluster.columns.to_frame(index=False).groupby([group], sort=False).size()

    y_boundaries = np.cumsum(y_groups)
    x_boundaries = np.cumsum(x_groups)

    for y_loc in y_boundaries[:-1]:
        ax.hlines(y_loc, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], colors="white", linewidth=5)
    for x_loc in x_boundaries[:-1]:
        ax.vlines(x_loc, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], colors="white", linewidth=5)

    # Add group tick labels
    prev_y_loc = 0
    for y_loc, y_group in zip(y_boundaries, y_groups.index):
        ax.text(-0.5, (y_loc+prev_y_loc)/2.0, y_group, ha="right", fontsize=16)
        prev_y_loc = y_loc

    prev_x_loc = 0
    for x_loc, x_group in zip(x_boundaries, x_groups.index):
        if group == "Diamond AMR":
            ax.text((x_loc+prev_x_loc)/2.0, -0.5, x_group, rotation=45, rotation_mode="anchor", va="center", fontsize=14)
        else:
            ax.text((x_loc+prev_x_loc)/2.0, y_boundaries.iloc[-1] + 0.5, x_group, rotation=45, rotation_mode="anchor", va="center", ha="right", fontsize=14)
        prev_x_loc = x_loc

    # Add group axis labels
    ax.text(-7, y_boundaries.iloc[-1]/2.0, "(Model, Alignment Identity)", rotation=90, va="center", fontsize=20)
    if group == "Diamond AMR":
        ax.text(
            x_boundaries.iloc[-1]/2.0, -2.5, r'$\underline{\textbf{\textsf{\Huge{Diamond\ AMR}}}}$', 
            rotation=0, ha="center", color="#3dac87", usetex=True)
    else:
        ax.text(
            x_boundaries.iloc[-1]/2.0, y_boundaries.iloc[-1] + 3.5, group, 
            rotation=0, ha="center", fontsize=20, color="#ffaa5b", weight="bold", style="italic")

    fig.savefig("{a}_nonclustered_amr_switch_distribution.png".format(a = group.split(' ')[0].lower()))
