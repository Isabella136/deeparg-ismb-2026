from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib.axes
import seaborn as sn
import pandas as pd
import numpy as np
import matplotlib
import sys

FEATURE_DATA = "../../database/v2_feature_data.csv"
HIT_COUNT = "../samples/deeparg_hit_count.tsv"
SEQ_COUNT = "../samples/sequence_count.tsv"

matplotlib.rcParams["mathtext.fontset"] = "dejavusans"
matplotlib.rcParams["font.family"] = "DejaVu Sans"

  
abbrev_file = open("amr_abbrev.csv", "r")
amr_abbrev = dict()
for line in abbrev_file.readlines():
    row = line.strip().split(",")
    amr_abbrev.update({row[0]: row[1]})
abbrev_file.close()

# Get alignment label count 
label_df = pd.read_csv("label_counts.tsv", sep='\t', header=0).rename(columns={"bitscore" : "count"})

# Populate hit_count_amr_df for bar graphs
seq_count_amr_df = pd.read_csv(SEQ_COUNT, sep="\t", header=0)
hit_count_amr_df = pd.read_csv(HIT_COUNT, sep="\t", header=0, index_col=0)

# This column is for top bar graphs
hit_count_amr_df.insert(5, "deeparg seq percentage", 0)
hit_count_amr_df["deeparg seq percentage"] = hit_count_amr_df.apply(
    lambda x: float(x["deeparg hit count"]) / float(x["diamond hit count"]), axis=1)

# This column is to later calculate values for bottom bar graphs
hit_count_amr_df.insert(6, "diff count", 0)
hit_count_amr_df["diff count"] = hit_count_amr_df.apply(
        lambda x: label_df
        .loc[(
            (label_df["Sample ID"] == x["sample id"]) &
            (label_df["Alignment Identity"] == x["alignment identity"]) &
            (label_df["Model"] == x["model"]) &
            (label_df["Is DeepARG Class"] != label_df["Is Most Frequent Class"]))][
                ["Sample ID", "Alignment Identity", "Model", "Query"]]
        .drop_duplicates()
        .groupby(by=["Sample ID", "Alignment Identity", "Model"])["Query"]
        .count()
        .loc[(x["sample id"], x["alignment identity"], x["model"])] if len(label_df
            .loc[((label_df["Sample ID"] == x["sample id"]) &
                (label_df["Alignment Identity"] == x["alignment identity"]) &
                (label_df["Model"] == x["model"]) &
                (label_df["Is DeepARG Class"] != label_df["Is Most Frequent Class"]))].index) > 0 else 0, axis=1)

# This column is for the bottom bar graphs
hit_count_amr_df = hit_count_amr_df.sort_values(by=["model", "sample id", "alignment identity"])
hit_count_amr_df.insert(7, "diff percentage", 0)
hit_count_amr_df["diff percentage"] = hit_count_amr_df.apply(
    lambda x: float(x["diff count"]) / float(hit_count_amr_df.loc[
        (hit_count_amr_df["sample id"] == x["sample id"]) & 
        (hit_count_amr_df["alignment identity"] == x["alignment identity"]) & 
        (hit_count_amr_df["model"] == x["model"])].iat[0, 3]), axis=1)

cb_palette = sn.color_palette("colorblind")

# Make figure 3
plt.figure(figsize=(40, 8))
axes: list[matplotlib.axes.Axes] = [
    plt.axes((0.05,0.12,0.44,0.8)), plt.axes((0.54,0.12,0.44,0.8))]

# Figure 3A
deeparg = sn.barplot(
    x=range(6),
    y=hit_count_amr_df.groupby(
        ["model", "alignment identity"])["deeparg seq percentage"].mean() ,
    color=cb_palette[1],
    ax=axes[0])
axes[0].set_xticks(
    ticks=range(6), 
    labels=["LS-30%", "LS-50%", "LS-80%", "SS-30%", "SS-50%", "SS-80%"],
    rotation_mode="anchor", ha='center', fontsize=25)
axes[0].set_yticks(
    ticks=[0,0.2,0.4,0.6,0.8,1], 
    labels=["0.0","20.0","40.0","60.0","80.0","100.0"], fontsize=30)
axes[0].set_ylabel(f"Average Percentage", rotation=90, fontsize=35, va='center')
plt.gcf().text(0.01, 0.96, "A", fontsize=35, weight='bold')
axes[0].grid(visible=True, which="major", axis='y', color="0.75")
axes[0].tick_params(length=0, pad=8)
    
# Figure 3B
different = sn.barplot(
    x=range(6),
    y=hit_count_amr_df.groupby(
        ["model", "alignment identity"])["diff percentage"].mean(),
    color=cb_palette[0],
    ax=axes[1])

axes[1].set_xticks(
    ticks=range(6), 
    labels=["LS-30%", "LS-50%", "LS-80%", "SS-30%", "SS-50%", "SS-80%"],
    rotation_mode="anchor", ha='center', fontsize=25)
axes[1].set_ylim(top=0.06)
axes[1].set_yticks(
    ticks=[0,0.01,0.02,0.03, 0.04, 0.05,0.06], 
    labels=[0.0,1.0,2.0,3.0,4.0,5.0,6.0], fontsize=30)
axes[1].set_ylabel(f"Average Percentage", rotation=90, fontsize=35)
plt.gcf().text(0.5, 0.96, "B", fontsize=35, weight='bold', va='center')
axes[1].grid(visible=True, which="major", axis='y', color="0.75")
axes[1].tick_params(length=0, pad=8)
    
plt.gcf().text(0.515, 0.009, "Model-Alignment Identity Condition", ha="center", fontsize=35)
plt.savefig("barplot.png")

# Only mapping a specific identity and model
model = sys.argv[1]
ident = int(sys.argv[2])
label_df = (label_df.loc[
    (label_df["Model"] == model) &
    (label_df["Alignment Identity"] == ident)].reset_index(drop=True))

# Rename AMR columns based on abbrev
label_df["amr"] = label_df["amr"].apply(
    lambda x: amr_abbrev[x])
label_df["Diamond Class"] = label_df["Diamond Class"].apply(
    lambda x: amr_abbrev[x])
label_df["Most Frequent Class"] = label_df["Most Frequent Class"].apply(
    lambda x: amr_abbrev[x])
label_df["DeepARG Class"] = label_df["DeepARG Class"].apply(
    lambda x: amr_abbrev[x])

# Create a custom color palette
custom_colors = [
'#D20A2E', # first color
'#d9d9d9', # middle color
'#0F52BA' # last color
]
custom_cmap = LinearSegmentedColormap.from_list("custom_gradient", custom_colors)
custom_cmap.set_bad()

for share in [True, False]:
    # Must get the queries of interest:
    pair_label_df = (label_df
        .loc[label_df["Is Most Frequent Class"] != label_df["Is DeepARG Class"]]
        [["Sample ID", "Query", "Most Frequent Class", "DeepARG Class"]]
        .drop_duplicates())
    
    # For the heatmaps, we will make X-axis be Diamond and y-axis be DeepARG.
    heatmap_x_axis = pd.Index(
        pair_label_df["Most Frequent Class"].drop_duplicates().sort_values().to_list())
    heatmap_y_axis = pd.Index(
        pair_label_df["DeepARG Class"].drop_duplicates().sort_values().to_list())
    mask_df = pd.DataFrame(
        data=np.full(
            shape=(heatmap_x_axis.shape[0], heatmap_y_axis.shape[0]),
            fill_value=True),
        index=heatmap_x_axis,
        columns=heatmap_y_axis)
    switch_df = pd.DataFrame(
        data=np.full(
            shape=(heatmap_x_axis.shape[0], heatmap_y_axis.shape[0]),
            fill_value=np.nan),
        index=heatmap_x_axis,
        columns=heatmap_y_axis)
    
    # Aggregate pair switch
    pair_df = (pair_label_df
        .loc[label_df["Is Most Frequent Class"] != label_df["Is DeepARG Class"]]
        [["Sample ID", "Most Frequent Class", "DeepARG Class", "Query"]]
        .drop_duplicates()
        .groupby(by=["Most Frequent Class", "DeepARG Class"])[["Query"]].count())
    
    if not share:
        amr_count_df = (label_df[["Sample ID", "Query", "amr", "count"]]
            .groupby(["Sample ID", "Query", "amr"])
            .sum())

        exp_df = pd.DataFrame(
            data=np.full(
                shape=(heatmap_x_axis.shape[0], heatmap_y_axis.shape[0]),
                fill_value=np.nan),
            index=heatmap_x_axis,
            columns=heatmap_y_axis)

        x_is_best_y_in_alignment = label_df.loc[pd.merge(
            left = label_df[["Most Frequent Class", "amr"]].rename(
                columns={"amr": "DeepARG Class"}),
            right = pair_label_df[[
                "Most Frequent Class", "DeepARG Class"]].drop_duplicates(),
            how = 'left',
            indicator='exist')["exist"] == 'both']
            
        y_alignments_counts = (x_is_best_y_in_alignment
            [["Sample ID", "Most Frequent Class", "Query", "amr"]]
            .drop_duplicates())
        y_alignments_counts["count"] = y_alignments_counts.apply(
            lambda x: amr_count_df.at[
                (x["Sample ID"], x["Query"], x["amr"]), "count"], axis=1)
        y_alignments_counts.set_index(
            ["Sample ID", "Most Frequent Class", "Query", "amr"],
            inplace=True)
        
        all_alignment_counts = (label_df[["Sample ID", "Query", "amr"]]
            .drop_duplicates())
        all_alignment_counts["count"] = all_alignment_counts.apply(
            lambda x: amr_count_df.at[
                (x["Sample ID"], x["Query"], x["amr"]), "count"], axis=1)
        all_alignment_counts = (all_alignment_counts
            .groupby(by=["Sample ID", "Query"])[["count"]]
            .sum())
        
        all_alignment_counts = y_alignments_counts.reset_index().apply(
            lambda x: 
                pd.Series(data={
                    "Sample ID": x["Sample ID"],
                    "Most Frequent Class": x["Most Frequent Class"], 
                    "Query": x["Query"], 
                    "amr": x["amr"],
                    "count" : all_alignment_counts.loc[
                        (x["Sample ID"], x["Query"])].at["count"]}), axis=1).set_index(
                            ["Sample ID", "Most Frequent Class", "Query", "amr"])
        
        probabilities = y_alignments_counts.div(all_alignment_counts)
        probabilities = (probabilities
            .groupby(by=["Most Frequent Class", "amr"])["count"]
            .sum())
        
        for row in pair_df.iterrows():
            # For x- and y-axis
            deeparg_class = row[0][1]
            diamond_class = row[0][0]
            
            # For switch_df
            switch = row[1]["Query"]

            # For exp_db            
            expected = probabilities.at[(diamond_class, deeparg_class)]

            # Insert in dfs
            switch_df.at[diamond_class, deeparg_class] = float(switch)
            exp_df.at[diamond_class, deeparg_class] = float(expected)
            mask_df.at[diamond_class, deeparg_class] = False

        # Now let's get the heatmap values
        heatmap_df = (switch_df
            .add(1, fill_value=np.nan)
            .div(exp_df.add(1), fill_value=np.nan)
            .map(lambda x: np.log(x), na_action='ignore')
            .sort_index())

    else:        
        share_df = pd.DataFrame(
            data=np.full(
                shape=(heatmap_x_axis.shape[0], heatmap_y_axis.shape[0]),
                fill_value=0),
            index=heatmap_x_axis,
            columns=heatmap_y_axis)

        # Aggregate pair switch
        pair_df = (pair_label_df
            .groupby(by=["Most Frequent Class", "DeepARG Class"])
            [["Query"]].count())
        
        y_alignments = (label_df
            .loc[(label_df["Is DeepARG Class"] &
                ~label_df["Is Most Frequent Class"])]
            [["Sample ID", "Most Frequent Class", "DeepARG Class", "Query", "clstr", "dom", "super"]]
            .drop_duplicates()
            .set_index(
                ["Sample ID", "Most Frequent Class", "DeepARG Class", "Query"], 
                append=True))
        x_alignments = (label_df
            .loc[(~label_df["Is DeepARG Class"] &
                label_df["Is Most Frequent Class"])]
            [["Sample ID", "Most Frequent Class", "DeepARG Class", "Query", "clstr", "dom", "super"]]
            .drop_duplicates()
            .set_index(
                ["Sample ID", "Most Frequent Class", "DeepARG Class", "Query"],
                append=True))
        
        for row in pair_label_df.iterrows():
            # For x- and y-axis
            deeparg = row[1]["DeepARG Class"]
            diamond = row[1]["Most Frequent Class"]

            # Check for annotation sharing  
            sample = row[1]["Sample ID"]
            query = row[1]["Query"]
            annot_share = 0
            x_query_alignment = x_alignments.loc[
                (slice(None), sample, diamond, deeparg, query)][["clstr", "dom", "super"]]
            y_query_alignment = y_alignments.loc[
                (slice(None), sample, diamond, deeparg, query)][["clstr", "dom", "super"]]
            if True in (x_query_alignment["clstr"].isin(y_query_alignment["clstr"]) | 
                        x_query_alignment["dom"].isin(y_query_alignment["dom"])| 
                        x_query_alignment["super"].isin(y_query_alignment["super"])).to_list():
                annot_share = 1
            mask_df.at[diamond, deeparg] = False
            
            # For switch_df
            switch = pair_df.at[(diamond, deeparg), "Query"]

            switch_df.at[diamond, deeparg] = float(switch)
            share_df.at[diamond, deeparg] += float(annot_share)
        
        heatmap_df = share_df.div(switch_df, fill_value=np.nan)
    
    # And now, heatmaps!
    fig = plt.figure(figsize=[18, 18])
    ax_left = fig.add_axes((0.12, 0.12, 0.75, 0.84))
    cbar = fig.add_axes((0.92, 0.24, 0.02, 0.6))
    sn.heatmap(
        data = heatmap_df.transpose(), 
        mask=mask_df.transpose(), 
        cmap=custom_cmap, 
        center=0, 
        vmin=0 if share else -6,
        vmax=1 if share else 3,
        ax=ax_left, 
        cbar_ax=cbar)

    cbar.tick_params(labelsize=30)

    # Label x-axis
    ax_left.set_xlabel("Diamond Class label", fontsize=35)

    # Label y-axis
    ax_left.set_ylabel(
        ylabel="DeepARG class label",
        fontsize=35)

    # Label x-ticks
    ax_left.set_xticks(
        ticks=ax_left.get_xticks(),
        labels=heatmap_df.index.drop_duplicates(),
        rotation_mode='anchor',
        rotation=45,
        ha='right',
        va='top',
        fontsize=30)

    # Label y-ticks
    ax_left.set_yticks(
        ticks=ax_left.get_yticks(),
        labels=heatmap_df.columns.drop_duplicates(),
        rotation_mode='anchor',
        rotation=0,
        ha='right',
        va='center',
        fontsize=30)

    # # Seperate by group
    # y_groups = heatmap_df.columns.to_frame(index=False).groupby(["Class"]).size()
    # x_groups = heatmap_df.index.to_frame(index=False).groupby(["Class"]).size()

    y_boundaries = range(heatmap_df.columns.size)
    x_boundaries = range(heatmap_df.index.size)

    for y_loc in y_boundaries[1:]:
        ax_left.hlines(
            y_loc, 
            xmin=ax_left.get_xlim()[0], 
            xmax=ax_left.get_xlim()[1], 
            colors="white", 
            linewidth=5)
    for x_loc in x_boundaries[1:]:
        ax_left.vlines(
            x_loc, 
            ymin=ax_left.get_ylim()[0], 
            ymax=ax_left.get_ylim()[1], 
            colors="white", 
            linewidth=5)
        
    # Add heatmap cell key:
    # fig.text(
    #     x=0.92, y=0.96,
    #     s="Heatmap Cell:",
    #     fontsize=25,
    #     ha='center',
    #     va='bottom')
    # matplotlib.rc('text', usetex=True)
    # fig.text(
    #     x=0.92, y=0.95,
    #     s=r'''\begin{tabular}{ c | c | c |} & LS & SS \\ \hline 30 & & \\ \hline 50 & & \\ \hline 80 & & \end{tabular}''',
    #     fontsize=25,
    #     ha='center',
    #     va='top')
    # matplotlib.rc('text', usetex=False)
        
    if not share:
        fig.savefig(f"most_freq_amr_switch_relative_to_alignment_{model}_{ident}.png")
    else:
        fig.savefig(f"most_freq_share_label_with_deeparg_{model}_{ident}.png")
        