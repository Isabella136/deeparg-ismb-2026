from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.axes
import seaborn as sn
import pandas as pd
import numpy as np
import matplotlib
import sys

TABLE_LOCATION = "differential_distribution_output"
FEATURE_DATA = "../../database/v2_feature_data.csv"
HIT_COUNT = "../samples/deeparg_hit_count.tsv"
SEQ_COUNT = "../samples/sequence_count.tsv"

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

# Get amr diff count from Deep to Diamond
amr_df = pd.read_csv(f"{TABLE_LOCATION}/switch_pair_abundance_amr.tsv", sep='\t', header=0)
diff_count = amr_df[[
    "sample", "alignment identity", "model", "diamond best-hit label",
    "deeparg hit label", "switch pair count"]].drop_duplicates()
diff_count = diff_count.groupby(by=["sample", "alignment identity", "model"])[["switch pair count"]].sum()

# Get alignment label count 
label_df = pd.read_csv(f"{TABLE_LOCATION}/label_counts.tsv", sep='\t', header=0).rename(columns={"bitscore" : "count"})

# Populate hit_count_amr_df for bar graphs
seq_count_amr_df = pd.read_csv(SEQ_COUNT, sep="\t", header=0)
hit_count_amr_df = pd.read_csv(HIT_COUNT, sep="\t", header=0, index_col=0)

# This column is for top bar graphs
hit_count_amr_df.insert(5, "deeparg seq percentage", 0)
hit_count_amr_df["deeparg seq percentage"] = hit_count_amr_df.apply(
    lambda x: 1.0 - float(x["deeparg hit count"]) / float(x["diamond hit count"]), axis=1)

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
plt.figure(figsize=(40, 15))
axes: list[matplotlib.axes.Axes] = [
    plt.axes((0.05,0.1,0.46,0.4)), plt.axes((0.52,0.1,0.46,0.4)),
    plt.axes((0.05,0.55,0.46,0.4)), plt.axes((0.52,0.55,0.46,0.4))]

# Figure 3A
for model, ax in zip(["SS", "LS"], axes[2:]):
    filtered_hit_count_amr_df = hit_count_amr_df.loc[hit_count_amr_df["model"] == model]
    deeparg = sn.barplot(
        x=range(39),
        y=np.insert(
            arr=filtered_hit_count_amr_df["deeparg seq percentage"].values,
            obj=slice(3,30,3),
            values=np.full(9, 0)), 
        color=cb_palette[1],
        ax=ax)
    
    ax.set_xticks(
        ticks=list(range(39)), 
        labels=np.full(39, ""))
    ax.set_ylim(top=0.6)
    if model == "SS":
        ax.set_yticks(
            ticks=[0,0.1,0.2,0.3,0.4,0.5,0.6], 
            labels=["0.0","10.0","20.0","30.0","40.0","50.0", "60.0"], fontsize=30)
        plt.gcf().text(0.008, 0.75, f"Aligned queries (%)", rotation=90, fontsize=35, va='center')
        plt.gcf().text(0.01, 0.97, "A", fontsize=35, weight='bold', va='center')
    else:
        ax.set_yticks(ticks=[0,0.1,0.2,0.3,0.4,0.5], labels="")
    ax.grid(visible=True, which="major", axis='y', color="0.75")
    ax.tick_params(length=0, pad=8)
    ax.set_title(f"DeepARG-{model}", fontsize=35, loc="left")
    
# Figure 3B
for model, ax in zip(["SS", "LS"], axes[:2]):
    filtered_hit_count_amr_df = hit_count_amr_df.loc[hit_count_amr_df["model"] == model]
    different = sn.barplot(
        x=range(39),
        y=np.insert(
            arr=filtered_hit_count_amr_df["diff percentage"].values,
            obj=slice(3,30,3),
            values=np.full(9, 0)), 
        color=cb_palette[0],
        ax=ax)

    ax.set_xticks(
        ticks=list(range(39)), 
        labels=np.insert(
            arr=filtered_hit_count_amr_df["alignment identity"].values.astype('<U2'),
            obj=slice(3,30,3),
            values=np.full(9, "")),
        rotation_mode="anchor", ha='center', fontsize=25)
    for loc, sample in enumerate(filtered_hit_count_amr_df["sample id"].drop_duplicates().values):
        print(sample)
        ax.text(
            loc*4+1, -0.018, f"#{loc+1}", ha="center", va="center", fontsize=30)
    ax.set_ylim(top=0.15)
    if model == "SS":
        ax.set_yticks(
            ticks=[0,0.025,0.05,0.075,0.1,0.125,0.15], 
            labels=[0.0,2.5,5.0,7.5,10.0,12.5,15.0], fontsize=30)
        plt.gcf().text(0.008, 0.28, f"DeepARG predictions (%)", rotation=90, fontsize=35, va='center')
        plt.gcf().text(0.01, 0.52, "B", fontsize=35, weight='bold', va='center')
    else:
        ax.set_yticks(ticks=[0,0.025,0.05,0.075,0.1,0.125,0.15], labels="")
    ax.grid(visible=True, which="major", axis='y', color="0.75")
    ax.tick_params(length=0, pad=8)
    
plt.gcf().text(0.515, 0.009, "Run (upper value is alignment identity, lower value is sample)", ha="center", fontsize=35)
plt.savefig("barplot.png")

# Rename AMR columns based on abbrev
for row in range(amr_df.shape[0]):
    amr_df.at[row, "diamond best-hit label"] = amr_abbrev[amr_df.at[row, "diamond best-hit label"]]
    amr_df.at[row, "deeparg hit label"] = amr_abbrev[amr_df.at[row, "deeparg hit label"]]

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
        [["Sample ID", "Alignment Identity", "Model", "Query", "Most Frequent Class", "DeepARG Class"]]
        .drop_duplicates())
    
    # For the heatmaps, we will make X-axis be Diamond, and y-axis be DeepARG.
    # For the purpose of space, we will aggregate all samples, meaning that we
    # only need 2x3 categories per pair.
    heatmap_x_axis = pd.MultiIndex.from_product(
        iterables=[
            pair_label_df["Most Frequent Class"].drop_duplicates().sort_values().to_list(),
            pair_label_df["Model"].drop_duplicates().sort_values().to_list()],
        names=["Class", "Model"])
    heatmap_y_axis = pd.MultiIndex.from_product(
        iterables=[
            pair_label_df["DeepARG Class"].drop_duplicates().sort_values().to_list(),
            pair_label_df["Alignment Identity"].drop_duplicates().sort_values().to_list()],
        names=["Class", "Identity"])
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
        [["Sample ID", "Alignment Identity", "Model", "Most Frequent Class", "DeepARG Class", "Query"]]
        .drop_duplicates()
        .groupby(by=[
            "Alignment Identity", "Model", "Most Frequent Class", "DeepARG Class"])
        [["Query"]].count())
    
    if not share:
        amr_count_df = (label_df[["Sample ID", "Alignment Identity", "Model", "Query", "amr", "count"]]
            .groupby(["Sample ID", "Alignment Identity", "Model", "Query", "amr"])
            .sum())

        exp_df = pd.DataFrame(
            data=np.full(
                shape=(heatmap_x_axis.shape[0], heatmap_y_axis.shape[0]),
                fill_value=np.nan),
            index=heatmap_x_axis,
            columns=heatmap_y_axis)

        x_is_best_y_in_alignment = label_df.loc[pd.merge(
            left = label_df[[
                "Alignment Identity", "Model", "Most Frequent Class", "amr"]].rename(
                columns={
                    "amr": "DeepARG Class"}),
            right = pair_label_df[[
                "Alignment Identity", "Model", "Most Frequent Class", "DeepARG Class"]].drop_duplicates(),
            how = 'left',
            indicator='exist')["exist"] == 'both']
        x_is_best_has_y_alignment_abundance = (x_is_best_y_in_alignment
            [["Sample ID", "Alignment Identity", "Model", "Most Frequent Class", "Query", "amr"]]
            .drop_duplicates()
            .groupby(by=["Alignment Identity", "Model", "Most Frequent Class", "amr"])[["Query"]]
            .count())
            
        y_alignments_counts = (x_is_best_y_in_alignment
            [["Sample ID", "Alignment Identity", "Model", 
            "Most Frequent Class", "Query", "amr"]]
            .drop_duplicates())
        y_alignments_counts["count"] = y_alignments_counts.apply(
            lambda x: amr_count_df.at[
                (x["Sample ID"], x["Alignment Identity"], x["Model"], x["Query"], x["amr"]), 
                "count"], axis=1)
        y_alignments_counts.set_index(
            ["Sample ID", "Alignment Identity", "Model", 
            "Most Frequent Class", "Query", "amr"],
            inplace=True)
        
        all_alignment_counts = (label_df
            [["Sample ID", "Alignment Identity", "Model", "Query", "amr"]]
            .drop_duplicates())
        all_alignment_counts["count"] = all_alignment_counts.apply(
            lambda x: amr_count_df.at[
                (x["Sample ID"], x["Alignment Identity"], x["Model"], x["Query"], x["amr"]), 
                "count"], axis=1)
        all_alignment_counts = (all_alignment_counts
            .groupby(by=["Sample ID", "Alignment Identity", "Model", "Query"])[["count"]]
            .sum())
        
        all_alignment_counts = y_alignments_counts.reset_index().apply(
            lambda x: 
                pd.Series(data={
                    "Sample ID": x["Sample ID"], 
                    "Alignment Identity": x["Alignment Identity"], 
                    "Model": x["Model"], 
                    "Most Frequent Class": x["Most Frequent Class"], 
                    "Query": x["Query"], 
                    "amr": x["amr"],
                    "count" : all_alignment_counts.loc[
                        (x["Sample ID"], x["Alignment Identity"], x["Model"], x["Query"])].at["count"]}), axis=1).set_index(
                            ["Sample ID", "Alignment Identity", "Model", 
                            "Most Frequent Class", "Query", "amr"])
        
        probabilities = y_alignments_counts.div(all_alignment_counts)
        probabilities = (probabilities
            .groupby(by=[
                "Alignment Identity", "Model", "Most Frequent Class", "amr"])["count"]
            .mean())
        
        for row in pair_df.iterrows():
            # For y-axis
            deeparg_class = row[0][3]
            identity = row[0][0]

            # For x-axis
            diamond_class = row[0][2]
            diamond_extra_label = row[0][2]
            model = row[0][1]
            
            # For switch_df
            switch = row[1]["Query"]

            # For exp_db
            x_switch_abundance = x_is_best_has_y_alignment_abundance["Query"].at[
                (identity, model, diamond_class, deeparg_class)]
            
            deeparg_probability = probabilities.at[
                (identity, model, diamond_class, deeparg_class)]

            # Insert in dfs
            switch_df.at[(diamond_class, model), (deeparg_class, identity)] = float(switch)
            exp_df.at[(diamond_class, model), (deeparg_class, identity)] = (
                float(deeparg_probability) * float(x_switch_abundance))
            mask_df.at[(diamond_class, model), (deeparg_class, identity)] = False

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
            .groupby(by=[
                "Alignment Identity", "Model", "Most Frequent Class", "DeepARG Class"])
            [["Query"]].count())
        
        y_alignments = (label_df
            .loc[(label_df["Is DeepARG Class"] &
                ~label_df["Is Most Frequent Class"])]
            [["Sample ID", "Alignment Identity", "Model", "Most Frequent Class", 
            "DeepARG Class", "Query", "clstr", "dom", "super"]]
            .drop_duplicates()
            .set_index(["Sample ID", "Alignment Identity", "Model", "Most Frequent Class", 
                "DeepARG Class", "Query"], append=True))
        x_alignments = (label_df
            .loc[(~label_df["Is DeepARG Class"] &
                label_df["Is Most Frequent Class"])]
            [["Sample ID", "Alignment Identity", "Model", "Most Frequent Class", 
            "DeepARG Class", "Query", "clstr", "dom", "super"]]
            .drop_duplicates()
            .set_index(["Sample ID", "Alignment Identity", "Model", "Most Frequent Class", 
                "DeepARG Class", "Query"], append=True))  
        
        for row in pair_label_df.iterrows():
            # For y-axis
            deeparg = row[1]["DeepARG Class"]
            identity = row[1]["Alignment Identity"]

            # For x-axis
            diamond = row[1]["Most Frequent Class"]
            model = row[1]["Model"]

            #for alignment crossing
            sample = row[1]["Sample ID"]
            query = row[1]["Query"]
            share = 0

            x_query_alignment = x_alignments.loc[
                (slice(None), sample, identity, model, diamond, deeparg, query)][["clstr", "dom", "super"]]
            y_query_alignment = y_alignments.loc[
                (slice(None), sample, identity, model, diamond, deeparg, query)][["clstr", "dom", "super"]]
            if True in (x_query_alignment["clstr"].isin(y_query_alignment["clstr"]) | 
                        x_query_alignment["dom"].isin(y_query_alignment["dom"])| 
                        x_query_alignment["super"].isin(y_query_alignment["super"])).to_list():
                share = 1
            mask_df.at[(diamond, model), (deeparg, identity)] = False
            
            # For switch_df
            switch = pair_df.at[(identity, model, diamond, deeparg), "Query"]

            switch_df.at[(diamond, model), (deeparg, identity)] = float(switch)
            share_df.at[(diamond, model), (deeparg, identity)] += float(share)
        
        heatmap_df = share_df.div(switch_df, fill_value=np.nan)
    
    # And now, heatmaps!
    fig = plt.figure(figsize=[18, 18])
    ax_left = fig.add_axes((0.1, 0.08, 0.75, 0.84))
    cbar = fig.add_axes((0.93, 0.2, 0.02, 0.6))
    sn.heatmap(
        data = heatmap_df.transpose(), 
        mask=mask_df.transpose(), 
        cmap=custom_cmap, 
        center=0, 
        ax=ax_left, 
        cbar_ax=cbar,
        xticklabels=heatmap_df.index.get_level_values("Class"),
        yticklabels=heatmap_df.columns.get_level_values("Class"))

    cbar.tick_params(labelsize=20)

    # Label x-axis
    ax_left.set_xlabel("Diamond Class label", fontsize=35)

    # Label y-axis
    ax_left.set_ylabel(
        ylabel="DeepARG class label",
        fontsize=35)

    # Label x-ticks
    ax_left.set_xticks(
        ticks=np.array(range(1, len(heatmap_df.index.get_level_values("Class")), 2)),
        labels=heatmap_df.index.get_level_values("Class").drop_duplicates(),
        rotation_mode='anchor',
        rotation=45,
        ha='right',
        va='top',
        fontsize=20)

    # Label y-ticks
    ax_left.set_yticks(
        ticks=np.array(range(1, len(heatmap_df.columns.get_level_values("Class")), 3)) + 0.5,
        labels=heatmap_df.columns.get_level_values("Class").drop_duplicates(),
        rotation_mode='anchor',
        fontsize=20)

    # Seperate by group
    y_groups = heatmap_df.columns.to_frame(index=False).groupby(["Class"]).size()
    x_groups = heatmap_df.index.to_frame(index=False).groupby(["Class"]).size()

    y_boundaries = np.cumsum(y_groups)
    x_boundaries = np.cumsum(x_groups)

    for y_loc in y_boundaries[:-1]:
        ax_left.hlines(
            y_loc, 
            xmin=ax_left.get_xlim()[0], 
            xmax=ax_left.get_xlim()[1], 
            colors="white", 
            linewidth=5)
    for x_loc in x_boundaries[:-1]:
        ax_left.vlines(
            x_loc, 
            ymin=ax_left.get_ylim()[0], 
            ymax=ax_left.get_ylim()[1], 
            colors="white", 
            linewidth=5)
        
    # Add heatmap cell key:
    fig.text(
        x=0.92, y=0.96,
        s="Heatmap Cell:",
        fontsize=25,
        ha='center',
        va='bottom')
    matplotlib.rc('text', usetex=True)
    fig.text(
        x=0.92, y=0.95,
        s=r'''\begin{tabular}{ c | c | c |} & LS & SS \\ \hline 30 & & \\ \hline 50 & & \\ \hline 80 & & \end{tabular}''',
        fontsize=25,
        ha='center',
        va='top')
    matplotlib.rc('text', usetex=False)
        
    if not share:
        fig.savefig(f"most_freq_amr_switch_relative_to_alignment.png")
    else:
        fig.savefig(f"most_freq_share_label_with_deeparg.png")
        