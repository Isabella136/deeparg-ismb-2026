from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.axes
import seaborn as sn
import pandas as pd
import numpy as np
import matplotlib

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
df = pd.read_csv(f"{TABLE_LOCATION}/switch_pair_abundance_amr.tsv", sep='\t', header=0)
diff_count = df[[
    "sample", "alignment identity", "model", "diamond best-hit label",
    "deeparg hit label", "switch pair count"]].drop_duplicates()
diff_count = diff_count.groupby(by=["sample", "alignment identity", "model"])[["switch pair count"]].sum()

# Populate hit_count_df for bar graphs
seq_count_df = pd.read_csv(SEQ_COUNT, sep="\t", header=0)
hit_count_df = pd.read_csv(HIT_COUNT, sep="\t", header=0, index_col=0)

bar_graph = False

if bar_graph:
    # Those two columns are for the top bar graphs
    hit_count_df.insert(5, "deeparg seq percentage", 0)
    hit_count_df["deeparg seq percentage"] = hit_count_df.apply(
        lambda x: float(x["deeparg hit count"]) / float(seq_count_df.loc[
            seq_count_df["sample id"] == x["sample id"]].iat[
                0, 1 if x["model"] == 'SS' else 2]), axis=1)
    hit_count_df.insert(6, "diamond seq percentage", 0)
    hit_count_df["diamond seq percentage"] = hit_count_df.apply(
        lambda x: float(x["diamond hit count"]) / float(seq_count_df.loc[
            seq_count_df["sample id"] == x["sample id"]].iat[
                0, 1 if x["model"] == 'SS' else 2]), axis=1)

    # This column is to later calculate values for bottom bar graphs
    hit_count_df.insert(7, "diff count", 0)
    hit_count_df["diff count"] = hit_count_df.apply(
        lambda x: get_hit_count(diff_count, x), axis=1)

    # Those two columns are for the bottom bar graphs
    hit_count_df = hit_count_df.sort_values(by=["model", "sample id", "alignment identity"])
    hit_count_df.insert(8, "hit percentage", 0)
    hit_count_df["hit percentage"] = hit_count_df.apply(
        lambda x: float(x["deeparg hit count"]) / float(hit_count_df.loc[
            (hit_count_df["sample id"] == x["sample id"]) & 
            (hit_count_df["alignment identity"] == 30) & 
            (hit_count_df["model"] == x["model"])].iat[0, 3]), axis=1)
    hit_count_df.insert(9, "diff percentage", 0)
    hit_count_df["diff percentage"] = hit_count_df.apply(
        lambda x: float(x["diff count"]) / float(hit_count_df.loc[
            (hit_count_df["sample id"] == x["sample id"]) & 
            (hit_count_df["alignment identity"] == 30) & 
            (hit_count_df["model"] == x["model"])].iat[0, 3]), axis=1)

    cb_palette = sn.color_palette("colorblind")

    # Make bar graph
    plt.figure(figsize=(40, 22))
    axes: list[matplotlib.axes.Axes] = [
        plt.axes((0.04,0.06,0.46,0.43)), plt.axes((0.52,0.06,0.46,0.43)),
        plt.axes((0.04,0.53,0.46,0.43)), plt.axes((0.52,0.53,0.46,0.43))]

    # We are making the overall diamond vs deeparg results
    for model, ax in zip(["SS", "LS"], axes[2:]):
        filtered_hit_count_df = hit_count_df.loc[hit_count_df["model"] == model]
        diamond = sn.barplot(
            x=range(39),
            y=np.insert(
                arr=filtered_hit_count_df["diamond seq percentage"].values,
                obj=slice(3,30,3),
                values=np.full(9, 0)),
            color=cb_palette[2],
            ax=ax)
        deeparg = sn.barplot(
            x=range(39),
            y=np.insert(
                arr=filtered_hit_count_df["deeparg seq percentage"].values,
                obj=slice(3,30,3),
                values=np.full(9, 0)), 
            color=cb_palette[3],
            ax=ax)
        
        if model == "LS":
            top_bar = mpatches.Patch(color=cb_palette[2], label='Not kept by DeepARG')
            bottom_bar = mpatches.Patch(color=cb_palette[3], label='Kept by DeepARG')
            ax.legend(handles=[top_bar, bottom_bar], fontsize=24)
        ax.set_xticks(
            ticks=list(range(39)), 
            labels=np.full(39, ""))
        ax.set_ylim(top=0.05)
        if model == "SS":
            ax.set_yticks(
                ticks=[0,0.01,0.02,0.03,0.04,0.05], 
                labels=["0","1","2","3","4","5"], fontsize=24)
            plt.gcf().text(0.015, 0.745, f"Percentage of input sequences", rotation=90, fontsize=30, va='center')
            plt.gcf().text(0.015, 0.98, "A", fontsize=40, weight='bold', va='center')
        else:
            ax.set_yticks(ticks=[0,0.01,0.02,0.03,0.04,0.05], labels="")
        ax.grid(visible=True, which="major", axis='y', color="0.75")
        ax.tick_params(length=0, pad=8)
        ax.set_title(f"Percentage of {"Reads" if model=="SS" else "ORFs"} Aligned By Diamond and Kept by DeepARG", fontsize=36, loc="left")

    # We are now making the switch vs no switch bar graph
    for model, ax in zip(["SS", "LS"], axes[:2]):
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

        if model == "LS":
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
                loc*4+1, -0.07, f"#{loc+1}", ha="center", va="center", fontsize=28)
        ax.set_ylim(top=1.0)
        if model == "SS":
            ax.set_yticks(
                ticks=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1], 
                labels=list(range(0, 110, 10)), fontsize=24)
            plt.gcf().text(0.015, 0.275, f"Percentage of total DeepARG 30% hits", rotation=90, fontsize=30, va='center')
            plt.gcf().text(0.015, 0.51, "B", fontsize=40, weight='bold', va='center')
        else:
            ax.set_yticks(ticks=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1], labels="")
        ax.grid(visible=True, which="major", axis='y', color="0.75")
        ax.tick_params(length=0, pad=8)
        ax.set_title(f"DeepARG-{model} Similarity to Diamond Best Hit", fontsize=36, loc="left")

    plt.gcf().text(0.51, 0.009, "Run (upper value is alignment identity, lower value is sample)", ha="center", fontsize=30)
    plt.savefig("barplot.png")

# Rename AMR columns based on abbrev
for row in range(df.shape[0]):
    df.at[row, "diamond best-hit label"] = amr_abbrev[df.at[row, "diamond best-hit label"]]
    df.at[row, "deeparg hit label"] = amr_abbrev[df.at[row, "deeparg hit label"]]

# For the heatmaps, we will make X-axis be Diamond, and y-axis be DeepARG.
# For the purpose of space, we will aggregate all samples, meaning that we
# only need 2x3 categories per pair.

heatmap_x_axis = pd.MultiIndex.from_product(
    iterables=[
        df["diamond best-hit label"].drop_duplicates().sort_values().to_list(),
        df["model"].drop_duplicates().sort_values().to_list()],
    names=["Class", "Model"])
heatmap_y_axis = pd.MultiIndex.from_product(
    iterables=[
        df["deeparg hit label"].drop_duplicates().sort_values().to_list(),
        df["alignment identity"].drop_duplicates().sort_values().to_list()],
    names=["Class", "Identity"])
switch_df = pd.DataFrame(
    data=np.full(
        shape=(heatmap_x_axis.shape[0], heatmap_y_axis.shape[0]),
        fill_value=0.0),
    index=heatmap_x_axis,
    columns=heatmap_y_axis)
mask_df = pd.DataFrame(
    data=np.full(
        shape=(heatmap_x_axis.shape[0], heatmap_y_axis.shape[0]),
        fill_value=True),
    index=heatmap_x_axis,
    columns=heatmap_y_axis)

# Make a new heatmap for expected number of amr switch pair X, Y given X and
# the DeepARG distribution of Y. We will need the proportion of non-class X
# queries that are DeepARG hits of class Y, and the number of AMR switches
# from class X.

exp_query_y_from_x = pd.DataFrame(
    data=np.full(
        shape=(heatmap_x_axis.shape[0], heatmap_y_axis.shape[0]),
        fill_value=1.0),
    index=heatmap_x_axis,
    columns=heatmap_y_axis)

# Now, we also need a new heatmap for expected number of amr switch pair X, Y
# given Y and the Diamond distribution of Y. We will need the proportion of
# non-class Y queries that are Diamond hits of class X, and the number of AMR
# switches to class Y

exp_query_x_to_y = pd.DataFrame(
    data=np.full(
        shape=(heatmap_x_axis.shape[0], heatmap_y_axis.shape[0]),
        fill_value=1.0),
    index=heatmap_x_axis,
    columns=heatmap_y_axis)

# Finally, we need a new heatmap for expected number of amr switch pair X, Y 
# given X and the database distribution of Y. We will need the proportion of 
# non-class X features that are class Y, and the number of AMR switches from class
# X.

exp_db = pd.DataFrame(
    data=np.full(
        shape=(heatmap_x_axis.shape[0], heatmap_y_axis.shape[0]),
        fill_value=1.0),
    index=heatmap_x_axis,
    columns=heatmap_y_axis)

# Aggregate indiv switch and abundance counts
df_pair = df[[
    "sample", "alignment identity", "model", "diamond best-hit label",
    "diamond best-hit state a abundance", "diamond best-hit state b abundance", 
    "deeparg hit label", "deeparg hit state a abundance",
    "deeparg hit state b abundance", "switch pair count"]].drop_duplicates()
df_pair = (df_pair
    .groupby(by=[
        "alignment identity", "model", "diamond best-hit label", "deeparg hit label"])
    [["diamond best-hit state a abundance", "diamond best-hit state b abundance", 
      "deeparg hit state a abundance", "deeparg hit state b abundance", 
      "switch pair count"]].sum())

# Get class count for DeepARG database by retrieving feature data
feature = pd.read_csv(FEATURE_DATA, header=0)
abbrev_file = open("../../database/amr_abbrev.csv", "r")
amr_abbrev = dict()
for line in abbrev_file.readlines():
    row = line.strip().split(",")
    amr_abbrev.update({row[0]: row[1]})
abbrev_file.close()
class_count_df = feature.drop_duplicates(subset=["amr class", "amr class v2 count"])
class_count_df["amr class"] = class_count_df["amr class"].apply(
    lambda x: amr_abbrev[x])
class_count_df.set_index(keys="amr class", inplace=True)

# Aggregate switch from x counts
df_from_x = (df
    [["sample", "alignment identity", "model", "diamond best-hit label",
        "deeparg hit label", "switch pair count"]]
    .groupby(by=[
        "alignment identity", "model", "diamond best-hit label"])
    [["switch pair count"]].sum())

# Aggregate switch to y counts
df_to_y = (df
    [["sample", "alignment identity", "model", "diamond best-hit label",
        "deeparg hit label", "switch pair count"]]
    .groupby(by=[
        "alignment identity", "model", "deeparg hit label"])
    [["switch pair count"]].sum())

# Aggregate query count
hit_count_agg = (hit_count_df
    [["sample id", "alignment identity", "model", "deeparg hit count"]]
    .groupby(by=["alignment identity", "model"])
    [["deeparg hit count"]].sum())


for row in df_pair.iterrows():
    # For x-axis
    deeparg = row[0][3]
    identity = row[0][0]

    # For y-axis
    diamond = row[0][2]
    model = row[0][1]
    
    # For switch_df
    switch = row[1]["switch pair count"]

    # For exp_query_y_from_x
    y_b_abundance = row[1]["deeparg hit state b abundance"]
    x_b_abundance = row[1]["diamond best-hit state b abundance"]

    # For exp_query_y_from_x and exp_db
    x_switch_abundance = df_from_x.loc[
        (identity, model, diamond), "switch pair count"]
    y_ref_abundance = class_count_df.loc[
        deeparg, "amr class v2 count"]
    not_x_ref_abundance = 12279 - class_count_df.loc[
        diamond, "amr class v2 count"]
    
    # For exp_query_x_to_y
    x_a_abundance = row[1]["diamond best-hit state a abundance"]
    y_a_abundance = row[1]["deeparg hit state a abundance"]
    y_switch_abundance = df_to_y.loc[
        (identity, model, deeparg), "switch pair count"]

    # For both exp_query
    query_count = hit_count_agg.loc[
        (identity, model), "deeparg hit count"]
    
    # Insert in dfs
    switch_df.at[(diamond, model), (deeparg, identity)] = float(switch)
    exp_query_y_from_x.at[(diamond, model), (deeparg, identity)] = (
        float(x_switch_abundance) / float(query_count - x_b_abundance) * float(y_b_abundance))
    exp_query_x_to_y.at[(diamond, model), (deeparg, identity)] = (
        float(y_switch_abundance) / float(query_count - y_a_abundance) * float(x_a_abundance))
    exp_db.at[(diamond, model), (deeparg, identity)] = (
        float(y_ref_abundance) / float(not_x_ref_abundance) * float(x_switch_abundance))
    mask_df.at[(diamond, model), (deeparg, identity)] = False

# Now let's get the heatmap values
df_heatmap_y_from_x = (switch_df
    .add(1)
    .div(exp_query_y_from_x.add(1))
    .map(lambda x: np.log(x))
    .sort_index())
df_heatmap_x_to_y = (switch_df
    .add(1)
    .div(exp_query_x_to_y.add(1))
    .map(lambda x: np.log(x))
    .sort_index())
df_heatmap_db = (switch_df
    .add(1)
    .div(exp_db.add(1))
    .map(lambda x: np.log(x))
    .sort_index())

max_val = max(
    df_heatmap_y_from_x.max().max(), 
    #df_heatmap_x_to_y.max().max(),
    df_heatmap_db.max().max())
min_val = min(
    df_heatmap_y_from_x.min().min(), 
    #df_heatmap_x_to_y.min().min(),
    df_heatmap_db.min().min())

# Create a custom color palette
custom_colors = [
'#D20A2E', # first color
'#d9d9d9', # middle color
'#0F52BA' # last color
]
custom_cmap = LinearSegmentedColormap.from_list("custom_gradient", custom_colors)
custom_cmap.set_bad()

# And now, heatmaps!
fig = plt.figure(figsize=[30, 18])
ax_left = fig.add_axes((0.07, 0.08, 0.4, 0.84))
ax_right = fig.add_axes((0.48, 0.08, 0.4, 0.84))
#ax_left = fig.add_axes((0.1, 0.08, 0.75, 0.84))
cbar = fig.add_axes((0.93, 0.2, 0.02, 0.6))
sn.heatmap(
    data = df_heatmap_y_from_x.transpose(), 
    mask=mask_df.transpose(), 
    cmap=custom_cmap, 
    center=0, 
    vmax=max_val, 
    vmin=min_val,
    ax=ax_left, 
    cbar_ax=cbar,
    xticklabels=heatmap_x_axis.get_level_values("Class"),
    yticklabels=heatmap_y_axis.get_level_values("Class"))

sn.heatmap(
    data = df_heatmap_x_to_y.transpose(), 
    #data = df_heatmap_db.transpose(),
    mask=mask_df.transpose(), 
    cmap=custom_cmap, 
    center=0, 
    vmax=max_val, 
    vmin=min_val,
    ax=ax_right, 
    cbar_ax=cbar,
    xticklabels=heatmap_x_axis.get_level_values("Class"),
    yticklabels=heatmap_y_axis.get_level_values("Class"))

cbar.tick_params(labelsize=20)

#Label x-axis
fig.text(
    x=0.475, y=0.01,
    s="Diamond class label",
    ha="center",
    fontsize=35)
# ax_left.set_xlabel("Diamond class label", fontsize=35)
ax_left.set_xlabel("")
ax_right.set_xlabel("")

# Label y-axis
ax_left.set_ylabel(
    ylabel="DeepARG class label",
    fontsize=35)
ax_right.set_ylabel("")

# Label x-ticks
ax_left.set_xticks(
    ticks=np.array(range(1, 44, 2)),
    labels=heatmap_x_axis.get_level_values("Class").drop_duplicates(),
    rotation_mode='anchor',
    rotation=45,
    ha='right',
    va='top',
    fontsize=20)
ax_right.set_xticks(
    ticks=np.array(range(1, 44, 2)),
    labels=heatmap_x_axis.get_level_values("Class").drop_duplicates(),
    rotation_mode='anchor',
    rotation=45,
    ha='right',
    va='top',
    fontsize=20)

# Label y-ticks
ax_left.set_yticks(
    ticks=np.array(range(1, 69, 3)) + 0.5,
    labels=heatmap_y_axis.get_level_values("Class").drop_duplicates(),
    rotation_mode='anchor',
    fontsize=20)
ax_right.set_yticks(
    ticks=np.array(range(1, 69, 3)) + 0.5,
    labels=np.full(shape=23, fill_value=""))

# Title
ax_left.set_title(
    label="$\\mathbb{E}[s_{X,Y}] = \\frac{|Y.B|}{|Q|-|X.B|}\\times\\sum_{K=1}^{N}s_{X,K}$",
    loc='left',
    va='bottom',
    fontsize=33)
ax_right.set_title(
    label="B. $\\mathbb{E}[s_{X,Y}] = \\frac{|X.A|}{|Q|-|Y.A|}\\times\\sum_{K=1}^{N}s_{K,Y}$",
    loc='left',
    va='bottom',
    fontsize=35)
# ax_right.set_title(
#     label="B. $\\mathbb{E}[s_{X,Y}] = \\frac{|Y|}{\\sum_{K\\neq X}|K|}\\times\\sum_{K=1}^{N}s_{X,K}$",
#     loc='left',
#     va='bottom',
#     fontsize=33)

# Add label to colorbar
fig.text(
    x=0.94, y=0.16,
    s="$\\ln\\left(\\frac{s_{X,Y}+1}{\\mathbb{E}[s_{X,Y}]+1}\\right)$",
    ha='center',
    va='center',
    fontsize=30)

# Seperate by group
y_groups = heatmap_y_axis.to_frame(index=False).groupby(["Class"]).size()
x_groups = heatmap_x_axis.to_frame(index=False).groupby(["Class"]).size()

y_boundaries = np.cumsum(y_groups)
x_boundaries = np.cumsum(x_groups)

for y_loc in y_boundaries[:-1]:
    ax_left.hlines(
        y_loc, 
        xmin=ax_left.get_xlim()[0], 
        xmax=ax_left.get_xlim()[1], 
        colors="white", 
        linewidth=5)
    ax_right.hlines(
        y_loc, 
        xmin=ax_right.get_xlim()[0], 
        xmax=ax_right.get_xlim()[1], 
        colors="white", 
        linewidth=5)
for x_loc in x_boundaries[:-1]:
    ax_left.vlines(
        x_loc, 
        ymin=ax_left.get_ylim()[0], 
        ymax=ax_left.get_ylim()[1], 
        colors="white", 
        linewidth=5)
    ax_right.vlines(
        x_loc, 
        ymin=ax_right.get_ylim()[0], 
        ymax=ax_right.get_ylim()[1], 
        colors="white", 
        linewidth=5)
    
# Add heatmap cell key:
fig.text(
    x=0.94, y=0.96,
    s="Heatmap Cell:",
    fontsize=25,
    ha='center',
    va='bottom')
matplotlib.rc('text', usetex=True)
fig.text(
    x=0.94, y=0.95,
    s=r'''\begin{tabular}{ c | c | c |} & LS & SS \\ \hline 30 & & \\ \hline 50 & & \\ \hline 80 & & \end{tabular}''',
    fontsize=25,
    ha='center',
    va='top')
    
fig.savefig("amr_switch_relative_to_query.png")