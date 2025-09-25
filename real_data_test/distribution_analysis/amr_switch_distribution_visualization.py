import matplotlib as plt
import seaborn as sn
import pandas as pd
import numpy as np
from io import StringIO


distribution_file = open("amr_switch_distribution.txt", "r")

# Get class count for DeepARG database, located in first 33 lines of input
class_count = dict()
for line_num in range(33):
    line = distribution_file.readline()
    if line_num == 0 :
        continue;
    amr_class = line.split("\t")[0]
    count = int(line.split("\t")[1])
    class_count.update({amr_class:count})

# Read rest of data into DataFrame
df: pd.DataFrame = pd.read_csv(StringIO(distribution_file.read()), sep='\t', header=0)

# Sort DataFrame 
df.sort_values(
    by=["Model", "Alignment Identity", "Sample", "Diamond AMR", "DeepARG AMR"], 
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
        by=["Diamond AMR", "DeepARG AMR"], key=lambda x: x.str.lower())["Diamond AMR"].values,
     df.get(["Diamond AMR", "DeepARG AMR"]).drop_duplicates().sort_values(
        by=["Diamond AMR", "DeepARG AMR"], key=lambda x: x.str.lower())["DeepARG AMR"].values],
    names=["Diamond AMR", "DeepARG AMR"])

print(headers)

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

    heatmap_df.at[(model, ident, sampl), (diamo, deepa)] = np.log(
        (1 + row[1]["Count"])/(
            1 + (class_count[deepa] / (sum(class_count.values()) - class_count[diamo]) * sum(
                df[df["Model"] == model][df["Alignment Identity"] == ident][df["Sample"] == sampl][df["Diamond AMR"] == diamo]["Count"].values))))
    mask_df.at[(model, ident, sampl), (diamo, deepa)] = False
    max_val = max(max_val, heatmap_df.at[(model, ident, sampl), (diamo, deepa)])
    min_val = min(min_val, heatmap_df.at[(model, ident, sampl), (diamo, deepa)])

from matplotlib.colors import LinearSegmentedColormap

# Create a custom color palette to go with my theme(TM)
custom_colors = [
   '#109368', # first color
   '#d9d9d9', # middle color
   '#DD7719' # last color
]
custom_cmap = LinearSegmentedColormap.from_list("custom_gradient", custom_colors)
custom_cmap.set_bad()

fig = plt.figure.Figure(figsize=[30, 15])
ax = fig.add_axes((0.1, 0.15, 0.8, 0.7))
cbar = fig.add_axes((0.94, 0.2, 0.02, 0.6))

sn.heatmap(
    heatmap_df, mask=mask_df, cmap=custom_cmap, center=0, vmax=max_val, vmin=min_val, 
    square=True, ax=ax, cbar_ax=cbar,
    xticklabels=heatmap_df.columns.get_level_values("DeepARG AMR"),
    yticklabels=heatmap_df.index.get_level_values("Sample"))

ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation_mode="anchor", rotation=45, ha='right')

# Seperate by group
y_groups = heatmap_df.index.to_frame(index=False).groupby(["Model", "Alignment Identity"]).size()
x_groups = heatmap_df.columns.to_frame(index=False).groupby(["Diamond AMR"], sort=False).size()

y_boundaries = np.cumsum(y_groups)
x_boundaries = np.cumsum(x_groups)

for y_loc in y_boundaries[:-1]:
    ax.hlines(y_loc, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], colors="white", linewidth=5)
for x_loc in x_boundaries[:-1]:
    ax.vlines(x_loc, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], colors="white", linewidth=5)

# Add group labels
prev_y_loc = 0
for y_loc, y_group in zip(y_boundaries, y_groups.index):
    ax.text(x_boundaries.iloc[-1] + 0.5, (y_loc+prev_y_loc)/2.0, y_group, va="center")
    prev_y_loc = y_loc

prev_x_loc = 0
for x_loc, x_group in zip(x_boundaries, x_groups.index):
    ax.text((x_loc+prev_x_loc)/2.0, -0.5, x_group, rotation=45)
    prev_x_loc = x_loc

fig.savefig("test.png")