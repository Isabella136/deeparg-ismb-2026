DB40 = "db40.clstr"
DB_V2_40 = "db_v2_40.clstr"

# First parse through DB40, which contains clusters for DeepARG DB v2
clusters_og: dict[int, frozenset[str]] = dict()
with open(DB40) as handle:
    curr_cluster_index = -1
    curr_cluster: set[str] = set()
    for line in handle:
        if line[0] == ">":
            if curr_cluster_index >= 0:
                clusters_og[curr_cluster_index] = frozenset(curr_cluster)
                curr_cluster.clear()
            curr_cluster_index += 1
        else:
            curr_cluster.add(line.split('>')[1].split(' ')[0][:-3] + "|v2")
            
# Then parse through DB_ALL_40, which contains clusters for DeepARG DB v2 but has "v2" concat to seq id
clusters_new: dict[int, frozenset[str]] = dict()
with open(DB_V2_40) as handle:
    curr_cluster_index = -1
    curr_cluster: set[str] = set()
    for line in handle:
        if line[0] == ">":
            if curr_cluster_index >= 0:
                clusters_new[curr_cluster_index] = frozenset(curr_cluster)
                curr_cluster.clear()
            curr_cluster_index += 1
        else:
            curr_cluster.add(line.split('>')[1].split(' ')[0][:-3])

# Now let's check to see if clusters are the same
# Will take O(n^2) time since we don't know if clusters are ordered the same
# But can save time by removing verified clusters from consideration
for curr_cluster_v2 in clusters_og.items():
    verified_cluster_index = None
    dispersed_cluster_indices = list()
    for curr_cluster_all in clusters_new.items():
        if len(curr_cluster_v2[1].symmetric_difference(curr_cluster_all[1])) == 0:
            verified_cluster_index = curr_cluster_all[0]
            break
        elif len(curr_cluster_v2[1].symmetric_difference(curr_cluster_all[1])) != (len(curr_cluster_v2[1]) + len(curr_cluster_all[1])):
            dispersed_cluster_indices.append(f"{curr_cluster_all[0]}")
    if verified_cluster_index != None:
        clusters_new.pop(verified_cluster_index)
    else:
        print(f"OG cluster {curr_cluster_v2[0]} is dispersed in NEW clusters {','.join(dispersed_cluster_indices)}")
        dispersed_cluster_indices.clear()

if len(clusters_new) == 0:
    print("HOORAY!")
                