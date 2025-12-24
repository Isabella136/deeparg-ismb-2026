from Bio.SeqIO import FastaIO
import pandas as pd

CDD_V1_DIR = "CDD_features_v1/"
CDD_V2_DIR = "CDD_features_v2/"
DB_V1 = "v1_features.fasta"
DB_V2 = "v2_features.fasta"
CLSTR_V1 = "database_clustering/db_v1_40.clstr"
CLSTR_V2 = "database_clustering/db_v2_40.clstr"
CLSTR_ALL = "database_clustering/db_all_final_40.clstr"

# Define column names for v1 table
v1_feature_data_cols = pd.Index([
    "original db",
    "amr class",
    "amr class v1 count",
    "amr class v2 count",
    "v1-only cluster index",
    "v1-only cluster size",
    "v1-only cluster|amr count",
    "v1-v2 cluster index",
    "v1-v2 cluster v1 size",
    "v1-v2 cluster|amr v1 count",
    "corresponding v2-only cluster(s) index/indices",
    "corresponding v2-only cluster(s) size(s)",
    "corresponding v2-only cluster(s)|amr counts(s)",
    "conserved domain(s) id(s)",
    "conserved domain(s) v1 count(s)",
    "conserved domain(s)|amr v1 count(s)",
    "conserved domain(s) v2 count(s)",
    "conserved domain(s)|amr v2 count(s)",
    "superfamily(ies) id(s)",
    "superfamily(ies) v1 count(s)",
    "superfamily(ies)|amr v1 count(s)",
    "superfamily(ies) v2 count(s)",
    "superfamily(ies)|amr v2 count(s)"])

# Create dictionary of v1 features
# Save amr class v1 distribution information
v1_feature_data_series: dict[str, pd.Series] = dict()
v1_amr_class_counts: dict[str, int] = dict()
with open(DB_V1) as handle:
    for record in FastaIO.FastaIterator(handle):
        feature_name = record.id[:-3]
        original_db = feature_name.split("|")[2]
        amr_class = feature_name.split("|")[3]
        match amr_class:
            case "beta_lactam": amr_class = "beta-lactam"
            case "macrolide-lincosamide-streptogramin": amr_class = "MLS"
            case "chloramphenicol": amr_class = "chloramphenicol-phenicol"
            case "quinolone": amr_class = "quinolone-fluoroquinolone"
            case "rifampin": amr_class = "rifampin-rifamycin"
            case "fusidic_acid": amr_class = "fusidic-acid"
        feature_series = pd.Series(
            data={
                "original db" : original_db, 
                "amr class" : amr_class, 
                "amr class v1 count" : 0, 
                "amr class v2 count" : 0, 
                "v1-only cluster index" : -1,
                "v1-only cluster size" : 0, 
                "v1-only cluster|amr count" : 0, 
                "v1-v2 cluster index" : -1,
                "v1-v2 cluster v1 size" : 0, 
                "v1-v2 cluster|amr v1 count" : 0, 
                "corresponding v2-only cluster(s) index/indices" : "",
                "corresponding v2-only cluster(s) size(s)" : "", 
                "corresponding v2-only cluster(s)|amr counts(s)" : "", 
                "conserved domain(s) id(s)" : "",
                "conserved domain(s) v1 count(s)" : "", 
                "conserved domain(s)|amr v1 count(s)" : "",
                "conserved domain(s) v2 count(s)" : "",
                "conserved domain(s)|amr v2 count(s)" : "", 
                "superfamily(ies) id(s)" : "",
                "superfamily(ies) v1 count(s)" : "",
                "superfamily(ies)|amr v1 count(s)" : "", 
                "superfamily(ies) v2 count(s)" : "",
                "superfamily(ies)|amr v2 count(s)" : ""},
            index=v1_feature_data_cols)
        if amr_class in v1_amr_class_counts:
            v1_amr_class_counts[amr_class] += 1
        else:
            v1_amr_class_counts[amr_class] = 1
        v1_feature_data_series[feature_name] = feature_series

# Extract v1-only cluster index for each v1 feature 
# Save v1-only cluster and v1-only cluster|amr distr. info
v1_cluster_sizes: dict[int, int] = dict()
v1_cluster_amr_counts: dict[str, int] = dict()
with open(CLSTR_V1) as handle:
    curr_clstr_idx = -1
    for line in handle:
        if line[0] == '>':
            curr_clstr_idx += 1
            v1_cluster_sizes[curr_clstr_idx] = 0
            continue
        seq_id = line.split('>')[1].split(' ')[0][:-6]
        v1_feature_data_series[seq_id].at["v1-only cluster index"] = curr_clstr_idx
        v1_cluster_sizes[curr_clstr_idx] += 1
        clstr_amr = f"{curr_clstr_idx}|{v1_feature_data_series[seq_id].at["amr class"]}" 
        if clstr_amr in v1_cluster_amr_counts:
            v1_cluster_amr_counts[clstr_amr] += 1
        else:
            v1_cluster_amr_counts[clstr_amr] = 1

# Extract domain(s) and superfamily(ies) for each v1 feature 
# Save domain, domain|amr, superfamily, and superfamily|amr v1 distr. info
v1_domain_counts: dict[str, int] = dict()
v1_domain_amr_counts: dict[str, int] = dict()
v1_superfamily_counts: dict[str, int] = dict()
v1_superfamily_amr_counts: dict[str, int] = dict()
for part in range(1,31):
    with open(f"{CDD_V1_DIR}Part{part}_hitdata.txt") as handle:
        prev_seq = ""
        for line in handle:
            if '>' not in line:
                continue
            row = line.split('>')[1]
            row_fields = row.split()
            seq_id = row_fields[0]
            dom = row_fields[7]
            dom_amr = f"{dom}|{v1_feature_data_series[seq_id].at["amr class"]}" 
            sup = row_fields[10] if row_fields[1] == "specific" else row_fields[7]
            sup_amr = f"{sup}|{v1_feature_data_series[seq_id].at["amr class"]}" 
            if prev_seq == seq_id:
                v1_feature_data_series[seq_id].at["conserved domain(s) id(s)"] = '$'.join(
                    [v1_feature_data_series[seq_id].at["conserved domain(s) id(s)"], dom])
                v1_feature_data_series[seq_id].at["superfamily(ies) id(s)"] = '$'.join(
                    [v1_feature_data_series[seq_id].at["superfamily(ies) id(s)"], sup])
            else:
                prev_seq = seq_id
                v1_feature_data_series[seq_id].at["conserved domain(s) id(s)"] = dom
                v1_feature_data_series[seq_id].at["superfamily(ies) id(s)"] = sup
            if dom in v1_domain_counts:
                v1_domain_counts[dom] += 1
            else:
                v1_domain_counts[dom] = 1
            if dom_amr in v1_domain_amr_counts:
                v1_domain_amr_counts[dom_amr] += 1
            else:
                v1_domain_amr_counts[dom_amr] = 1
            if sup in v1_superfamily_counts:
                v1_superfamily_counts[sup] += 1
            else:
                v1_superfamily_counts[sup] = 1
            if sup_amr in v1_superfamily_amr_counts:
                v1_superfamily_amr_counts[sup_amr] += 1
            else:
                v1_superfamily_amr_counts[sup_amr] = 1

# Define column names for v2 table
v2_feature_data_cols = pd.Index([
    "amr class",
    "amr class v2 count",
    "amr class v1 count",
    "v2-only cluster index",
    "v2-only cluster size",
    "v2-only cluster|amr count",
    "v1-v2 cluster index",
    "v1-v2 cluster v2 size",
    "v1-v2 cluster|amr v2 count",
    "corresponding v1-only cluster(s) index/indices",
    "corresponding v1-only cluster(s) size(s)",
    "corresponding v1-only cluster(s)|amr counts(s)",
    "conserved domain(s) id(s)",
    "conserved domain(s) v2 count(s)",
    "conserved domain(s)|amr v2 count(s)",
    "conserved domain(s) v1 count(s)",
    "conserved domain(s)|amr v1 count(s)",
    "superfamily(ies) id(s)",
    "superfamily(ies) v2 count(s)",
    "superfamily(ies)|amr v2 count(s)",
    "superfamily(ies) v1 count(s)",
    "superfamily(ies)|amr v1 count(s)"])

# Create dictionary of v2 features and extract amr class v1 distr info
# Save amr class v2 distribution information
v2_feature_data_series: dict[str, pd.Series] = dict()
v2_amr_class_counts: dict[str, int] = dict()
with open(DB_V2) as handle:
    for record in FastaIO.FastaIterator(handle):
        feature_name = record.id[:-3]
        amr_class = feature_name.split("|")[3]
        match amr_class:
            case "phenicol": amr_class = "chloramphenicol-phenicol"
            case "fluoroquinolone": amr_class = "quinolone-fluoroquinolone"
            case "rifamycin": amr_class = "rifampin-rifamycin"
            case "tetracenomycin_C": amr_class = "tetracenomycin"
        feature_series = pd.Series(
            data={
                "amr class" : amr_class, 
                "amr class v2 count" : 0, 
                "amr class v1 count" : v1_amr_class_counts.get(amr_class, 0), 
                "v2-only cluster index" : -1,
                "v2-only cluster size" : 0,
                "v2-only cluster|amr count" : 0,
                "v1-v2 cluster index" : -1,
                "v1-v2 cluster v2 size" : 0,
                "v1-v2 cluster|amr v2 count" : 0,
                "corresponding v1-only cluster(s) index/indices" : "",
                "corresponding v1-only cluster(s) size(s)" : "",
                "corresponding v1-only cluster(s)|amr counts(s)" : "",
                "conserved domain(s) id(s)" : "",
                "conserved domain(s) v2 count(s)" : "",
                "conserved domain(s)|amr v2 count(s)" : "",
                "conserved domain(s) v1 count(s)" : "",
                "conserved domain(s)|amr v1 count(s)" : "",
                "superfamily(ies) id(s)" : "",
                "superfamily(ies) v2 count(s)" : "",
                "superfamily(ies)|amr v2 count(s)" : "",
                "superfamily(ies) v1 count(s)" : "",
                "superfamily(ies)|amr v1 count(s)" : ""},
            index=v2_feature_data_cols)
        if amr_class in v2_amr_class_counts:
            v2_amr_class_counts[amr_class] += 1
        else:
            v2_amr_class_counts[amr_class] = 1
        v2_feature_data_series[feature_name] = feature_series

# Extract v2-only cluster index for each v2 feature 
# Save v2-only cluster and v2-only cluster|amr distr. info
v2_cluster_sizes: dict[int, int] = dict()
v2_cluster_amr_counts: dict[str, int] = dict()
with open(CLSTR_V2) as handle:
    curr_clstr_idx = -1
    for line in handle:
        if line[0] == '>':
            curr_clstr_idx += 1
            v2_cluster_sizes[curr_clstr_idx] = 0
            continue
        seq_id = line.split('>')[1].split(' ')[0][:-6]
        v2_feature_data_series[seq_id].at["v2-only cluster index"] = curr_clstr_idx
        v2_cluster_sizes[curr_clstr_idx] += 1
        clstr_amr = f"{curr_clstr_idx}|{v2_feature_data_series[seq_id].at["amr class"]}" 
        if clstr_amr in v2_cluster_amr_counts:
            v2_cluster_amr_counts[clstr_amr] += 1
        else:
            v2_cluster_amr_counts[clstr_amr] = 1

# Extract domain(s) and superfamily(ies) for each v2 feature 
# Extract domain, domain|amr, superfamily, and superfamily|amr v1 distr. info
# Save domain, domain|amr, superfamily, and superfamily|amr v2 distr. info
v2_domain_counts: dict[str, int] = dict()
v2_domain_amr_counts: dict[str, int] = dict()
v2_superfamily_counts: dict[str, int] = dict()
v2_superfamily_amr_counts: dict[str, int] = dict()
for part in range(1,26):
    with open(f"{CDD_V2_DIR}Part{part}_hitdata.txt") as handle:
        prev_seq = ""
        for line in handle:
            if '>' not in line:
                continue
            row = line.split('>')[1]
            row_fields = row.split()
            seq_id = row_fields[0]
            dom = row_fields[7]
            dom_amr = f"{dom}|{v2_feature_data_series[seq_id].at["amr class"]}" 
            sup = row_fields[10] if row_fields[1] == "specific" else row_fields[7]
            sup_amr = f"{sup}|{v2_feature_data_series[seq_id].at["amr class"]}" 
            if prev_seq == seq_id:
                v2_feature_data_series[seq_id].at["conserved domain(s) id(s)"] = '$'.join(
                    [v2_feature_data_series[seq_id].at["conserved domain(s) id(s)"], dom])
                v2_feature_data_series[seq_id].at["conserved domain(s) v1 count(s)"] = '$'.join([
                    v2_feature_data_series[seq_id].at["conserved domain(s) v1 count(s)"],
                    str(v1_domain_counts.get(dom, 0))])
                v2_feature_data_series[seq_id].at["conserved domain(s)|amr v1 count(s)"] = '$'.join([
                    v2_feature_data_series[seq_id].at["conserved domain(s)|amr v1 count(s)"],
                    str(v1_domain_amr_counts.get(dom_amr, 0))])
                v2_feature_data_series[seq_id].at["superfamily(ies) id(s)"] = '$'.join(
                    [v2_feature_data_series[seq_id].at["superfamily(ies) id(s)"], sup])
                v2_feature_data_series[seq_id].at["superfamily(ies) v1 count(s)"] = '$'.join([
                    v2_feature_data_series[seq_id].at["superfamily(ies) v1 count(s)"],
                    str(v1_superfamily_counts.get(sup, 0))])
                v2_feature_data_series[seq_id].at["superfamily(ies)|amr v1 count(s)"] = '$'.join([
                    v2_feature_data_series[seq_id].at["superfamily(ies)|amr v1 count(s)"],
                    str(v1_superfamily_amr_counts.get(sup_amr, 0))])
            else:
                prev_seq = seq_id
                v2_feature_data_series[seq_id].at[
                    "conserved domain(s) id(s)"] = dom
                v2_feature_data_series[seq_id].at[
                    "conserved domain(s) v1 count(s)"] = str(v1_domain_counts.get(dom, 0))
                v2_feature_data_series[seq_id].at[
                    "conserved domain(s)|amr v1 count(s)"] = str(v1_domain_amr_counts.get(dom_amr, 0))
                v2_feature_data_series[seq_id].at[
                    "superfamily(ies) id(s)"] = sup
                v2_feature_data_series[seq_id].at[
                    "superfamily(ies) v1 count(s)"] = str(v1_superfamily_counts.get(sup, 0))
                v2_feature_data_series[seq_id].at[
                    "superfamily(ies)|amr v1 count(s)"] = str(v1_superfamily_amr_counts.get(sup_amr, 0))
            if dom in v2_domain_counts:
                v2_domain_counts[dom] += 1
            else:
                v2_domain_counts[dom] = 1
            if dom_amr in v2_domain_amr_counts:
                v2_domain_amr_counts[dom_amr] += 1
            else:
                v2_domain_amr_counts[dom_amr] = 1
            if sup in v2_superfamily_counts:
                v2_superfamily_counts[sup] += 1
            else:
                v2_superfamily_counts[sup] = 1
            if sup_amr in v2_superfamily_amr_counts:
                v2_superfamily_amr_counts[sup_amr] += 1
            else:
                v2_superfamily_amr_counts[sup_amr] = 1

# Extract v1-v2 cluster index for each feature 
# Save v1-v2 cluster and v1-v2 cluster|amr v1 distr. info
# Save v1-v2 cluster and v1-v2 cluster|amr v2 distr. info
v1_v2_cluster_v1_sizes: dict[int, int] = dict()
v1_v2_cluster_amr_v1_counts: dict[str, int] = dict()
v1_v2_cluster_v2_sizes: dict[int, int] = dict()
v1_v2_cluster_amr_v2_counts: dict[str, int] = dict()
with open(CLSTR_ALL) as handle:
    curr_clstr_idx = -1
    for line in handle:
        if line[0] == '>':
            curr_clstr_idx += 1
            v1_v2_cluster_v1_sizes[curr_clstr_idx] = 0
            v1_v2_cluster_v2_sizes[curr_clstr_idx] = 0
            continue
        
        seq_id_ver = line.split('>')[1].split(' ')[0][:-3]
        ver = seq_id_ver.split('|')[-1]
        seq_id = seq_id_ver[:-3]

        if ver == "v1":
            amr = v1_feature_data_series[seq_id].at['amr class']
            v1_feature_data_series[seq_id].at["v1-v2 cluster index"] = curr_clstr_idx
            v1_v2_cluster_v1_sizes[curr_clstr_idx] += 1
            clstr_amr = f"{curr_clstr_idx}|{amr}" 
            if clstr_amr in v1_v2_cluster_amr_v1_counts:
                v1_v2_cluster_amr_v1_counts[clstr_amr] += 1
            else:
                v1_v2_cluster_amr_v1_counts[clstr_amr] = 1
        else:
            amr = v2_feature_data_series[seq_id].at['amr class']
            v2_feature_data_series[seq_id].at["v1-v2 cluster index"] = curr_clstr_idx
            v1_v2_cluster_v2_sizes[curr_clstr_idx] += 1
            clstr_amr = f"{curr_clstr_idx}|{amr}" 
            if clstr_amr in v1_v2_cluster_amr_v2_counts:
                v1_v2_cluster_amr_v2_counts[clstr_amr] += 1
            else:
                v1_v2_cluster_amr_v2_counts[clstr_amr] = 1

# Make v1 and v2 feature dataframes, with seq_id being the index:
v1_feature_df = pd.DataFrame(
    data=list(v1_feature_data_series.values()),
    index=list(v1_feature_data_series.keys()),
    columns=v1_feature_data_cols)

v2_feature_df = pd.DataFrame(
    data=list(v2_feature_data_series.values()),
    index=list(v2_feature_data_series.keys()),
    columns=v2_feature_data_cols)

# As of now, "amr class v1 count" is undefined in the v1 feature dataframe
# And "amr class v2 count" is undefined in both the v1 and v2 feature dataframes
# Can retrieve amr class counts from v1_amr_class_counts and v2_amr_class_counts
print("Get amr class counts")
v1_feature_df["amr class v1 count"] = v1_feature_df.apply(
    lambda x: v1_amr_class_counts.get(x["amr class"], 0), axis=1)
v1_feature_df["amr class v2 count"] = v1_feature_df.apply(
    lambda x: v2_amr_class_counts.get(x["amr class"], 0), axis=1)
v2_feature_df["amr class v2 count"] = v2_feature_df.apply(
    lambda x: v2_amr_class_counts.get(x["amr class"], 0), axis=1)

# Now, "v1-only cluster size" is undefined in the v1 feature dataframe
# And "v2-only cluster size" is undefined in the v2 feature dataframe
# Can retrieve size from v1_cluster_sizes and v2_cluster_sizes
print("Get ver-only cluster size")
v1_feature_df["v1-only cluster size"] = v1_feature_df.apply(
    lambda x: v1_cluster_sizes[x["v1-only cluster index"]], axis=1)
v2_feature_df["v2-only cluster size"] = v2_feature_df.apply(
    lambda x: v2_cluster_sizes[x["v2-only cluster index"]], axis=1)

# Now, "v1-only cluster|amr count" is undefined in the v1 feature dataframe
# And "v2-only cluster|amr count" is undefined in the v2 feature dataframe
# Can retrieve count from v1_cluster_amr_counts and v2_cluster_amr_counts
print("Get ver-only cluster|amr counts")
v1_feature_df["v1-only cluster|amr count"] = v1_feature_df.apply(
    lambda x: v1_cluster_amr_counts['|'.join(
        (str(x["v1-only cluster index"]), x["amr class"]))], axis=1)
v2_feature_df["v2-only cluster|amr count"] = v2_feature_df.apply(
    lambda x: v2_cluster_amr_counts['|'.join(
        (str(x["v2-only cluster index"]), x["amr class"]))], axis=1)

# Now, "v1-v2 cluster v1 size" is undefined in the v1 feature dataframe
# And "v1-v2 cluster v2 size" is undefined in the v2 feature dataframe
# Can retrieve size from v1_v2_cluster_v1_sizes and v1_v2_cluster_v2_sizes
print("Get all-ver cluster size")
v1_feature_df["v1-v2 cluster v1 size"] = v1_feature_df.apply(
    lambda x: v1_v2_cluster_v1_sizes[x["v1-v2 cluster index"]], axis=1)
v2_feature_df["v1-v2 cluster v2 size"] = v2_feature_df.apply(
    lambda x: v1_v2_cluster_v2_sizes[x["v1-v2 cluster index"]], axis=1)

# Now, "v1-v2 cluster|amr v1 count" is undefined in the v1 feature dataframe
# And "v1-v2 cluster|amr v2 count" is undefined in the v2 feature dataframe
# Can retrieve count from v1_v2_cluster_amr_v1_counts and v1_v2_cluster_amr_v2_counts
print("Get all-ver cluster|amr counts")
v1_feature_df["v1-v2 cluster|amr v1 count"] = v1_feature_df.apply(
    lambda x: v1_v2_cluster_amr_v1_counts['|'.join(
        (str(x["v1-v2 cluster index"]), x["amr class"]))], axis=1)
v2_feature_df["v1-v2 cluster|amr v2 count"] = v2_feature_df.apply(
    lambda x: v1_v2_cluster_amr_v2_counts['|'.join(
        (str(x["v1-v2 cluster index"]), x["amr class"]))], axis=1)

# Now, "corresponding v2-only cluster(s) index/indices" is undefined in the v1 feature dataframe
# And "corresponding v1-only cluster(s) index/indices" is undefined in the v2 feature dataframe
# Can retrieve index from other feature dataframe
print("Get corresponding cluster")
v1_feature_df["corresponding v2-only cluster(s) index/indices"] = v1_feature_df.apply(
    lambda x: '$'.join(str(idx) for idx in v2_feature_df.loc[
        v2_feature_df["v1-v2 cluster index"] == x["v1-v2 cluster index"]
    ]["v2-only cluster index"].drop_duplicates().to_list()), axis=1)
v2_feature_df["corresponding v1-only cluster(s) index/indices"] = v2_feature_df.apply(
    lambda x: '$'.join(str(idx) for idx in v1_feature_df.loc[
        v1_feature_df["v1-v2 cluster index"] == x["v1-v2 cluster index"]
    ]["v1-only cluster index"].drop_duplicates().to_list()), axis=1)

# Now, "corresponding v2-only cluster(s) size(s)" is undefined in the v1 feature dataframe
# And "corresponding v1-only cluster(s) size(s)" is undefined in the v2 feature dataframe
# Can retrieve size from v1_cluster_sizes and v2_cluster_sizes
print("Get corresponding cluster size")
v1_feature_df["corresponding v2-only cluster(s) size(s)"] = v1_feature_df.apply(
    lambda x: '$'.join(
        str(v2_cluster_sizes.get(int("-1" if idx == "" else idx), 0)) for idx 
        in x["corresponding v2-only cluster(s) index/indices"].split('$')), axis=1)
v2_feature_df["corresponding v1-only cluster(s) size(s)"] = v2_feature_df.apply(
    lambda x: '$'.join(
        str(v1_cluster_sizes.get(int("-1" if idx == "" else idx), 0)) for idx 
        in x["corresponding v1-only cluster(s) index/indices"].split('$')), axis=1)

# Now, "corresponding v2-only cluster(s)|amr counts(s)" is undefined in the v1 feature dataframe
# And "corresponding v1-only cluster(s)|amr counts(s)" is undefined in the v2 feature dataframe
# Can retrieve count from v1_cluster_amr_counts and v2_cluster_amr_counts
print("Get corresponding cluster|amr count")
v1_feature_df["corresponding v2-only cluster(s)|amr counts(s)"] = v1_feature_df.apply(
    lambda x: '$'.join(
        str(v2_cluster_amr_counts.get('|'.join((idx, x["amr class"])), 0)) for idx 
        in x["corresponding v2-only cluster(s) index/indices"].split('$')), axis=1)
v2_feature_df["corresponding v1-only cluster(s)|amr counts(s)"] = v2_feature_df.apply(
    lambda x: '$'.join(
        str(v1_cluster_amr_counts.get('|'.join((idx, x["amr class"])), 0)) for idx 
        in x["corresponding v1-only cluster(s) index/indices"].split('$')), axis=1)

# Now, "conserved domain(s) v1 count(s)" is undefined in the v1 feature dataframe
# And "conserved domain(s) v2 count(s)" is undefined in both feature dataframes
# Can retrieve count from v1_domain_counts and v2_domain_counts
print("Get domain count")
v1_feature_df["conserved domain(s) v1 count(s)"] = v1_feature_df.apply(
    lambda x: '$'.join(str(v1_domain_counts.get(id, 0)) 
        for id in x["conserved domain(s) id(s)"].split('$')), axis=1)
v1_feature_df["conserved domain(s) v2 count(s)"] = v1_feature_df.apply(
    lambda x: '$'.join(str(v2_domain_counts.get(id, 0)) 
        for id in x["conserved domain(s) id(s)"].split('$')), axis=1)
v2_feature_df["conserved domain(s) v2 count(s)"] = v2_feature_df.apply(
    lambda x: '$'.join(str(v2_domain_counts.get(id, 0)) 
        for id in x["conserved domain(s) id(s)"].split('$')), axis=1)

# Now, "conserved domain(s)|amr v1 count(s)" is undefined in the v1 feature dataframe
# And "conserved domain(s)|amr v2 count(s)" is undefined in both feature dataframes
# Can retrieve count from v1_domain_amr_counts and v2_domain_amr_counts
print("Get domain|amr count")
v1_feature_df["conserved domain(s)|amr v1 count(s)"] = v1_feature_df.apply(
    lambda x: '$'.join(str(v1_domain_amr_counts.get('|'.join((id, x["amr class"])), 0))
        for id in x["conserved domain(s) id(s)"].split('$')), axis=1)
v1_feature_df["conserved domain(s)|amr v2 count(s)"] = v1_feature_df.apply(
    lambda x: '$'.join(str(v2_domain_amr_counts.get('|'.join((id, x["amr class"])), 0)) 
        for id in x["conserved domain(s) id(s)"].split('$')), axis=1)
v2_feature_df["conserved domain(s)|amr v2 count(s)"] = v2_feature_df.apply(
    lambda x: '$'.join(str(v2_domain_amr_counts.get('|'.join((id, x["amr class"])), 0)) 
        for id in x["conserved domain(s) id(s)"].split('$')), axis=1)

# Now, "superfamily(ies) v1 count(s)" is undefined in the v1 feature dataframe
# And "superfamily(ies) v2 count(s)" is undefined in both feature dataframes
# Can retrieve count from v1_superfamily_counts and v2_superfamily_counts
print("Get superfamily count")
v1_feature_df["superfamily(ies) v1 count(s)"] = v1_feature_df.apply(
    lambda x: '$'.join(str(v1_superfamily_counts.get(id, 0)) 
        for id in x["superfamily(ies) id(s)"].split('$')), axis=1)
v1_feature_df["superfamily(ies) v2 count(s)"] = v1_feature_df.apply(
    lambda x: '$'.join(str(v2_superfamily_counts.get(id, 0)) 
        for id in x["superfamily(ies) id(s)"].split('$')), axis=1)
v2_feature_df["superfamily(ies) v2 count(s)"] = v2_feature_df.apply(
    lambda x: '$'.join(str(v2_superfamily_counts.get(id, 0)) 
        for id in x["superfamily(ies) id(s)"].split('$')), axis=1)

# Now, "superfamily(ies)|amr v1 count(s)" is undefined in the v1 feature dataframe
# And "superfamily(ies)|amr v2 count(s)" is undefined in both feature dataframes
# Can retrieve count from v1_superfamily_amr_counts and v2_superfamily_amr_counts
print("Get superfamily|amr count")
v1_feature_df["superfamily(ies)|amr v1 count(s)"] = v1_feature_df.apply(
    lambda x: '$'.join(str(v1_superfamily_amr_counts.get('|'.join((id, x["amr class"])), 0))
        for id in x["superfamily(ies) id(s)"].split('$')), axis=1)
v1_feature_df["superfamily(ies)|amr v2 count(s)"] = v1_feature_df.apply(
    lambda x: '$'.join(str(v2_superfamily_amr_counts.get('|'.join((id, x["amr class"])), 0)) 
        for id in x["superfamily(ies) id(s)"].split('$')), axis=1)
v2_feature_df["superfamily(ies)|amr v2 count(s)"] = v2_feature_df.apply(
    lambda x: '$'.join(str(v2_superfamily_amr_counts.get('|'.join((id, x["amr class"])), 0)) 
        for id in x["superfamily(ies) id(s)"].split('$')), axis=1)

# Everything should now be written to dataframes, so we are good to output
v1_feature_df.to_csv("v1_feature_data.csv")
v2_feature_df.to_csv("v2_feature_data.csv")