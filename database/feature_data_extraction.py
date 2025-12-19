from Bio.SeqIO import FastaIO
import pandas as pd
import numpy as np

CDD_V1_DIR = "CDD_features_v1/"
CDD_V2_DIR = "CDD_features_v2/"
DB_V1 = "v1_features.fasta"
DB_V2 = "v2_features.fasta"
CLSTR_V1 = "database_clustering/db_v1_40.clstr"
CLSTR_V2 = "database_clustering/db_v2_40.clstr"
CLSTR_ALL = "database_clustering/db_all_final_40.clstr"

v1_feature_data_cols = pd.Index([
    "id",                                   "original db",                                  "amr class",
    "amr class v1 count",                   "amr class v2 count",                           "v1-only cluster index",
    "v1-only cluster size",                 "v1-only cluster|amr count",                    "v1-v2 cluster index",
    "v1-v2 cluster v1-only size",           "v1-v2 cluster|amr v1-only count",              "corresponding v2 cluster(s) index/indices",
    "corresponding v2 cluster(s) size(s)",  "corresponding v2 cluster(s)|amr counts(s)",    "conserved domain(s) id(s)",
    "conserved domain(s) v1 count(s)",      "conserved domain(s)|amr v1 count(s)",          "conserved domain(s) v2 count(s)",
    "conserved domain(s)|amr v2 count(s)",  "superfamily(ies) id(s)",                       "superfamily(ies) v1 count(s)",
    "superfamily(ies)|amr v1 count(s)",     "superfamily(ies) v2 count(s)",                 "superfamily(ies)|amr v2 count(s)"])

v1_feature_data_series: dict[str, pd.Series] = dict()
v1_amr_class_counts: dict[str, int] = dict()

with open(DB_V1) as handle:
    for record in FastaIO.FastaIterator(handle):
        feature_name = record.id
        original_db = feature_name.split("|")[2]
        amr_class = feature_name.split("|")[3]
        feature_series = pd.Series(
            data=[
                feature_name, original_db, amr_class, 
                -1, -1, -1,
                -1, -1, -1,
                -1, -1, "",
                "", "", "",
                "", "", "",
                "", "", "",
                "", "", "",],
            index=v1_feature_data_cols)
        if amr_class in v1_amr_class_counts:
            v1_amr_class_counts[amr_class] += 1
        else:
            v1_amr_class_counts[amr_class] = 1
        v1_feature_data_series[feature_name[:-3]] = feature_series

v1_cluster_sizes: dict[int, int] = dict()
v1_cluster_amr_counts: dict[str, int] = dict()

with open(CLSTR_V1) as handle:
    curr_clstr_idx = -1
    for line in handle:
        if line[0] == '>':
            if curr_clstr_idx >= 0:
                curr_clstr_idx += 1
            v1_cluster_sizes[curr_clstr_idx] = 0
            continue
        seq_id = line.split('>')[1].split(' ')[0][:-6]
        v1_feature_data_series[seq_id].at["v1-only cluster index"] = curr_clstr_idx
        v1_cluster_sizes[curr_clstr_idx] += 1
        clstr_amr = f"{curr_clstr_idx}|{v1_feature_data_series[seq_id].at["amr class"]}" 
        if clstr_amr in v1_amr_class_counts:
            v1_amr_class_counts[clstr_amr] += 1
        else:
            v1_amr_class_counts[clstr_amr] = 1

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
            dom = row[7]
            dom_amr = f"{dom}|{v1_feature_data_series[seq_id].at["amr class"]}" 
            sup = row[10] if row[1] == "specific" else row[7]
            sup_amr = f"{sup}|{v1_feature_data_series[seq_id].at["amr class"]}" 
            if prev_seq == row[0]:
                v1_feature_data_series[seq_id].at["conserved domain(s) id(s)"] = '$'.join(
                    [v1_feature_data_series[seq_id].at["conserved domain(s) id(s)"], dom])
                v1_feature_data_series[seq_id].at["superfamily(ies) id(s)"] = '$'.join(
                    [v1_feature_data_series[seq_id].at["superfamily(ies) id(s)"], sup])
            else:
                prev_seq = row[0]
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

v2_feature_data_cols = pd.Index([
    "id",                                           "amr class",                                    "amr class v2 count",
    "amr class v1 count",                           "v2-only cluster index",                        "v2-only cluster size",
    "v2-only cluster|amr count",                    "v1-v2 cluster index",                          "v1-v2 cluster v2-only size",
    "v1-v2 cluster|amr v2-only count",              "corresponding v1 cluster(s) index/indices",    "corresponding v1 cluster(s) size(s)",
    "corresponding v1 cluster(s)|amr counts(s)",    "conserved domain(s) id(s)",                    "conserved domain(s) v2 count(s)",
    "conserved domain(s)|amr v2 count(s)",          "conserved domain(s) v1 count(s)",              "conserved domain(s)|amr v1 count(s)",
    "superfamily(ies) id(s)",                       "superfamily(ies) v2 count(s)",                 "superfamily(ies)|amr v2 count(s)",
    "superfamily(ies) v1 count(s)",                 "superfamily(ies)|amr v1 count(s)"])

v2_feature_data_series: dict[str, pd.Series] = dict()
v2_amr_class_counts: dict[str, int] = dict()

with open(DB_V2) as handle:
    for record in FastaIO.FastaIterator(handle):
        feature_name = record.id
        amr_class = feature_name.split("|")[3]
        feature_series = pd.Series(
            data=[
                feature_name, amr_class, -1,
                v1_amr_class_counts.get(amr_class, 0), -1, -1, 
                -1, -1, -1, 
                -1, "", "",
                "", "", "", 
                "", "", "", 
                "", "", "", 
                "", ""],
            index=v2_feature_data_cols)
        if amr_class in v2_amr_class_counts:
            v2_amr_class_counts[amr_class] += 1
        else:
            v2_amr_class_counts[amr_class] = 1
        v2_feature_data_series[feature_name[:-3]] = feature_series

v2_cluster_sizes: dict[int, int] = dict()
v2_cluster_amr_counts: dict[str, int] = dict()

with open(CLSTR_V2) as handle:
    curr_clstr_idx = -1
    for line in handle:
        if line[0] == '>':
            if curr_clstr_idx >= 0:
                curr_clstr_idx += 1
            v2_cluster_sizes[curr_clstr_idx] = 0
            continue
        seq_id = line.split('>')[1].split(' ')[0][:-6]
        v2_feature_data_series[seq_id].at["v2-only cluster index"] = curr_clstr_idx
        v2_cluster_sizes[curr_clstr_idx] += 1
        clstr_amr = f"{curr_clstr_idx}|{v2_feature_data_series[seq_id].at["amr class"]}" 
        if clstr_amr in v2_amr_class_counts:
            v2_amr_class_counts[clstr_amr] += 1
        else:
            v2_amr_class_counts[clstr_amr] = 1

v2_domain_counts: dict[str, int] = dict()
v2_domain_amr_counts: dict[str, int] = dict()
v2_superfamily_counts: dict[str, int] = dict()
v2_superfamily_amr_counts: dict[str, int] = dict()

for part in range(1,27):
    with open(f"{CDD_V2_DIR}Part{part}_hitdata.txt") as handle:
        prev_seq = ""
        for line in handle:
            if '>' not in line:
                continue
            row = line.split('>')[1]
            row_fields = row.split()
            dom = row[7]
            dom_amr = f"{dom}|{v2_feature_data_series[seq_id].at["amr class"]}" 
            sup = row[10] if row[1] == "specific" else row[7]
            sup_amr = f"{sup}|{v2_feature_data_series[seq_id].at["amr class"]}" 
            if prev_seq == row[0]:
                v2_feature_data_series[seq_id].at["conserved domain(s) id(s)"] = '$'.join(
                    [v2_feature_data_series[seq_id].at["conserved domain(s) id(s)"], dom])
                v2_feature_data_series[seq_id].at["conserved domain(s) v1 count(s)"] = '$'.join([
                    v2_feature_data_series[seq_id].at["conserved domain(s) v1 count(s)"],
                    str(v1_domain_counts[dom])])
                v2_feature_data_series[seq_id].at["conserved domain(s)|amr v1 count(s)"] = '$'.join([
                    v2_feature_data_series[seq_id].at["conserved domain(s)|amr v1 count(s)"],
                    str(v1_domain_amr_counts[dom_amr])])
                v2_feature_data_series[seq_id].at["superfamily(ies) id(s)"] = '$'.join(
                    [v2_feature_data_series[seq_id].at["superfamily(ies) id(s)"], sup])
                v2_feature_data_series[seq_id].at["superfamily(ies) v1 count(s)"] = '$'.join([
                    v2_feature_data_series[seq_id].at["superfamily(ies) v1 count(s)"],
                    str(v1_superfamily_counts[sup])])
                v2_feature_data_series[seq_id].at["superfamily(ies)|amr v1 count(s)"] = '$'.join([
                    v2_feature_data_series[seq_id].at["superfamily(ies)|amr v1 count(s)"],
                    str(v1_superfamily_amr_counts[sup_amr])])
            else:
                prev_seq = row[0]
                v2_feature_data_series[seq_id].at["conserved domain(s) id(s)"] = dom
                v2_feature_data_series[seq_id].at["conserved domain(s) v1 count(s)"] = str(v1_domain_counts[dom])
                v2_feature_data_series[seq_id].at["conserved domain(s)|amr v1 count(s)"] = str(v1_domain_amr_counts[dom_amr])
                v2_feature_data_series[seq_id].at["superfamily(ies) id(s)"] = sup
                v2_feature_data_series[seq_id].at["superfamily(ies) v1 count(s)"] = str(v1_superfamily_counts[sup])
                v2_feature_data_series[seq_id].at["superfamily(ies)|amr v1 count(s)"] = str(v1_superfamily_amr_counts[sup_amr])
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

v1_v2_cluster_v1_sizes: dict[int, int] = dict()
v1_v2_cluster_amr_v1_counts: dict[str, int] = dict()
v1_v2_cluster_v2_sizes: dict[int, int] = dict()
v1_v2_cluster_amr_v2_counts: dict[str, int] = dict()

with open(CLSTR_ALL) as handle:
    curr_clstr_idx = -1
    for line in handle:
        if line[0] == '>':
            if curr_clstr_idx >= 0:
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