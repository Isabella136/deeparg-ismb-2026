from classes.reference import Reference
from classes.query import Query
import pandas as pd
import pickle
import os

CDD_DIR = "../../database/CDD_features_v2/"
REF_LOC = "../../data/database/v2/features.fasta"
CLSTR_LOC = "../../database/database_clustering/db_v2_40.clstr"
DEEPARG_HIT_FILE = "X.mapping.ARG"
ALIGNMENT_FILE = "X.align.daa.tsv"
SAMPLE_ID_FILE = "../real_samples.txt"
PATH_TO_SS = "deeparg_results"
PATH_TO_LS = "spades/deeparg_results"

GRAPH = False

# Domains are sorted by bitscore
def sorted_domains_combo(dom_list: list[str]) -> list[str]:
    if len(dom_list) == 1: return dom_list
    first_dom = dom_list[0]
    recursive_combo_list = sorted_domains_combo(dom_list[1:])
    extended_combo_list = [
        "$".join((first_dom, curr_combo)) for curr_combo in recursive_combo_list]
    extended_combo_list.append(first_dom)
    extended_combo_list.extend(recursive_combo_list)
    return extended_combo_list
    

def main():
    # Read content of real_samples.txt to find biosample ID
    with open(SAMPLE_ID_FILE, "r") as sample_id_buf:
        sample_id_list = sample_id_buf.read().split('\n')
    first = True
    for model in ["LS", "SS"]:    
        # Open metadata to get ARG features and output
        ref_dict: dict[str, Reference] = dict()
        metadata_IO_DL = pickle.load(open("../../data/model/v2/metadata_"+model+".pkl", "rb"))
        for ref in metadata_IO_DL["features"]:
            ref_dict[ref] = Reference(ref)

        # Look through CD features to add domain info
        for part in range(1,26):
            with open(f"{CDD_DIR}Part{part}_hitdata.txt", "r") as cdd_info_buf:
                cdd_info_list = cdd_info_buf.read().split('\n')
            for cdd_info_line in cdd_info_list[1:-1]:
                row = cdd_info_line.split('\t')
                name = row[0].split('>')[-1]
                if name in ref_dict.keys():
                    ref_dict[name].add_domain(row)
            cdd_info_list.clear()

        # Add length information
        with open(REF_LOC, "r") as ref_buf:
            ref_list = ref_buf.read().split('\n')
        for ref_index in range((len(ref_list) - 1) // 2):
            name = ref_list[ref_index*2][1:]
            if name not in ref_dict.keys():
                continue
            ref_dict[name].add_length_info(len(ref_list[ref_index*2+1]))
        ref_list.clear()

        # Add cluster information
        with open(CLSTR_LOC, "r") as clstr_buf:
            clstr_list = clstr_buf.read().split('\n')
        cluster = -1
        for line in clstr_list[:-1]:
            if '>' == line[0]:
                cluster += 1
            else:
                name = line.split('>')[-1].split('...')[0][:-3]
                if name not in ref_dict.keys():
                    continue
                ref_dict[name].define_cluster(cluster)

        # Count the time each label appears in reference
        reference_clstr_amr_count: dict[str, int] = dict()
        reference_arg_amr_count: dict[str, int] = dict()
        reference_dom_amr_count: dict[str, int] = dict()
        reference_super_amr_count: dict[str, int] = dict()
        reference_amr_count: dict[str, int] = dict()
        reference_all_groups_count: dict[str, int] = dict()
        for ref in ref_dict.values():
            (clstr, arg, dom, super, amr) = ref.get_annotations()

            # Count amr, clstr|amr, and arg|amr
            clstr_amr = '|'.join([clstr, amr])
            arg_amr = '|'.join([arg, amr])
            if amr in reference_amr_count.keys():
                reference_amr_count[amr] += 1
            else:
                reference_amr_count[amr] = 1
            if clstr_amr in reference_clstr_amr_count.keys():
                reference_clstr_amr_count[clstr_amr] += 1
            else:
                reference_clstr_amr_count[clstr_amr] = 1
            if arg_amr in reference_arg_amr_count.keys():
                reference_arg_amr_count[arg_amr] += 1
            else:
                reference_arg_amr_count[arg_amr] = 1

            # We may have numerous domains in gene; should count each potential combination as well
            dom_combo_list = sorted_domains_combo(dom.split('$'))
            super_combo_list = sorted_domains_combo(super.split('$'))

            # We also need to count instances of arg as domain
            if dom != arg:
                dom_combo_list.append(arg)
                super_combo_list.append(arg)

            # Count each dom|amr, super|amr, and clstr|arg|dom|super|amr            
            for i in range(len(dom_combo_list)):
                dom_amr = "|".join([dom_combo_list[i], amr])
                super_amr = "|".join([super_combo_list[i], amr])
                all_groups = '|'.join([clstr, arg, dom_combo_list[i], super_combo_list[i], amr])
                
                if dom_amr in reference_dom_amr_count.keys():
                    reference_dom_amr_count[dom_amr] += 1
                else:
                    reference_dom_amr_count[dom_amr] = 1
                if super_amr in reference_super_amr_count.keys():
                    reference_super_amr_count[super_amr] += 1
                else:
                    reference_super_amr_count[super_amr] = 1
                if all_groups not in reference_all_groups_count:
                    reference_all_groups_count[all_groups] = 1
                else:
                    reference_all_groups_count[all_groups] += 1
        
        # Go through each run one at a time
        for sample_id in sample_id_list:
            for identity in [30, 50, 80]:                
                # Get Diamond alignment information
                query_dict: dict[str, Query] = dict()
                with open("/".join((
                        "..", "samples", sample_id, PATH_TO_LS if model=="LS" else PATH_TO_SS,
                        f"arg_alignment_identity_{identity}", ALIGNMENT_FILE)), "r") as alignment_buf:
                    alignment_list = alignment_buf.read().split('\n')
                for alignment in alignment_list[:-1]:
                    row = alignment.split('\t')
                    if (len(query_dict) > 0) and (row[0] in query_dict.keys()):
                        query_dict[row[0]].add_alignment(row, ref_dict, model == "LS")
                    else:
                        query_dict[row[0]] = Query(row, ref_dict, model == "LS", sample_id)
                alignment_list.clear()

                # Get DeepARG hit information
                with open("/".join((
                        "..", "samples", sample_id, PATH_TO_LS if model=="LS" else PATH_TO_SS,
                        f"arg_alignment_identity_{identity}", DEEPARG_HIT_FILE)), "r") as deeparg_buf:
                    deeparg_list = deeparg_buf.read().split('\n')
                for hit in deeparg_list[1:-1]:
                    row = hit.split('\t')
                    
                    # Mark queries that match with read_id (row[3]) as deeparg_hit
                    # Save alignment that matches with best hit (row[5]) in top_deeparg_hit
                    query_dict[row[3]].add_deeparg_hit(row[5])

                    # Count amr switch
                    diamond_amr = query_dict[row[3]].get_top_diamond_classification()
                    deeparg_amr = query_dict[row[3]].get_top_deeparg_classification()
                    if diamond_amr != deeparg_amr:
                        key = "\t".join((diamond_amr, deeparg_amr))
                deeparg_list.clear()

                # Get info on query alignments amr and all groups proportion
                print(f"Sample {sample_id} Alignment identity {identity} Model {model}")
                all_labels_df = pd.DataFrame(
                    columns=[
                        "clstr", "arg", "dom", "super", "amr", "count", 
                        "amr ref count", "clstr|amr ref count",
                        "dom|amr ref count", "super|amr ref count",
                        "Is Diamond Best Hit Label", 
                        "Is Diamond Best-Hit Class",
                        "Is Most Frequent Class",
                        "Is DeepARG Class",
                        "Diamond Class",
                        "Diamond clstr",
                        "Diamond dom",
                        "Diamond super",
                        "Most Frequent Class",
                        "DeepARG Class",
                        "Query"])
                query_count = 0
                for query in query_dict.values():
                    if not query.is_deeparg_hit():
                        continue
                    if not query.has_multiple_classes():
                        continue
                    if query_count % 1000 == 0:
                        print(f"Query # {query_count}")
                    query_count+=1
                    query_vector = query.create_query_vector(
                        reference_clstr_amr_count,
                        reference_dom_amr_count,
                        reference_super_amr_count,
                        reference_amr_count)
                    if not query_vector.has_multiple_possible_classes():
                        continue
                    all_labels_df = pd.concat(
                        [all_labels_df, query_vector.get_label_counts()],
                        ignore_index=True)


                all_labels_df.insert(
                    loc=0, column='Model', value=model)
                all_labels_df.insert(
                    loc=0, column='Alignment Identity', value=identity)
                all_labels_df.insert(
                    loc=0, column='Sample ID', value=sample_id)
                
                if first:
                    all_labels_df.to_csv(
                        path_or_buf=f"label_counts.tsv",
                        sep="\t", float_format='{:.4f}'.format, index=False)
                else:
                    all_labels_df.to_csv(
                        path_or_buf=f"label_counts.tsv",
                        sep="\t", float_format='{:.4f}'.format, mode='a', header=False, index=False)
                    first = False

main()