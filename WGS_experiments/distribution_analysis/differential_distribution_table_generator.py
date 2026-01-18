from differential_distribution_classes.graph_items import ClstrVertex, DomainVertex, SuperVertex, AmrVertex, ArgVertex, MultiVertex
from differential_distribution_classes.reference import Reference
from differential_distribution_classes.graph import Graph
from differential_distribution_classes.query import Query
import pandas as pd
import numpy as np
import pickle
import sys
import os

CDD_DIR = "../../database/CDD_features_v2/"
REF_LOC = "../../data/database/v2/features.fasta"
CLSTR_LOC = "../../database/database_clustering/db_v2_40.clstr"
DEEPARG_HIT_FILE = "X.mapping.ARG"
ALIGNMENT_FILE = "X.align.daa.tsv"
SAMPLE_ID_FILE = "../real_samples.txt"
PATH_TO_SS = "deeparg_results"
PATH_TO_LS = "spades/deeparg_results"

AMR_SWITCH_DISTRIBUTION = "amr_switch_distribution.txt"
OUTPUT_LOC = "differential_distribution_output"
CLUSTER_OUTPUT = "cluster.tsv"
ARG_OUTPUT = "arg.tsv"
AMR_OUTPUT = "amr.tsv"
DOMAIN_OUTPUT = "domain.tsv"
SUPER_OUTPUT = "super.tsv"

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
            (clstr, arg, dom, super, amr) = ref.get_groupings()

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

        # Check to see if output directory exist, or make one
        if not os.path.exists(OUTPUT_LOC):
            os.mkdir(OUTPUT_LOC)
        
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
                        path_or_buf=f"{OUTPUT_LOC}/label_counts.tsv",
                        sep="\t", float_format='{:.4f}'.format, index=False)
                else:
                    all_labels_df.to_csv(
                        path_or_buf=f"{OUTPUT_LOC}/label_counts.tsv",
                        sep="\t", float_format='{:.4f}'.format, mode='a', header=False, index=False)
                if not GRAPH:
                    first = False
                
                # Create vertices and edges to represent labels and pairs, respectively
                if GRAPH:
                    amr_vertices: dict[str, AmrVertex] = dict()
                    clstr_vertices: dict[str, ClstrVertex] = dict()
                    arg_vertices: dict[str, ArgVertex] = dict()
                    domain_vertices: dict[str, DomainVertex] = dict()
                    super_vertices: dict[str, SuperVertex] = dict()
                    all_groups_vertices: dict[str, MultiVertex] = dict()
                    for query in query_dict.values():
                        if query.is_deeparg_hit():
                            (clstr_diamond, arg_diamond, dom_diamond, super_diamond, amr_diamond) = (
                                query.get_top_diamond_alignment().get_groupings())
                            is_alignment_domain_less = len(query.get_top_diamond_alignment().get_domains()) == 0
                            
                            all_groups_diamond = '|'.join(query.get_top_diamond_alignment().get_groupings())
                            clstr_amr_diamond = '|'.join([clstr_diamond, amr_diamond])
                            arg_amr_diamond = '|'.join([arg_diamond, amr_diamond])
                            dom_amr_diamond = '|'.join([dom_diamond, amr_diamond])
                            super_amr_diamond = '|'.join([super_diamond, amr_diamond])

                            if all_groups_diamond not in all_groups_vertices.keys():
                                all_groups_vertices[all_groups_diamond] = MultiVertex(
                                    all_groups_diamond, 
                                    reference_all_groups_count[all_groups_diamond], 
                                    is_alignment_domain_less)
                            
                            if amr_diamond not in amr_vertices.keys():
                                amr_vertices[amr_diamond] = AmrVertex(
                                    amr_diamond, reference_amr_count[amr_diamond])
                            
                            if clstr_amr_diamond not in clstr_vertices.keys():
                                clstr_vertices[clstr_amr_diamond] = ClstrVertex(
                                    clstr_amr_diamond, reference_clstr_amr_count[clstr_amr_diamond])
                            
                            if arg_amr_diamond not in arg_vertices.keys():
                                arg_vertices[arg_amr_diamond] = ArgVertex(
                                    arg_amr_diamond, reference_arg_amr_count[arg_amr_diamond])
                            
                            if dom_amr_diamond not in domain_vertices.keys():
                                domain_vertices[dom_amr_diamond] = DomainVertex(
                                    dom_amr_diamond, 
                                    reference_dom_amr_count[dom_amr_diamond], 
                                    is_alignment_domain_less)
                            
                            if super_amr_diamond not in super_vertices.keys():
                                super_vertices[super_amr_diamond] = SuperVertex(
                                    super_amr_diamond, 
                                    reference_super_amr_count[super_amr_diamond], 
                                    is_alignment_domain_less)

                            if not query.are_diamond_and_deeparg_the_same():
                                (clstr_deeparg, arg_deeparg, dom_deeparg, super_deeparg, amr_deeparg) = (
                                    query.get_top_deeparg_hit().get_groupings())
                                is_alignment_domain_less = len(query.get_top_deeparg_hit().get_domains()) == 0
                                
                                all_groups_deeparg = '|'.join(query.get_top_deeparg_alignment().get_groupings())
                                clstr_amr_deeparg = '|'.join([clstr_deeparg, amr_deeparg])
                                arg_amr_deeparg = '|'.join([arg_deeparg, amr_deeparg])
                                dom_amr_deeparg = '|'.join([dom_deeparg, amr_deeparg])
                                super_amr_deeparg = '|'.join([super_deeparg, amr_deeparg])

                                if all_groups_deeparg not in all_groups_vertices.keys():
                                    all_groups_vertices[all_groups_deeparg] = MultiVertex(
                                        all_groups_deeparg, 
                                        reference_all_groups_count[all_groups_deeparg], 
                                        is_alignment_domain_less)
                                all_groups_vertices[all_groups_deeparg].add_edge_to_b(
                                    all_groups_vertices[all_groups_diamond], 
                                    all_groups_vertices[all_groups_diamond].get_edge_from_a(
                                        all_groups_vertices[all_groups_deeparg]))

                                if amr_deeparg not in amr_vertices.keys():
                                    amr_vertices[amr_deeparg] = AmrVertex(
                                        amr_deeparg, reference_amr_count[amr_deeparg])
                                amr_vertices[amr_deeparg].add_edge_to_b(
                                    amr_vertices[amr_diamond], 
                                    amr_vertices[amr_diamond].get_edge_from_a(
                                        amr_vertices[amr_deeparg]))

                                if clstr_amr_deeparg not in clstr_vertices.keys():
                                    clstr_vertices[clstr_amr_deeparg] = ClstrVertex(
                                        clstr_amr_deeparg, reference_clstr_amr_count[clstr_amr_deeparg])
                                clstr_vertices[clstr_amr_deeparg].add_edge_to_b(
                                    clstr_vertices[clstr_amr_diamond], 
                                    clstr_vertices[clstr_amr_diamond].get_edge_from_a(
                                        clstr_vertices[clstr_amr_deeparg]))

                                if arg_amr_deeparg not in arg_vertices.keys():
                                    arg_vertices[arg_amr_deeparg] = ArgVertex(
                                        arg_amr_deeparg, reference_arg_amr_count[arg_amr_deeparg])
                                arg_vertices[arg_amr_deeparg].add_edge_to_b(
                                    arg_vertices[arg_amr_diamond], 
                                    arg_vertices[arg_amr_diamond].get_edge_from_a(
                                        arg_vertices[arg_amr_deeparg]))

                                if dom_amr_deeparg not in domain_vertices.keys():
                                    domain_vertices[dom_amr_deeparg] = DomainVertex(
                                        dom_amr_deeparg, 
                                        reference_dom_amr_count[dom_amr_deeparg], 
                                        is_alignment_domain_less)
                                domain_vertices[dom_amr_deeparg].add_edge_to_b(
                                    domain_vertices[dom_amr_diamond], 
                                    domain_vertices[dom_amr_diamond].get_edge_from_a(
                                        domain_vertices[dom_amr_deeparg]))

                                if super_amr_deeparg not in super_vertices.keys():
                                    super_vertices[super_amr_deeparg] = SuperVertex(
                                        super_amr_deeparg, 
                                        reference_super_amr_count[super_amr_deeparg], 
                                        is_alignment_domain_less)
                                super_vertices[super_amr_deeparg].add_edge_to_b(
                                    super_vertices[super_amr_diamond], 
                                    super_vertices[super_amr_diamond].get_edge_from_a(
                                        super_vertices[super_amr_deeparg]))
                            
                            else:
                                amr_vertices[amr_diamond].increment_states_but_not_pair()
                                arg_vertices[arg_amr_diamond].increment_states_but_not_pair()
                                clstr_vertices[clstr_amr_diamond].increment_states_but_not_pair()
                                domain_vertices[dom_amr_diamond].increment_states_but_not_pair()
                                super_vertices[super_amr_diamond].increment_states_but_not_pair()
                                all_groups_vertices[all_groups_diamond].increment_states_but_not_pair()


                    # Create graphs then build switch pair and connected subgraph differential distribution tables
                    amr_graph = Graph(amr_vertices)
                    try:
                        # We only need switch pair differential distribution table for amr
                        amr_pair_abundance = amr_graph.get_pair_table(sample_id, identity, model)
                        if first:
                            amr_pair_abundance.to_csv(
                                path_or_buf=f"{OUTPUT_LOC}/switch_pair_abundance_{AMR_OUTPUT}",
                                sep="\t", float_format='{:.4f}'.format, index=False)
                        else:
                            amr_pair_abundance.to_csv(
                                path_or_buf=f"{OUTPUT_LOC}/switch_pair_abundance_{AMR_OUTPUT}",
                                sep="\t", float_format='{:.4f}'.format, mode='a', header=False, index=False)
                            
                        arg_graph = Graph(arg_vertices)
                        arg_pair_abundance = arg_graph.get_pair_table(sample_id, identity, model)
                        if first:
                            arg_pair_abundance.to_csv(
                                path_or_buf=f"{OUTPUT_LOC}/switch_pair_abundance_{ARG_OUTPUT}",
                                sep="\t", float_format='{:.4f}'.format, index=False)
                        else:
                            arg_pair_abundance.to_csv(
                                path_or_buf=f"{OUTPUT_LOC}/switch_pair_abundance_{ARG_OUTPUT}",
                                sep="\t", float_format='{:.4f}'.format, mode='a', header=False, index=False)
                            
                        all_groups_graph = Graph(all_groups_vertices)
                        all_groups_pair_abundance = all_groups_graph.get_pair_table(sample_id, identity, model)
                        if first:
                            all_groups_pair_abundance.to_csv(
                                path_or_buf=f"{OUTPUT_LOC}/switch_pair_abundance.tsv",
                                sep="\t", float_format='{:.4f}'.format, index=False)
                        else:
                            all_groups_pair_abundance.to_csv(
                                path_or_buf=f"{OUTPUT_LOC}/switch_pair_abundance.tsv",
                                sep="\t", float_format='{:.4f}'.format, mode='a', header=False, index=False)

                        clstr_graph = Graph(clstr_vertices)
                        clstr_connected_abundance = clstr_graph.get_connected_table(sample_id, identity, model)
                        clstr_pair_abundance = clstr_graph.get_pair_table(sample_id, identity, model)
                        if first:
                            clstr_connected_abundance.to_csv(
                                path_or_buf=f"{OUTPUT_LOC}/connected_subgraph_abundance_{CLUSTER_OUTPUT}",
                                sep="\t", float_format='{:.4f}'.format, index=False)
                            clstr_pair_abundance.to_csv(
                                path_or_buf=f"{OUTPUT_LOC}/switch_pair_abundance_{CLUSTER_OUTPUT}",
                                sep="\t", float_format='{:.4f}'.format, index=False)
                        else:
                            clstr_connected_abundance.to_csv(
                                path_or_buf=f"{OUTPUT_LOC}/connected_subgraph_abundance_{CLUSTER_OUTPUT}",
                                sep="\t", float_format='{:.4f}'.format, mode='a', header=False, index=False)
                            clstr_pair_abundance.to_csv(
                                path_or_buf=f"{OUTPUT_LOC}/switch_pair_abundance_{CLUSTER_OUTPUT}",
                                sep="\t", float_format='{:.4f}'.format, mode='a', header=False, index=False)

                        domain_graph = Graph(domain_vertices)
                        domain_connected_abundance = domain_graph.get_connected_table(sample_id, identity, model)
                        domain_pair_abundance = domain_graph.get_pair_table(sample_id, identity, model)
                        if first:
                            domain_connected_abundance.to_csv(
                                path_or_buf=f"{OUTPUT_LOC}/connected_subgraph_abundance_{DOMAIN_OUTPUT}",
                                sep="\t", float_format='{:.4f}'.format, index=False)
                            domain_pair_abundance.to_csv(
                                path_or_buf=f"{OUTPUT_LOC}/switch_pair_abundance_{DOMAIN_OUTPUT}",
                                sep="\t", float_format='{:.4f}'.format, index=False)
                        else:
                            domain_connected_abundance.to_csv(
                                path_or_buf=f"{OUTPUT_LOC}/connected_subgraph_abundance_{DOMAIN_OUTPUT}",
                                sep="\t", float_format='{:.4f}'.format, mode='a', header=False, index=False)
                            domain_pair_abundance.to_csv(
                                path_or_buf=f"{OUTPUT_LOC}/switch_pair_abundance_{DOMAIN_OUTPUT}",
                                sep="\t", float_format='{:.4f}'.format, mode='a', header=False, index=False)

                        super_graph = Graph(super_vertices)
                        super_connected_abundance = super_graph.get_connected_table(sample_id, identity, model)
                        super_pair_abundance = super_graph.get_pair_table(sample_id, identity, model)
                        if first:
                            super_connected_abundance.to_csv(
                                path_or_buf=f"{OUTPUT_LOC}/connected_subgraph_abundance_{SUPER_OUTPUT}",
                                sep="\t",  float_format='{:.4f}'.format, index=False)
                            super_pair_abundance.to_csv(
                                path_or_buf=f"{OUTPUT_LOC}/switch_pair_abundance_{SUPER_OUTPUT}",
                                sep="\t", float_format='{:.4f}'.format, index=False)
                        else:
                            super_connected_abundance.to_csv(
                                path_or_buf=f"{OUTPUT_LOC}/connected_subgraph_abundance_{SUPER_OUTPUT}",
                                sep="\t",  float_format='{:.4f}'.format, mode='a', header=False, index=False)
                            super_pair_abundance.to_csv(
                                path_or_buf=f"{OUTPUT_LOC}/switch_pair_abundance_{SUPER_OUTPUT}",
                                sep="\t", float_format='{:.4f}'.format, mode='a', header=False, index=False)
                        
                        first = False
                    except:
                        continue

main()