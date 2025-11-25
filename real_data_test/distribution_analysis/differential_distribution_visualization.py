from differential_distribution_classes.graph_items import TrioVertex, ClstrVertex, DomainVertex, SuperVertex
from differential_distribution_classes.reference import Reference
from differential_distribution_classes.graph import Graph
from differential_distribution_classes.query import Query
import pandas as pd
import numpy as np
import os

CDD_DIR = "../../CDD_features/"
REF_LOC = "../../data/database/v2/features.fasta"
CLSTR_LOC = "../database_clustering/db40.clstr"
DEEPARG_HIT_FILE = "X.mapping.ARG"
ALIGNMENT_FILE = "X.align.daa.tsv"
SAMPLE_ID_FILE = "../real_samples.txt"
PATH_TO_SS = "deeparg_results"
PATH_TO_LS = "spades/deeparg_results"

AMR_SWITCH_DISTRIBUTION = "amr_switch_distribution.txt"
CLR_LOC = "center_log_ratio_transform"
CLUSTER_OUTPUT = "cluster.tsv"
TRIO_OUTPUT = "trio.tsv"
DOMAIN_OUTPUT = "domain.tsv"
SUPER_OUTPUT = "super.tsv"

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
    
    # Get amr class distribution in database
    database_amr_count: dict[str, int] = dict()

    # Look through CD features and intialize Reference object
    ref_dict: dict[str, Reference] = dict()
    for part in range(1,26):
        with open(f"{CDD_DIR}Part{part}_hitdata.txt", "r") as cdd_info_buf:
            cdd_info_list = cdd_info_buf.read().split('\n')
        for cdd_info_line in cdd_info_list[1:-1]:
            row = cdd_info_line.split('\t')
            name = row[0].split('>')[-1]
            if (len(ref_dict) > 0) and (name in ref_dict.keys()):
                ref_dict[name].add_domain(row)
            else:
                ref_dict[name] = Reference(row)
                amr = ref_dict[name].get_classification()
                if amr in database_amr_count.keys():
                    database_amr_count[amr] += 1
                else:
                    database_amr_count[amr] = 1
        cdd_info_list.clear()

    # Add additional references that don't have domains + save length information
    with open(REF_LOC, "r") as ref_buf:
        ref_list = ref_buf.read().split('\n')
    for ref_index in range((len(ref_list) - 1) // 2):
        name = ref_list[ref_index*2][1:]
        if name not in ref_dict.keys():
            ref_dict[name] = Reference(name)
            amr = ref_dict[name].get_classification()
            if amr in database_amr_count.keys():
                database_amr_count[amr] += 1
            else:
                database_amr_count[amr] = 0
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
            name = line.split('>')[-1].split('...')[0]
            ref_dict[name].define_cluster(cluster)

    # Count the time each id appears in reference
    reference_clstr_id_count: dict[str, int] = dict()
    reference_trio_id_count: dict[str, int] = dict()
    reference_dom_id_count: dict[str, int] = dict()
    reference_super_id_count: dict[str, int] = dict()
    for ref in ref_dict.values():
        ref_ids = ref.get_domain_identifiers()
        (clstr_id, trio_id, dom_id, super_id) = ref_ids
        if clstr_id in reference_clstr_id_count.keys():
            reference_clstr_id_count[clstr_id] += 1
        else:
            reference_clstr_id_count[clstr_id] = 1

        trio_fields = trio_id.split('|')
        super_fields = super_id.split('|')
        dom_list = trio_fields[0].split('$')
        super_list = super_fields[0].split('$')
        # We may have numerous domains in gene; should count each potential combination as well
        dom_combo_list = sorted_domains_combo(dom_list)
        super_combo_list = sorted_domains_combo(super_list)
        # We also need to count instances of arg as domain
        if trio_fields[0] != trio_fields[1]:
            dom_combo_list.append(trio_fields[1])
            super_combo_list.append(trio_fields[1])
        for i in range(len(dom_combo_list)):
            trio_id = "|".join((dom_combo_list[i], trio_fields[1], trio_fields[2]))
            if trio_id in reference_trio_id_count.keys():
                reference_trio_id_count[trio_id] += 1
            else:
                reference_trio_id_count[trio_id] = 1
            dom_id = "|".join((dom_combo_list[i], trio_fields[2]))
            if dom_id in reference_dom_id_count.keys():
                reference_dom_id_count[dom_id] += 1
            else:
                reference_dom_id_count[dom_id] = 1
            super_id = "|".join((super_combo_list[i], super_fields[1]))
            if super_id in reference_super_id_count.keys():
                reference_super_id_count[super_id] += 1
            else:
                reference_super_id_count[super_id] = 1
        
    # Create list of AMR switch distributions to output later
    amr_switch_info = list()

    # Check to see if center log ratio transform directory exist, or make one
    if not os.path.exists(CLR_LOC):
        os.mkdir(CLR_LOC)
    
    # Go through each run one at a time
    for sample_id in sample_id_list:
        for identity in [30, 50, 80]:
            for model in ["LS", "SS"]:
                # Make vertices for all four ARG categorizing units
                trio_vertices: dict[str, TrioVertex] = dict()
                clstr_vertices: dict[str, ClstrVertex] = dict()
                domain_vertices: dict[str, DomainVertex] = dict()
                super_vertices: dict[str, SuperVertex] = dict()
                
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
                        query_dict[row[0]] = Query(row, ref_dict, model == "LS")
                alignment_list.clear()

                # Create dict for amr switches
                amr_switch_dict: dict[str, int] = dict()

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
                        if key in amr_switch_dict.keys():
                            amr_switch_dict[key] += 1
                        else:
                            amr_switch_dict[key] = 1
                deeparg_list.clear()

                # Save AMR switch distribution for current run
                for (dia_to_dee, count) in amr_switch_dict.items():
                    amr_switch_info.append("\t".join((
                        sample_id, str(identity), model, dia_to_dee, str(count))))
                
                # Create vertices nd edges to represent categories and pairs, respectively
                for query in query_dict.values():
                    if query.is_deeparg_hit():
                        alignment_a_ids = query.get_top_diamond_alignment_domain_identifiers()
                        (clstr_a, trio_a, dom_a, super_a) = alignment_a_ids

                        if clstr_a not in clstr_vertices.keys():
                            clstr_vertices[clstr_a] = ClstrVertex(clstr_a, reference_clstr_id_count[clstr_a])
                        if trio_a not in trio_vertices.keys():
                            trio_vertices[trio_a] = TrioVertex(trio_a, reference_trio_id_count[trio_a])
                        if dom_a not in domain_vertices.keys():
                            domain_vertices[dom_a] = DomainVertex(dom_a, reference_dom_id_count[dom_a])
                        if super_a not in super_vertices.keys():
                            super_vertices[super_a] = SuperVertex(super_a, reference_super_id_count[super_a])

                        if not query.are_diamond_and_deeparg_the_same():
                            alignment_b_ids = query.get_top_deeparg_hit_domain_identifiers()
                            (clstr_b, trio_b, dom_b, super_b) = alignment_b_ids

                            if clstr_b not in clstr_vertices.keys():
                                clstr_vertices[clstr_b] = ClstrVertex(clstr_b, reference_clstr_id_count[clstr_b])
                            clstr_edge = clstr_vertices[clstr_a].get_edge_from_a(clstr_vertices[clstr_b])
                            clstr_vertices[clstr_b].add_edge_to_b(clstr_vertices[clstr_a], clstr_edge)
                        
                            if trio_b not in trio_vertices.keys():
                                trio_vertices[trio_b] = TrioVertex(trio_b, reference_trio_id_count[trio_b])
                            trio_edge = trio_vertices[trio_a].get_edge_from_a(trio_vertices[trio_b])
                            trio_vertices[trio_b].add_edge_to_b(trio_vertices[trio_a], trio_edge)

                            if dom_b not in domain_vertices.keys():
                                domain_vertices[dom_b] = DomainVertex(dom_b, reference_dom_id_count[dom_b])
                            dom_edge = domain_vertices[dom_a].get_edge_from_a(domain_vertices[dom_b])
                            domain_vertices[dom_b].add_edge_to_b(domain_vertices[dom_a], dom_edge)

                            if super_b not in super_vertices.keys():
                                super_vertices[super_b] = SuperVertex(super_b, reference_super_id_count[super_b])
                            super_edge = super_vertices[super_a].get_edge_from_a(super_vertices[super_b])
                            super_vertices[super_b].add_edge_to_b(super_vertices[super_a], super_edge)
                        
                        else:
                            clstr_vertices[clstr_a].increment_states_but_not_pair()
                            trio_vertices[trio_a].increment_states_but_not_pair()
                            domain_vertices[dom_a].increment_states_but_not_pair()
                            super_vertices[super_a].increment_states_but_not_pair()

                # Create graphs then build center log ratio transform tables
                clstr_graph = Graph(clstr_vertices)
                clstr_connected_clr = clstr_graph.get_connected_clr_transform()
                clstr_connected_clr.to_csv(
                    path_or_buf=f"{CLR_LOC}/{sample_id}_{identity}_{model}_connected_{CLUSTER_OUTPUT}",
                    sep="\t", index_label="cluster|amr", float_format='{:.4f}'.format)
                clstr_relative_abundance_clr = clstr_graph.get_relative_abundance()
                clstr_relative_abundance_clr.to_csv(
                    path_or_buf=f"{CLR_LOC}/{sample_id}_{identity}_{model}_relative_abundance_{CLUSTER_OUTPUT}",
                    sep="\t", index_label="cluster|amr", float_format='{:.4f}'.format)
                clstr_adjacent_clr = clstr_graph.get_adjacent_clr_transform()
                clstr_adjacent_clr.to_csv(
                    path_or_buf=f"{CLR_LOC}/{sample_id}_{identity}_{model}_adjacent_{CLUSTER_OUTPUT}",
                    sep="\t", index_label="cluster|amr", float_format='{:.4f}'.format)

                trio_graph = Graph(trio_vertices)
                trio_connected_clr = trio_graph.get_connected_clr_transform()
                trio_connected_clr.to_csv(
                    path_or_buf=f"{CLR_LOC}/{sample_id}_{identity}_{model}_connected_{TRIO_OUTPUT}",
                    sep="\t", index_label="domain|arg|amr", float_format='{:.4f}'.format)
                trio_relative_abundance_clr = trio_graph.get_relative_abundance()
                trio_relative_abundance_clr.to_csv(
                    path_or_buf=f"{CLR_LOC}/{sample_id}_{identity}_{model}_relative_abundance_{TRIO_OUTPUT}",
                    sep="\t", index_label="domain|arg|amr", float_format='{:.4f}'.format)
                trio_adjacent_clr = trio_graph.get_adjacent_clr_transform()
                trio_adjacent_clr.to_csv(
                    path_or_buf=f"{CLR_LOC}/{sample_id}_{identity}_{model}_adjacent_{TRIO_OUTPUT}",
                    sep="\t", index_label="domain|arg|amr", float_format='{:.4f}'.format)

                domain_graph = Graph(domain_vertices)
                domain_connected_clr = domain_graph.get_connected_clr_transform()
                domain_connected_clr.to_csv(
                    path_or_buf=f"{CLR_LOC}/{sample_id}_{identity}_{model}_connected_{DOMAIN_OUTPUT}",
                    sep="\t", index_label="domain|amr", float_format='{:.4f}'.format)
                domain_relative_abundance_clr = domain_graph.get_relative_abundance()
                domain_relative_abundance_clr.to_csv(
                    path_or_buf=f"{CLR_LOC}/{sample_id}_{identity}_{model}_relative_abundance_{DOMAIN_OUTPUT}",
                    sep="\t", index_label="domain|amr", float_format='{:.4f}'.format)
                domain_adjacent_clr = domain_graph.get_adjacent_clr_transform()
                domain_adjacent_clr.to_csv(
                    path_or_buf=f"{CLR_LOC}/{sample_id}_{identity}_{model}_adjacent_{DOMAIN_OUTPUT}",
                    sep="\t", index_label="domain|amr", float_format='{:.4f}'.format)

                super_graph = Graph(super_vertices)
                super_connected_clr = super_graph.get_connected_clr_transform()
                super_connected_clr.to_csv(
                    path_or_buf=f"{CLR_LOC}/{sample_id}_{identity}_{model}_connected_{SUPER_OUTPUT}",
                    sep="\t", index_label="super|amr", float_format='{:.4f}'.format)
                super_relative_abundance_clr = super_graph.get_relative_abundance()
                super_relative_abundance_clr.to_csv(
                    path_or_buf=f"{CLR_LOC}/{sample_id}_{identity}_{model}_relative_abundance_{SUPER_OUTPUT}",
                    sep="\t", index_label="super|amr", float_format='{:.4f}'.format)
                super_adjacent_clr = super_graph.get_adjacent_clr_transform()
                super_adjacent_clr.to_csv(
                    path_or_buf=f"{CLR_LOC}/{sample_id}_{identity}_{model}_adjacent_{SUPER_OUTPUT}",
                    sep="\t", index_label="super|amr", float_format='{:.4f}'.format)

                # And now let's output three mega tables where each combos of cluster and 
                # trio or domain|amr or super|amr is connected to the appropriate state_i_clr
                # a_to_b_sign
                all_combo_tuples: set[tuple[str,str,str,str]] = set()
                for query in query_dict.values():
                    if query.is_deeparg_hit():
                        alignment_a_ids = query.get_top_diamond_alignment_domain_identifiers()
                        all_combo_tuples.add(alignment_a_ids)
                        if not query.are_diamond_and_deeparg_the_same():
                            alignment_b_ids = query.get_top_deeparg_hit_domain_identifiers()
                            all_combo_tuples.add(alignment_b_ids)
                all_combo_indices = pd.MultiIndex.from_tuples(all_combo_tuples, names=[
                    "cluster|amr", "trio", "domain|amr", "super|amr"])
                all_combo_table = pd.DataFrame(
                    np.empty(shape=[len(all_combo_tuples), 28]),
                    index=all_combo_indices, 
                    columns=[
                        "cluster|amr state I connected clr", "cluster|amr state I adjacent clr", 
                        "cluster|amr state A connected clr", "cluster|amr state A adjacent clr", 
                        "cluster|amr state B connected clr", "cluster|amr state B adjacent clr", 
                        "cluster|amr state A to B sign", "trio state I connected clr", 
                        "trio state I adjacent clr", "trio state A connected clr", 
                        "trio state A adjacent clr", "trio state B connected clr", 
                        "trio state B adjacent clr", "trio state A to B sign",
                        "domain|amr state I connected clr", "domain|amr state I adjacent clr", 
                        "domain|amr state A connected clr", "domain|amr state A adjacent clr",
                        "domain|amr state B connected clr", "domain|amr state B adjacent clr",
                        "domain|amr state A to B sign", "super|amr state I connected clr", 
                        "super|amr state I adjacent clr", "super|amr state A connected clr", 
                        "super|amr state A adjacent clr", "super|amr state B connected clr", 
                        "super|amr state B adjacent clr", "super|amr state A to B sign"])
                
                all_combo_table["cluster|amr state I connected clr"] = [
                    clstr_connected_clr.at[x, "state I clr"]
                    for x in list(all_combo_indices.get_level_values("cluster|amr").values)]
                all_combo_table["cluster|amr state I adjacent clr"] = [
                    clstr_adjacent_clr.at[x, "state I clr"] 
                    for x in list(all_combo_indices.get_level_values("cluster|amr").values)]
                all_combo_table["cluster|amr state A connected clr"] = [
                    clstr_connected_clr.at[x, "state A clr"]
                    for x in list(all_combo_indices.get_level_values("cluster|amr").values)]
                all_combo_table["cluster|amr state A adjacent clr"] = [
                    clstr_adjacent_clr.at[x, "state A clr"] 
                    for x in list(all_combo_indices.get_level_values("cluster|amr").values)]
                all_combo_table["cluster|amr state B connected clr"] = [
                    clstr_connected_clr.at[x, "state B clr"]
                    for x in list(all_combo_indices.get_level_values("cluster|amr").values)]
                all_combo_table["cluster|amr state B adjacent clr"] = [
                    clstr_adjacent_clr.at[x, "state B clr"] 
                    for x in list(all_combo_indices.get_level_values("cluster|amr").values)]
                all_combo_table["cluster|amr state A to B sign"] = [
                    clstr_adjacent_clr.at[x, "state B - state A sign"] 
                    for x in list(all_combo_indices.get_level_values("cluster|amr").values)]
                
                all_combo_table["trio state I connected clr"] = [
                    trio_connected_clr.at[x, "state I clr"]
                    for x in list(all_combo_indices.get_level_values("trio").values)]
                all_combo_table["trio state I adjacent clr"] = [
                    trio_adjacent_clr.at[x, "state I clr"] 
                    for x in list(all_combo_indices.get_level_values("trio").values)]
                all_combo_table["trio state A connected clr"] = [
                    trio_connected_clr.at[x, "state A clr"]
                    for x in list(all_combo_indices.get_level_values("trio").values)]
                all_combo_table["trio state A adjacent clr"] = [
                    trio_adjacent_clr.at[x, "state A clr"] 
                    for x in list(all_combo_indices.get_level_values("trio").values)]
                all_combo_table["trio state B connected clr"] = [
                    trio_connected_clr.at[x, "state B clr"]
                    for x in list(all_combo_indices.get_level_values("trio").values)]
                all_combo_table["trio state B adjacent clr"] = [
                    trio_adjacent_clr.at[x, "state B clr"] 
                    for x in list(all_combo_indices.get_level_values("trio").values)]
                all_combo_table["trio state A to B sign"] = [
                    trio_adjacent_clr.at[x, "state B - state A sign"] 
                    for x in list(all_combo_indices.get_level_values("trio").values)]
                
                all_combo_table["domain|amr state I connected clr"] = [
                    domain_connected_clr.at[x, "state I clr"]
                    for x in list(all_combo_indices.get_level_values("domain|amr").values)]
                all_combo_table["domain|amr state I adjacent clr"] = [
                    domain_adjacent_clr.at[x, "state I clr"] 
                    for x in list(all_combo_indices.get_level_values("domain|amr").values)]
                all_combo_table["domain|amr state A connected clr"] = [
                    domain_connected_clr.at[x, "state A clr"]
                    for x in list(all_combo_indices.get_level_values("domain|amr").values)]
                all_combo_table["domain|amr state A adjacent clr"] = [
                    domain_adjacent_clr.at[x, "state A clr"] 
                    for x in list(all_combo_indices.get_level_values("domain|amr").values)]
                all_combo_table["domain|amr state B connected clr"] = [
                    domain_connected_clr.at[x, "state B clr"]
                    for x in list(all_combo_indices.get_level_values("domain|amr").values)]
                all_combo_table["domain|amr state B adjacent clr"] = [
                    domain_adjacent_clr.at[x, "state B clr"] 
                    for x in list(all_combo_indices.get_level_values("domain|amr").values)]
                all_combo_table["domain|amr state A to B sign"] = [
                    domain_adjacent_clr.at[x, "state B - state A sign"] 
                    for x in list(all_combo_indices.get_level_values("domain|amr").values)]
                
                all_combo_table["super|amr state I connected clr"] = [
                    super_connected_clr.at[x, "state I clr"]
                    for x in list(all_combo_indices.get_level_values("super|amr").values)]
                all_combo_table["super|amr state I adjacent clr"] = [
                    super_adjacent_clr.at[x, "state I clr"] 
                    for x in list(all_combo_indices.get_level_values("super|amr").values)]
                all_combo_table["super|amr state A connected clr"] = [
                    super_connected_clr.at[x, "state A clr"]
                    for x in list(all_combo_indices.get_level_values("super|amr").values)]
                all_combo_table["super|amr state A adjacent clr"] = [
                    super_adjacent_clr.at[x, "state A clr"] 
                    for x in list(all_combo_indices.get_level_values("super|amr").values)]
                all_combo_table["super|amr state B connected clr"] = [
                    super_connected_clr.at[x, "state B clr"]
                    for x in list(all_combo_indices.get_level_values("super|amr").values)]
                all_combo_table["super|amr state B adjacent clr"] = [
                    super_adjacent_clr.at[x, "state B clr"] 
                    for x in list(all_combo_indices.get_level_values("super|amr").values)]
                all_combo_table["super|amr state A to B sign"] = [
                    super_adjacent_clr.at[x, "state B - state A sign"] 
                    for x in list(all_combo_indices.get_level_values("super|amr").values)]
                
                trio_combo_table = all_combo_table.droplevel(["domain|amr","super|amr"])
                domain_combo_table = all_combo_table[[
                    "cluster|amr state I connected clr", "cluster|amr state I adjacent clr",
                    "cluster|amr state A connected clr", "cluster|amr state A adjacent clr",
                    "cluster|amr state A connected clr", "cluster|amr state A adjacent clr",
                    "cluster|amr state A to B sign", "domain|amr state I connected clr", 
                    "domain|amr state I adjacent clr", "domain|amr state A connected clr", 
                    "domain|amr state A adjacent clr", "domain|amr state B connected clr", 
                    "domain|amr state B adjacent clr", "domain|amr state A to B sign", 
                    "super|amr state I connected clr", "super|amr state I adjacent clr",
                    "super|amr state A connected clr", "super|amr state A adjacent clr",
                    "super|amr state B connected clr", "super|amr state B adjacent clr", 
                    "super|amr state A to B sign"]]
                domain_combo_table = domain_combo_table.droplevel(["trio","super|amr"])
                domain_combo_table = domain_combo_table.loc[
                    ~domain_combo_table.index.duplicated(keep='first'), :]
                super_combo_table = all_combo_table[[
                    "cluster|amr state I connected clr", "cluster|amr state I adjacent clr", 
                        "cluster|amr state A connected clr", "cluster|amr state A adjacent clr", 
                        "cluster|amr state B connected clr", "cluster|amr state B adjacent clr", 
                        "cluster|amr state A to B sign", "super|amr state I connected clr", 
                        "super|amr state I adjacent clr", "super|amr state A connected clr", 
                        "super|amr state A adjacent clr", "super|amr state B connected clr", 
                        "super|amr state B adjacent clr", "super|amr state A to B sign"]]
                super_combo_table = super_combo_table.droplevel(["trio","domain|amr"])
                super_combo_table = super_combo_table.loc[
                    ~super_combo_table.index.duplicated(keep='first'), :]
                
                trio_combo_table.to_csv(
                    path_or_buf=f"{CLR_LOC}/{sample_id}_{identity}_{model}_combo_cluster_{TRIO_OUTPUT}",
                    sep="\t", float_format='{:.4f}'.format)
                domain_combo_table.to_csv(
                    path_or_buf=f"{CLR_LOC}/{sample_id}_{identity}_{model}_combo_cluster_{DOMAIN_OUTPUT}",
                    sep="\t", float_format='{:.4f}'.format)
                super_combo_table.to_csv(
                    path_or_buf=f"{CLR_LOC}/{sample_id}_{identity}_{model}_combo_cluster_{SUPER_OUTPUT}",
                    sep="\t", float_format='{:.4f}'.format)


    # Output AMR switch distribution information
    with open(AMR_SWITCH_DISTRIBUTION, "w") as amr_switch_buf:
        amr_switch_buf.write("AMR classification\tCount\n")
        for (amr_class, count) in database_amr_count.items():
            amr_switch_buf.write(f"{amr_class}\t{count}\n")
        amr_switch_buf.write(
            "Sample\tAlignment Identity\tModel\tDiamond AMR\tDeepARG AMR\tCount\n")
        for row in amr_switch_info:
            amr_switch_buf.write(f"{row}\n")

main()