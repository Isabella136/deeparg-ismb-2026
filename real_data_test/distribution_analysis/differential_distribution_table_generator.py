from differential_distribution_classes.graph_items import ClstrVertex, DomainVertex, SuperVertex, AmrVertex
from differential_distribution_classes.reference import Reference
from differential_distribution_classes.graph import Graph
from differential_distribution_classes.query import Query
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

AMR_SWITCH_DISTRIBUTION = "amr_switch_distribution.txt"
OUTPUT_LOC = "differential_distribution_output"
CLUSTER_OUTPUT = "cluster.tsv"
TRIO_OUTPUT = "trio.tsv"
AMR_OUTPUT = "amr.tsv"
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

        # Count the time each id appears in reference
        reference_clstr_id_count: dict[str, int] = dict()
        reference_dom_id_count: dict[str, int] = dict()
        reference_super_id_count: dict[str, int] = dict()
        reference_amr_count: dict[str, int] = dict()
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
            if trio_fields[2] in reference_amr_count.keys():
                reference_amr_count[trio_fields[2]] += 1
            else:
                reference_amr_count[trio_fields[2]] = 1
            for i in range(len(dom_combo_list)):
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

        # Check to see if output directory exist, or make one
        if not os.path.exists(OUTPUT_LOC):
            os.mkdir(OUTPUT_LOC)
        
        # Go through each run one at a time
        for sample_id in sample_id_list:
            for identity in [30, 50, 80]:
                # Make vertices for all four labels
                amr_vertices: dict[str, AmrVertex] = dict()
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
                
                # Create vertices and edges to represent labels and pairs, respectively
                for query in query_dict.values():
                    if query.is_deeparg_hit():
                        diamond_ids = query.get_top_diamond_alignment_domain_identifiers()
                        (clstr_dia, trio_dia, dom_dia, super_dia) = diamond_ids
                        is_alignment_domain_less = len(query.get_top_diamond_alignment().get_domains()) == 0
                        amr_dia = query.get_top_diamond_classification()

                        if amr_dia not in amr_vertices.keys():
                            amr_vertices[amr_dia] = AmrVertex(amr_dia, reference_amr_count[amr_dia])
                        if clstr_dia not in clstr_vertices.keys():
                            clstr_vertices[clstr_dia] = ClstrVertex(clstr_dia, reference_clstr_id_count[clstr_dia])
                        if dom_dia not in domain_vertices.keys():
                            domain_vertices[dom_dia] = DomainVertex(
                                dom_dia, reference_dom_id_count[dom_dia], is_alignment_domain_less)
                        if super_dia not in super_vertices.keys():
                            super_vertices[super_dia] = SuperVertex(
                                super_dia, reference_super_id_count[super_dia], is_alignment_domain_less)

                        if not query.are_diamond_and_deeparg_the_same():
                            deeparg_ids = query.get_top_deeparg_hit_domain_identifiers()
                            (clstr_dee, trio_dee, dom_dee, super_dee) = deeparg_ids
                            is_alignment_domain_less = len(query.get_top_deeparg_hit().get_domains()) == 0
                            amr_dee = query.get_top_deeparg_classification()

                            if amr_dee not in amr_vertices.keys():
                                amr_vertices[amr_dee] = AmrVertex(amr_dee, reference_amr_count[amr_dee])
                            amr_edge = amr_vertices[amr_dia].get_edge_from_a(amr_vertices[amr_dee])
                            amr_vertices[amr_dee].add_edge_to_b(amr_vertices[amr_dia], amr_edge)

                            if clstr_dee not in clstr_vertices.keys():
                                clstr_vertices[clstr_dee] = ClstrVertex(clstr_dee, reference_clstr_id_count[clstr_dee])
                            clstr_edge = clstr_vertices[clstr_dia].get_edge_from_a(clstr_vertices[clstr_dee])
                            clstr_vertices[clstr_dee].add_edge_to_b(clstr_vertices[clstr_dia], clstr_edge)

                            if dom_dee not in domain_vertices.keys():
                                domain_vertices[dom_dee] = DomainVertex(
                                    dom_dee, reference_dom_id_count[dom_dee], is_alignment_domain_less)
                            dom_edge = domain_vertices[dom_dia].get_edge_from_a(domain_vertices[dom_dee])
                            domain_vertices[dom_dee].add_edge_to_b(domain_vertices[dom_dia], dom_edge)

                            if super_dee not in super_vertices.keys():
                                super_vertices[super_dee] = SuperVertex(
                                    super_dee, reference_super_id_count[super_dee], is_alignment_domain_less)
                            super_edge = super_vertices[super_dia].get_edge_from_a(super_vertices[super_dee])
                            super_vertices[super_dee].add_edge_to_b(super_vertices[super_dia], super_edge)
                        
                        else:
                            amr_vertices[amr_dia].increment_states_but_not_pair()
                            clstr_vertices[clstr_dia].increment_states_but_not_pair()
                            domain_vertices[dom_dia].increment_states_but_not_pair()
                            super_vertices[super_dia].increment_states_but_not_pair()


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