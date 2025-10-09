from differential_distribution_classes.graph_items import TrioVertex, ClstrVertex, DomainVertex, SuperVertex
from differential_distribution_classes.query import Query
from differential_distribution_classes.reference import Reference

CDD_DIR = "../../CDD_features/"
REF_LOC = "../../data/database/v2/features.fasta"
CLSTR_LOC = "../database_clustering/db40.clstr"
DEEPARG_HIT_FILE = "X.mapping.ARG"
ALIGNMENT_FILE = "X.align.daa.tsv"
SAMPLE_ID_FILE = "../real_samples.txt"
PATH_TO_SS = "deeparg_results"
PATH_TO_LS = "spades/deeparg_results"

AMR_SWITCH_DISTRIBUTION = "amr_switch_distribution.txt"

def main():
    # Read content of real_samples.txt to find biosample ID
    with open(SAMPLE_ID_FILE, "r") as sample_id_buf:
        sample_id_list = sample_id_buf.read().split('\n')
    
    # Get amr class distribution in database
    database_amr_count: dict[str, int] =  dict()

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
    for line in clstr_list:
        if '>' == line[0]:
            cluster += 1
        else:
            name = line.split('>')[-1].split('.')[0]
            ref_dict[name].define_cluster(cluster)

    # Make vertices for all four ARG categorizing units
    trio_vertices: dict[str, TrioVertex] = dict()
    clstr_vertices: dict[str, ClstrVertex] = dict()
    domain_vertices: dict[str, DomainVertex] = dict()
    super_vertices: dict[str, SuperVertex] = dict()

    amr_switch_info = list()
    # Go through each run one at a time
    for sample_id in sample_id_list:
        for identity in [30, 50, 80]:
            for model in ["LS", "SS"]:
                # Get Diamond alignment information
                query_dict: dict[str, Query] = dict()
                with open("/".join((
                        "..", "samples", sample_id, PATH_TO_LS if model=="LS" else PATH_TO_SS,
                        f"arg_alignment_identity_{identity}", ALIGNMENT_FILE)), "r") as alignment_buf:
                    alignment_list = alignment_buf.read().split('\n')
                for alignment in alignment_list[:-1]:
                    row = alignment.split('\t')
                    if (len(query_dict) > 0) and (row[0] in query_dict.keys()):
                        query_dict[row[0]].add_alignment(row, ref_dict)
                    else:
                        query_dict[row[0]] = Query(row, ref_dict)
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
                for (dia_to_dee, count) in amr_switch_dict.items():
                    amr_switch_info.append("\t".join((
                        sample_id, str(identity), model, dia_to_dee, str(count))))
                    
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