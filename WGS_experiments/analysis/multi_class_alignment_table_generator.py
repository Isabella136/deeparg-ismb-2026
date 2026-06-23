import pathlib
import pickle

import pandas as pd
from classes.query import Query
from classes.reference import Reference

CDD_DIR = "../../database/CDD_features_v2/"
REF_LOC = "../../data/database/v2/features.fasta"
CLSTR_LOC = "../../database/database_clustering/db_v2_40.clstr"
DEEPARG_HIT_FILE = "X.mapping.ARG"
ALIGNMENT_FILE = "X.align.daa.tsv"
SAMPLE_ID_FILE = "../real_samples.txt"
PATH_TO_SS = "deeparg_results"
PATH_TO_LS = "spades/deeparg_results"


# Domains are sorted by bitscore
def sorted_domains_combo(dom_list: list[str]) -> list[str]:  # noqa: D103
    if len(dom_list) == 1:
        return dom_list
    first_dom = dom_list[0]
    recursive_combo_list = sorted_domains_combo(dom_list[1:])
    extended_combo_list = [
        f"{first_dom}${curr_combo}" for curr_combo in recursive_combo_list
    ]
    extended_combo_list.append(first_dom)
    extended_combo_list.extend(recursive_combo_list)
    return extended_combo_list


# Read content of real_samples.txt to find biosample ID
sample_id_list = (
    pathlib.Path(SAMPLE_ID_FILE).read_text(encoding="utf-8").split("\n")
)
first = True
for model in ["LS", "SS"]:
    # Open metadata to get ARG features and output
    ref_dict: dict[str, Reference] = {}
    with pathlib.Path("../../data/model/v2/metadata_" + model + ".pkl").open(
        "rb"
    ) as metadata_pkl:
        metadata = pickle.load(metadata_pkl)
    for ref in metadata["features"]:
        ref_dict[ref] = Reference(ref)

    # Look through CD features to add domain info
    for part in range(1, 26):
        cdd_info_list = (
            pathlib
            .Path(f"{CDD_DIR}Part{part}_hitdata.txt")
            .read_text(encoding="utf-8")
            .split("\n")
        )
        for cdd_info_line in cdd_info_list[1:-1]:
            row = cdd_info_line.split("\t")
            name = row[0].split(">")[-1]
            if name in ref_dict:
                ref_dict[name].add_domain(row)
        cdd_info_list.clear()

    # Add length information
    ref_list = pathlib.Path(REF_LOC).read_text(encoding="utf-8").split("\n")
    for ref_index in range((len(ref_list) - 1) // 2):
        name = ref_list[ref_index * 2][1:]
        if name not in ref_dict:
            continue
        ref_dict[name].add_length_info(len(ref_list[ref_index * 2 + 1]))
    ref_list.clear()

    # Add cluster information
    clstr_list = pathlib.Path(CLSTR_LOC).read_text(encoding="utf-8").split("\n")
    cluster = -1
    for line in clstr_list[:-1]:
        if line[0] == ">":
            cluster += 1
        else:
            name = line.split(">")[-1].split("...")[0][:-3]
            if name not in ref_dict:
                continue
            ref_dict[name].define_cluster(cluster)

    # Count the time each label appears in reference
    reference_clstr_amr_count: dict[str, int] = {}
    reference_arg_amr_count: dict[str, int] = {}
    reference_dom_amr_count: dict[str, int] = {}
    reference_super_amr_count: dict[str, int] = {}
    reference_amr_count: dict[str, int] = {}
    reference_all_groups_count: dict[str, int] = {}
    for ref in ref_dict.values():
        (clstr, arg, dom, superf, amr) = ref.get_annotations()

        # Count amr, clstr|amr, and arg|amr
        clstr_amr = f"{clstr}|{amr}"
        arg_amr = f"{arg}|{amr}"
        if amr in reference_amr_count:
            reference_amr_count[amr] += 1
        else:
            reference_amr_count[amr] = 1
        if clstr_amr in reference_clstr_amr_count:
            reference_clstr_amr_count[clstr_amr] += 1
        else:
            reference_clstr_amr_count[clstr_amr] = 1
        if arg_amr in reference_arg_amr_count:
            reference_arg_amr_count[arg_amr] += 1
        else:
            reference_arg_amr_count[arg_amr] = 1

        # We may have numerous domains in gene; should count each potential
        # combination as well
        dom_combo_list = sorted_domains_combo(dom.split("$"))
        super_combo_list = sorted_domains_combo(superf.split("$"))

        # We also need to count instances of arg as domain
        if dom != arg:
            dom_combo_list.append(arg)
            super_combo_list.append(arg)

        # Count each dom|amr, superf|amr, and clstr|arg|dom|superf|amr
        for i in range(len(dom_combo_list)):
            dom_amr = "|".join([dom_combo_list[i], amr])
            super_amr = "|".join([super_combo_list[i], amr])
            all_groups = "|".join([
                clstr,
                arg,
                dom_combo_list[i],
                super_combo_list[i],
                amr,
            ])

            if dom_amr in reference_dom_amr_count:
                reference_dom_amr_count[dom_amr] += 1
            else:
                reference_dom_amr_count[dom_amr] = 1
            if super_amr in reference_super_amr_count:
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
            query_dict: dict[str, Query] = {}
            alignment_list = (
                pathlib
                .Path(
                    "/".join((
                        "..",
                        "samples",
                        sample_id,
                        PATH_TO_LS if model == "LS" else PATH_TO_SS,
                        f"arg_alignment_identity_{identity}",
                        ALIGNMENT_FILE,
                    ))
                )
                .read_text(encoding="utf-8")
                .split("\n")
            )
            for alignment in alignment_list[:-1]:
                row = alignment.split("\t")
                if (len(query_dict) > 0) and (row[0] in query_dict):
                    query_dict[row[0]].add_alignment(
                        row, ref_dict, model == "LS"
                    )
                else:
                    query_dict[row[0]] = Query(
                        row, ref_dict, model == "LS", sample_id
                    )
            alignment_list.clear()

            # Get DeepARG hit information
            deeparg_list = (
                pathlib
                .Path(
                    "/".join((
                        "..",
                        "samples",
                        sample_id,
                        PATH_TO_LS if model == "LS" else PATH_TO_SS,
                        f"arg_alignment_identity_{identity}",
                        DEEPARG_HIT_FILE,
                    ))
                )
                .read_text(encoding="utf-8")
                .split("\n")
            )
            for hit in deeparg_list[1:-1]:
                row = hit.split("\t")

                # Mark queries that match with read_id (row[3]) as deeparg_hit
                # Save alignment that matches with best hit (row[5]) in
                # top_deeparg_hit
                query_dict[row[3]].add_deeparg_hit(row[5])

                # Count amr switch
                diamond_amr = query_dict[
                    row[3]
                ].get_top_diamond_classification()
                deeparg_amr = query_dict[
                    row[3]
                ].get_top_deeparg_classification()
                if diamond_amr != deeparg_amr:
                    key = f"{diamond_amr}\t{deeparg_amr}"
            deeparg_list.clear()

            # Get info on query alignments amr and all groups proportion
            all_labels_df = pd.DataFrame(
                columns=[
                    "clstr",
                    "arg",
                    "dom",
                    "superf",
                    "amr",
                    "count",
                    "amr ref count",
                    "clstr|amr ref count",
                    "dom|amr ref count",
                    "superf|amr ref count",
                    "Is Diamond Best Hit Label",
                    "Is Diamond Best-Hit Class",
                    "Is Most Frequent Class",
                    "Is DeepARG Class",
                    "Is Class with Biggest Superfam",
                    "Is Class with Biggest Cluster",
                    "Diamond Class",
                    "Diamond clstr",
                    "Diamond dom",
                    "Diamond superf",
                    "Most Frequent Class",
                    "DeepARG Class",
                    "Class with Biggest Superfam",
                    "Class with Biggest Cluster",
                    "Query",
                ]
            )
            query_count = 0
            for query in query_dict.values():
                if not query.is_deeparg_hit():
                    continue
                if not query.has_multiple_classes():
                    continue
                query_count += 1
                query_vector = query.create_query_vector(
                    reference_clstr_amr_count,
                    reference_dom_amr_count,
                    reference_super_amr_count,
                    reference_amr_count,
                )
                if not query_vector.has_multiple_possible_classes():
                    continue
                all_labels_df = pd.concat(
                    [all_labels_df, query_vector.get_label_counts()],
                    ignore_index=True,
                )

            all_labels_df.insert(loc=0, column="Model", value=model)
            all_labels_df.insert(
                loc=0, column="Alignment Identity", value=identity
            )
            all_labels_df.insert(loc=0, column="Sample ID", value=sample_id)

            if first:
                all_labels_df.to_csv(
                    path_or_buf="label_counts.tsv",
                    sep="\t",
                    float_format="{:.4f}".format,
                    index=False,
                )
            else:
                all_labels_df.to_csv(
                    path_or_buf="label_counts.tsv",
                    sep="\t",
                    float_format="{:.4f}".format,
                    mode="a",
                    header=False,
                    index=False,
                )
            first = False
