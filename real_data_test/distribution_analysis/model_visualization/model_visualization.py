import pandas as pd
import numpy as np
import pickle
import keras
from keras import layers
from keras.optimizers import SGD
from keras.losses import CategoricalCrossentropy
import sys
import os

sys.path.append('../')

from differential_distribution_classes.reference import Reference
from differential_distribution_classes.query import Query, QueryVector
from helper_classes.feature_grouping import FeatureGroupBuilder

CDD_DIR = "../../../CDD_features/"
REF_LOC = "../../../data/database/v2/features.fasta"
CLSTR_LOC = "../../database_clustering/db40.clstr"
DEEPARG_HIT_FILE = "X.mapping.ARG"
ALIGNMENT_FILE = "X.align.daa.tsv"
SAMPLE_ID_FILE = "../../real_samples.txt"
PATH_TO_SS = "deeparg_results"
PATH_TO_LS = "spades/deeparg_results"

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
        cdd_info_list.clear()

    # Add additional references that don't have domains + save length information
    with open(REF_LOC, "r") as ref_buf:
        ref_list = ref_buf.read().split('\n')
    for ref_index in range((len(ref_list) - 1) // 2):
        name = ref_list[ref_index*2][1:]
        if name not in ref_dict.keys():
            ref_dict[name] = Reference(name)
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
    
    # Go through model one at a time
    for model in ["LS", "SS"]:
        query_vectors: list[QueryVector] = list()

        for sample_id in sample_id_list:
            # Get Diamond alignment information
            query_dict: dict[str, Query] = dict()
            with open("/".join((
                    "../..", "samples", sample_id, PATH_TO_LS if model=="LS" else PATH_TO_SS,
                    f"arg_alignment_identity_30", ALIGNMENT_FILE)), "r") as alignment_buf:
                alignment_list = alignment_buf.read().split('\n')
            for alignment in alignment_list[:-1]:
                row = alignment.split('\t')
                if (len(query_dict) > 0) and (row[0] in query_dict.keys()):
                    query_dict[row[0]].add_alignment(row, ref_dict, model == "LS")
                else:
                    query_dict[row[0]] = Query(row, ref_dict, model == "LS")
            alignment_list.clear()

            # Get DeepARG hit information
            with open("/".join((
                    "../..", "samples", sample_id, PATH_TO_LS if model=="LS" else PATH_TO_SS,
                    f"arg_alignment_identity_30", DEEPARG_HIT_FILE)), "r") as deeparg_buf:
                deeparg_list = deeparg_buf.read().split('\n')
            for hit in deeparg_list[1:-1]:
                row = hit.split('\t')
                
                # Mark queries that match with read_id (row[3]) as deeparg_hit
                # Save alignment that matches with best hit (row[5]) in top_deeparg_hit
                query_dict[row[3]].add_deeparg_hit(row[5])
            deeparg_list.clear()
            query_vectors.extend([query.create_query_vector() for query in query_dict.values() if query.passed_cov_threshold()])

        feature_group_builder = FeatureGroupBuilder(query_vectors)
        feature_group_builder.print_feature_matrix(model)
            

def model_visualizer(version: str):
    # Open metadata on input and output
    metadata_IO_DL = pickle.load(open("../../../data/model/v2/metadata_"+version+".pkl", "rb"))
    input_nodes = metadata_IO_DL['input_nodes']
    output_nodes = metadata_IO_DL['output_nodes']
    model_DL = pickle.load(open("../../../data/model/v2/model_"+version+".pkl", "rb"), encoding="latin1")
    input_features = metadata_IO_DL['features']
    output_class = metadata_IO_DL['Y_rev']

    # # Recreate model architecture into Keras 3.0
    # input_layer = layers.Input(shape=(input_nodes,))
    # hidden_1 = layers.Dense(2000)(input_layer)
    # drop_1 = layers.Dropout(0.5)(hidden_1)
    # hidden_2 = layers.Dense(1000)(drop_1)
    # drop_2 = layers.Dropout(0.5)(hidden_2)
    # hidden_3 = layers.Dense(500)(drop_2)
    # drop_3 = layers.Dropout(0.5)(hidden_3)
    # hidden_4 = layers.Dense(100)(drop_3)
    # output_layer = layers.Dense(output_nodes, activation="softmax")(hidden_4)
    # model = keras.Model(inputs=input_layer, outputs=output_layer)
    # model.compile(
    #     optimizer=SGD(learning_rate=0.01, momentum=0.9, nesterov=True),
    #     loss=CategoricalCrossentropy()
    # )

    
    # model.layers[1].set_weights(model_DL['dense1'])
    # model.layers[3].set_weights(model_DL['dense3'])
    # model.layers[5].set_weights(model_DL['dense5'])
    # model.layers[7].set_weights(model_DL['dense7'])
    # model.layers[8].set_weights(model_DL['dense8'])

    # model.save("model_"+version+".h5")
    print(input_features)
    print(output_class)

main()
