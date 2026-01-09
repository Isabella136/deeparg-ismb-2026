import pandas as pd
import numpy as np
import pickle
from sklearn import tree, feature_selection, model_selection
from matplotlib import pyplot as plt
from multiprocessing import Pool, shared_memory
import sys
import os

sys.path.append('../')

from differential_distribution_classes.reference import Reference
from differential_distribution_classes.query import Query, QueryDecisionVector #QueryVector
# from helper_classes.feature_grouping import FeatureGroupBuilder

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
    
if __name__ == '__main__':
    print(f"Pooling? We have {os.process_cpu_count()} cpus")
    with Pool(processes=os.process_cpu_count()-1) as p:
        print("Pooling!")
        for model in ['LS', 'SS']:
            # Open metadata to get ARG features and output
            ref_dict: dict[str, Reference] = dict()
            metadata_IO_DL = pickle.load(open("../../../data/model/v2/metadata_"+model+".pkl", "rb"))
            for ref in metadata_IO_DL["features"]:
                ref_dict[ref] = Reference(ref)
            amr_to_index: dict[str, int] = {amr: idx for idx, amr in metadata_IO_DL['Y_rev'].items()}
            amr_to_index.update({"none": len(metadata_IO_DL['Y_rev'])})

            # Create amr_to_index dict in shared memory
            amr_idx_array = np.array(list(amr_to_index.items()))
            try:
                amr_idx_memory = shared_memory.SharedMemory(name='amr')
                amr_idx_memory.unlink()
                amr_idx_memory = shared_memory.SharedMemory(
                    create=True, size=amr_idx_array.nbytes, name='amr')
            except:
                amr_idx_memory = shared_memory.SharedMemory(
                    create=True, size=amr_idx_array.nbytes, name='amr')
            shared_amr_idx_array = np.ndarray(
                amr_idx_array.shape, dtype=amr_idx_array.dtype, buffer=amr_idx_memory.buf)
            shared_amr_idx_array[:] = amr_idx_array[:]
            del amr_to_index, amr_idx_array
            
            # Read content of real_samples.txt to find biosample ID
            with open(SAMPLE_ID_FILE, "r") as sample_id_buf:
                sample_id_list = sample_id_buf.read().split('\n')

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
                    name = line.split('>')[-1].split('...')[0]
                    if name not in ref_dict.keys():
                        continue
                    ref_dict[name].define_cluster(cluster)

            clstr_amr_count: dict[str, int] = dict()
            dom_amr_count: dict[str, int] = dict()
            super_amr_count: dict[str, int] = dict()
            amr_count: dict[str, int] = dict()
            for ref in ref_dict.values():
                ref_ids = ref.get_domain_identifiers()
                (clstr_id, trio_id, dom_id, super_id) = ref_ids
                if clstr_id in clstr_amr_count.keys():
                    clstr_amr_count[clstr_id] += 1
                else:
                    clstr_amr_count[clstr_id] = 1

                trio_fields = trio_id.split('|')
                super_fields = super_id.split('|')
                dom_list = trio_fields[0].split('$')
                super_list = super_fields[0].split('$')
                if trio_fields[2] not in amr_count:
                    amr_count[trio_fields[2]] = 0
                amr_count[trio_fields[2]] += 1
                # We may have numerous domains in gene; should count each potential combination as well
                dom_combo_list = sorted_domains_combo(dom_list)
                super_combo_list = sorted_domains_combo(super_list)
                # We also need to count instances of arg as domain
                if trio_fields[0] != trio_fields[1]:
                    dom_combo_list.append(trio_fields[1])
                    super_combo_list.append(trio_fields[1])
                for i in range(len(dom_combo_list)):
                    dom_id = "|".join((dom_combo_list[i], trio_fields[2]))
                    if dom_id in dom_amr_count.keys():
                        dom_amr_count[dom_id] += 1
                    else:
                        dom_amr_count[dom_id] = 1
                    super_id = "|".join((super_combo_list[i], super_fields[1]))
                    if super_id in super_amr_count.keys():
                        super_amr_count[super_id] += 1
                    else:
                        super_amr_count[super_id] = 1

            # Create feature array in shared memory
            features = np.concat((
                list(ref_dict.keys()), list(clstr_amr_count.keys()), 
                ["dom:"+x for x in list(dom_amr_count.keys())], 
                ["super:"+x for x in list(super_amr_count.keys())],
                list(metadata_IO_DL['Y_rev'].values())), axis=None)
            try:
                features_memory = shared_memory.SharedMemory(name='features')
                features_memory.unlink()
                features_memory = shared_memory.SharedMemory(
                    create=True, size=features.nbytes, name='features')
            except:
                features_memory = shared_memory.SharedMemory(
                    create=True, size=features.nbytes, name='features')
            shared_features = np.ndarray(
                features.shape, dtype=features.dtype, buffer=features_memory.buf)
            shared_features[:] = features[:]

            # Create clstr amr count dict in shared memory
            clstr_array = np.array(list(clstr_amr_count.items()))
            try:
                clstr_memory = shared_memory.SharedMemory(name='clstr')
                clstr_memory.unlink()
                clstr_memory = shared_memory.SharedMemory(
                    create=True, size=clstr_array.nbytes, name='clstr')
            except:
                clstr_memory = shared_memory.SharedMemory(
                    create=True, size=clstr_array.nbytes, name='clstr')
            shared_clstr_array = np.ndarray(
                clstr_array.shape, dtype=clstr_array.dtype, buffer=clstr_memory.buf)
            shared_clstr_array[:] = clstr_array[:]
            del clstr_array

            # Create dom amr count dict in shared memory
            dom_array = np.array(list(dom_amr_count.items()))
            try:
                dom_memory = shared_memory.SharedMemory(name='dom')
                dom_memory.unlink()
                dom_memory = shared_memory.SharedMemory(
                    create=True, size=dom_array.nbytes, name='dom')
            except:
                dom_memory = shared_memory.SharedMemory(
                    create=True, size=dom_array.nbytes, name='dom')
            shared_dom_array = np.ndarray(
                dom_array.shape, dtype=dom_array.dtype, buffer=dom_memory.buf)
            shared_dom_array[:] = dom_array[:]
            del dom_array

            # Create super amr count dict in shared memory
            super_array = np.array(list(super_amr_count.items()))
            try:
                super_memory = shared_memory.SharedMemory(name='super')
                super_memory.unlink()
                super_memory = shared_memory.SharedMemory(
                    create=True, size=super_array.nbytes, name='super')
            except:
                super_memory = shared_memory.SharedMemory(
                    create=True, size=super_array.nbytes, name='super')
            shared_super_array = np.ndarray(
                super_array.shape, dtype=super_array.dtype, buffer=super_memory.buf)
            shared_super_array[:] = super_array[:]
            del super_array

            shared_arrays_data = [
                (shared_features.shape, shared_features.dtype),
                (shared_clstr_array.shape, shared_clstr_array.dtype),
                (shared_dom_array.shape, shared_dom_array.dtype),
                (shared_super_array.shape, shared_super_array.dtype),
                (shared_amr_idx_array.shape, shared_amr_idx_array.dtype)]
            
            query_decision_vectors: list[QueryDecisionVector] = list()
            for sample_id in sample_id_list:
                # Get Diamond alignment information
                print(f"Starting with sample {sample_id}")
                sys.stdout.flush()
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
                        query_dict[row[0]] = Query(row, ref_dict, model == "LS", sample_id)
                del alignment_list
                print("Alignment is done")
                sys.stdout.flush()
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
                del deeparg_list
                print("DeepARG is done")
                sys.stdout.flush()
                query_decision_vectors_temp = p.starmap(Query.create_query_decision_vector, [
                    (query, shared_arrays_data) for query in query_dict.values() if query.passed_cov_threshold()])
                print("Decision Vectors are made")
                sys.stdout.flush()
                del query_dict
                query_decision_vectors.extend(query_decision_vectors_temp)
                print("Decision Vectors are saved")
                sys.stdout.flush()

            # Close and unlink shared memory to avoid leaks
            features_memory.close()
            features_memory.unlink()
            amr_idx_memory.close()
            amr_idx_memory.unlink()
            clstr_memory.close()
            clstr_memory.unlink()
            super_memory.close()
            super_memory.unlink()
            dom_memory.close()
            dom_memory.unlink()

            print("Making matrix")
            sys.stdout.flush()
            query_decision_matrix = pd.concat([query.get_decision_vector() for query in query_decision_vectors])
            Y = query_decision_matrix["final class"]

            for curr_features, fig_name in [
                    (features, "all_features"), 
                    (list(clstr_amr_count.keys()), "cluster"), 
                    (["dom:"+x for x in list(dom_amr_count.keys())], "domain"), 
                    (["super:"+x for x in list(super_amr_count.keys())], "super")]:
                X = query_decision_matrix[curr_features]
                print("Remove useless features")
                sys.stdout.flush()
                selector = feature_selection.VarianceThreshold()
                print(f"No. of features before: {len(curr_features)}")
                sys.stdout.flush()
                X_new = selector.fit_transform(X, Y)
                remaining_features = selector.get_feature_names_out(curr_features)
                print(f"No. of features after variance selection: {len(remaining_features)}")
                sys.stdout.flush()
                Y_dedup = Y.drop_duplicates().sort_values()
                get_amr = np.vectorize(lambda i: dict(metadata_IO_DL['Y_rev']).get(int(i), "none"))
                amr_vals = get_amr(Y_dedup)
                # print("Set up train and test")
                # X_train, X_test, Y_train, Y_test = model_selection.train_test_split(
                #     X_new, Y, test_size=0.1)
                print("Making trees")
                sys.stdout.flush()
                decision_tree = tree.DecisionTreeClassifier(max_depth = 20 if model=="LS" else 10).fit(X_new, Y)
                print("Ploting tree")
                sys.stdout.flush()
                fig = plt.figure(figsize=(384,108), dpi=100)
                ax = fig.add_axes((0,0,1,1))
                tree.plot_tree(decision_tree, ax=ax, feature_names=remaining_features, class_names=amr_vals)
                print("Saving tree")
                sys.stdout.flush()
                fig.savefig(f"decision_tree_{fig_name}_{model}.png")
                sys.stdout.flush()

        print("DONE!")