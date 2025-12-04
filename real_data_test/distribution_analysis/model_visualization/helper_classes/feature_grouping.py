import sys

sys.path.append('../')

from differential_distribution_classes.query import QueryVector
from differential_distribution_classes.reference import Reference
import pandas as pd
import numpy as np

class QuerySubGroups:
    query_vectors: dict[str, QueryVector]
    subgroup_amr: str
    feature_matrix: pd.DataFrame

    def __init__(self, amr: str):
        self.query_vectors = dict()
        self.subgroup_amr = amr

    def add_vector(self, query_vector: tuple[str:QueryVector]):
        if query_vector[1].get_deeparg_class() != self.subgroup_amr:
            raise RuntimeError("Query's DeepARG AMR class should be same as current subgroup")
        self.query_vectors.update({query_vector[0]:query_vector[1]})

    def get_feature_matrix(self) -> pd.DataFrame:
        return self.feature_matrix
    
    def make_feature_matrix(self) -> pd.DataFrame: 
        ref_matrix_indices: dict[str, pd.Series] = dict()
        query_names: list[str] = list()
        for query_vector in self.query_vectors.items():
            ref_matrix_indices.update(query_vector[1].get_refs_and_matrix_indices())
            query_names.append(query_vector[0])
        new_indices_df: pd.DataFrame = pd.DataFrame(
            columns=["query", "clstr", "arg", "dom", "super", "amr", "names"])
        for name in query_names:
            batch = pd.DataFrame(
                ref_matrix_indices.values(), 
                columns=["clstr", "arg", "dom", "super", "amr", "names"])
            batch.insert(0, column="query", value=name)
            new_indices_df = pd.concat([new_indices_df, batch])
        new_indices_df.insert(0, column="str", value=self.subgroup_amr)
        new_indices: pd.MultiIndex = pd.MultiIndex.from_frame(new_indices_df)
        self.feature_matrix = pd.DataFrame(
            np.zeros(shape=[len(new_indices), 1]), index=new_indices, columns=["bit-score"])
        for query_vector in self.query_vectors.items():
            for row in query_vector[1].get_feature_matrix().iterrows():
                index_as_list = [self.subgroup_amr, query_vector[0]]
                index_as_list.extend(list(row[0]))
                self.feature_matrix.at[tuple(index_as_list),"bit-score"] = row[1][0]
        return self.get_feature_matrix()

class QueryGroup:
    query_vectors: dict[str, QueryVector]
    query_subgroups: dict[str, QuerySubGroups]

    def __init__(self, query_vectors: list[QueryVector]):
        self.query_vectors = {vec.get_query_name():vec for vec in query_vectors}

    def create_subgroups(self):
        self.query_subgroups = dict()
        for vec in self.query_vectors.items():
            amr = vec[1].get_deeparg_class()
            curr_subgroup = self.query_subgroups.get(amr, QuerySubGroups(amr))
            curr_subgroup.add_vector(vec)
            self.query_subgroups.update({amr:curr_subgroup})

    def print_feature_matrix(self, group_index: int, model: str):
        try:
            feature_matrix = list(self.query_subgroups.values())[0].get_feature_matrix()
        except:
            feature_matrix = list(self.query_subgroups.values())[0].make_feature_matrix()
        if len(self.query_subgroups.values()) > 1:
            for subgroup in list(self.query_subgroups.values())[1:]:
                try:
                    feature_matrix = pd.concat([feature_matrix, subgroup.get_feature_matrix()])
                except:
                    feature_matrix = pd.concat([feature_matrix, subgroup.make_feature_matrix()])
        feature_matrix.to_csv(f"{model}_query_group_{group_index}.csv", sep=',')

class FeatureGroup:
    def __init__(self):
        pass

class FeatureGroupBuilder:
    feature_groups: list[FeatureGroup]
    query_groups: list[QueryGroup]

    def __init__(self, query_vectors: list[QueryVector]):
        # key is tuple of indices to reference
        # value is list of queries that aligned to no references other than the ones in key
        multi_amr_queries : dict[set[str], list[QueryVector]] = dict()
        one_amr_queries : dict[set[str], list[QueryVector]] = dict()
        
        for vec in query_vectors:
            curr_ref_names = vec.get_reference_names()
            if vec.has_multiple_possible_classes():
                query_list_keys: list[set[str]] = [
                    key for key in multi_amr_queries.keys() if len(key.intersection(curr_ref_names)) > 0]
                big_vec_list = [vec]
                big_vec_list_key = curr_ref_names
                for key in query_list_keys:
                    big_vec_list_key.update(key)
                    big_vec_list.extend(multi_amr_queries.pop(key))
                multi_amr_queries[frozenset(big_vec_list_key)] = big_vec_list
            else:
                query_list_keys: list[set[str]] = [
                    key for key in one_amr_queries.keys() if len(key.intersection(curr_ref_names)) > 0]
                big_vec_list = [vec]
                big_vec_list_key = curr_ref_names
                for key in query_list_keys:
                    big_vec_list_key.update(key)
                    big_vec_list.extend(one_amr_queries.pop(key))
                one_amr_queries[frozenset(big_vec_list_key)] = big_vec_list

        self.query_groups = list()
        for query_vector_list in multi_amr_queries.values():
            self.query_groups.append(QueryGroup(query_vector_list))
            self.query_groups[-1].create_subgroups()
        for query_vector_list in one_amr_queries.values():
            self.query_groups.append(QueryGroup(query_vector_list))
            self.query_groups[-1].create_subgroups()

    def print_feature_matrix(self, model: str):
        for index, group in enumerate(self.query_groups):
            group.print_feature_matrix(index, model)