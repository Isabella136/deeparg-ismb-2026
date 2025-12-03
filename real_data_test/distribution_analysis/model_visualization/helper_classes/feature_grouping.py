from differential_distribution_classes.query import QueryVector
from differential_distribution_classes.reference import Reference

class QuerySubGroups:
    query_vectors: list[QueryVector]
    subgroup_amr = str

    def __init__(self, amr: str):
        self.query_vectors = list()
        self.subgroup_amr = amr

    def add_vector(self, query_vector: QueryVector):
        if query_vector.get_deeparg_class() != self.subgroup_amr:
            raise RuntimeError("Query's DeepARG AMR class should be same as current subgroup")
        self.query_vectors.append(query_vector)

class QueryGroup:
    query_vectors: list[QueryVector]
    query_subgroups: dict[str, QuerySubGroups]

    def __init__(self, query_vector: QueryVector):
        self.query_vectors = [query_vector]

    def create_subgroups(self):
        for vec in self.query_vectors:
            amr = vec.get_deeparg_class()
            curr_subgroup = self.query_subgroups.get(amr, QuerySubGroups(amr))
            curr_subgroup.add_vector(vec)

class FeatureGroup:
    def __init__(self):
        pass

class FeatureGroupBuilder:
    feature_groups: list[FeatureGroup]
    query_groups: list[QueryGroup]

    def __init__(self, query_vectors: list[QueryVector]):
        multi_amr_queries = list[QueryGroup]
        one_amr_queries = list[QueryGroup]
        for vec in query_vectors:
            curr_refs = vec.get_reference_names()
            if not vec.has_multiple_possible_classes():
                
                continue
            else:
                continue