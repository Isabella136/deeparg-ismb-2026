from differential_distribution_classes.query import QueryVector

class QuerySubGroups:
    query_vectors: list[QueryVector]

    def __init__(self):
        pass

class QueryGroup:
    query_vectors: list[QueryVector]
    query_subgroups: list[QuerySubGroups]

    def __init__(self):
        pass

class FeatureGroup:
    def __init__(self):
        pass

class FeatureGroupBuilder:
    feature_groups: list[FeatureGroup]
    query_groups: list[QueryGroup]

    def __init__(self, query_vectors: list[QueryVector], features: list[str]):
        pass