from differential_distribution_classes.alignment import Alignment
from differential_distribution_classes.reference import Reference
import pandas as pd
import numpy as np

class Query:
    deeparg_hit: bool
    top_deeparg_hit: int
    above_cov_thresh: bool
    top_diamond_alignment: int
    alignments: list[Alignment]
    name: str

    def __init__(self, row: list[str], ref_dict: dict[str, Reference], ls_model: bool):
        self.deeparg_hit = False
        self.name = row[0]
        new_alignment = Alignment(row, ref_dict)
        if (not ls_model) or (new_alignment.get_coverage() >= 0.8):
            self.above_cov_thresh = True
            self.alignments = [new_alignment]
            self.top_diamond_alignment = 0
        else:
            self.above_cov_thresh = False
            self.alignments = list()
        

    def add_alignment(self, row: list[str], ref_dict: dict[str, Reference], ls_model: bool):
        new_alignment = Alignment(row, ref_dict)
        if (not ls_model) or (new_alignment.get_coverage() >= 0.8):
            if not self.above_cov_thresh: 
                self.above_cov_thresh = True
                self.alignments.append(new_alignment)
                self.top_diamond_alignment = 0
            else:
                self.alignments.append(new_alignment)
                curr_top_alignment = self.get_top_diamond_alignment()
                if self.alignments[-1].get_bitscore() > curr_top_alignment.get_bitscore():
                    self.top_diamond_alignment = len(self.alignments) - 1

    def add_deeparg_hit(self, best_hit: str):
        self.deeparg_hit = True
        for index, alignment in enumerate(self.alignments):
            if alignment.get_name() == best_hit:
                self.top_deeparg_hit = index

                # If diamond != deeparg but bitscores are equal, set diamond equal to deeparg
                if not self.are_diamond_and_deeparg_the_same():
                    diamond_alignment = self.get_top_diamond_alignment()
                    if diamond_alignment.get_bitscore() == alignment.get_bitscore():
                        self.top_diamond_alignment = self.top_deeparg_hit
                        break
                else:
                    break
    
    def is_deeparg_hit(self) -> bool :
        return self.deeparg_hit
    
    def are_diamond_and_deeparg_the_same(self) -> bool :
        if not self.deeparg_hit:
            raise RuntimeError("You didn't check to see if query is a deeparg hit")
        
        return (self.top_diamond_alignment == self.top_deeparg_hit)
    
    def get_query_name(self) -> str :
        return self.name
    
    def get_all_alignments_names(self) -> list[str]:
        return [alignment.get_name() for alignment in self.alignments]
    
    def get_all_alignment_groupings(self) -> list[tuple[str]]:
        return [alignment.get_groupings() for alignment in self.alignments]
    
    def get_all_alignment_bitscores(self) -> list[float] :
        return [alignment.get_bitscore() for alignment in self.alignments]

    def get_top_diamond_alignment(self) -> Alignment :
        return self.alignments[self.top_diamond_alignment]
    
    def get_top_diamond_classification(self) -> str :
        return self.get_top_diamond_alignment().get_classification()

    def get_top_diamond_alignment_domain_identifiers(self) -> tuple[str, str, str, str] :
        """Returns clstr|amr, dom|arg|amr, dom|amr, and super|amr ids of best Diamond alignment, respectively"""
        return self.get_top_diamond_alignment().get_domain_identifiers()
    
    def get_top_deeparg_hit(self) -> Alignment :
        if not self.deeparg_hit:
            raise RuntimeError("You didn't check to see if query is a deeparg hit")
        
        return self.alignments[self.top_deeparg_hit]
    
    def get_top_deeparg_classification(self, error: bool = True) -> str : 
        if not error and not self.deeparg_hit:
            return "none"
        return self.get_top_deeparg_hit().get_classification()
    
    def get_top_deeparg_hit_domain_identifiers(self) -> tuple[str, str, str, str] :
        """Returns clstr|amr, dom|arg|amr, dom|amr, and super|amr ids of DeepARG hit, respectively"""
        return self.get_top_deeparg_hit().get_domain_identifiers()
    
    def passed_cov_threshold(self) -> bool :
        return len(self.alignments) > 0
    
    def create_query_vector(self) -> "QueryVector":
        return QueryVector(self)
    
class QueryVector:
    feature_matrix: pd.DataFrame
    deeparg_class: str
    name: str
    
    def __init__(self, query: Query):
        self.deeparg_class = query.get_top_deeparg_classification(False)
        names = np.array(query.get_all_alignments_names())
        groupings = np.array(query.get_all_alignment_groupings())
        index = pd.MultiIndex.from_arrays(arrays=[
            groupings[:,0], groupings[:,1], groupings[:,2],
            groupings[:,3], groupings[:,4], names], names=[
                "clstr", "arg", "dom", "super", "amr", "names"])
        self.feature_matrix = pd.DataFrame(
            data=query.get_all_alignment_bitscores(), index=index)
        self.name = query.get_query_name()
        
    def has_multiple_possible_classes(self) -> bool :
        return len(self.feature_matrix.index.get_level_values("amr").drop_duplicates()) > 1
    
    def is_deeparg_hit(self) -> bool :
        return self.deeparg_class != "none"
    
    def get_deeparg_class(self) -> str :
        return self.deeparg_class
    
    def get_refs_and_matrix_indices(self) -> dict[str, pd.Series] :
        ref_matrix_indices = dict()
        matrix_indices = self.feature_matrix.index.to_frame()
        for row in matrix_indices.iterrows():
            ref_matrix_indices.update({row[1]["names"]:row[1]})
        return ref_matrix_indices
    
    def get_query_name(self) -> str :
        return self.name

    def get_feature_matrix(self) -> pd.DataFrame :
        return self.feature_matrix
        
    def get_reference_names(self) -> set[str] :
        return set(self.feature_matrix.index.get_level_values("names").values)