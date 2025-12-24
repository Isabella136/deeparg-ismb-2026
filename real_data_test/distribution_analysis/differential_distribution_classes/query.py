from differential_distribution_classes.alignment import Alignment
from differential_distribution_classes.reference import Reference
from multiprocessing import shared_memory
import pandas as pd
import numpy as np

class Query:
    deeparg_hit: bool
    top_deeparg_hit: int
    above_cov_thresh: bool
    top_diamond_alignment: int
    alignments: list[Alignment]
    name: str

    def __init__(
            self, row: list[str], ref_dict: dict[str, Reference], ls_model: bool, sample: str):
        self.deeparg_hit = False
        self.name = f"{sample}|{row[0]}"
        # When using SS model, not all 12279 features are used; should skip alignments to those
        if row[1] not in ref_dict:
            self.above_cov_thresh = False
            self.alignments = list()
        else: 
            new_alignment = Alignment(row, ref_dict)
            if (not ls_model) or (new_alignment.get_coverage() >= 0.8):
                self.above_cov_thresh = True
                self.alignments = [new_alignment]
                self.top_diamond_alignment = 0
            # When using LS model, alignments with less than 80% coverage are thrown out
            else:
                self.above_cov_thresh = False
                self.alignments = list()
        
    # If applicable, add alignment to query object
    def add_alignment(self, row: list[str], ref_dict: dict[str, Reference], ls_model: bool):
        if row[1] in ref_dict:
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

    # Point out the DeepARG hit
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
    
    def get_all_alignments(self) -> list[Alignment] :
        return self.alignments
    
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
    
    def create_query_decision_vector(
            self, shared_arrays_data: list[tuple]) -> "QueryDecisionVector":
        return QueryDecisionVector(self, shared_arrays_data)
    
class QueryDecisionVector:
    deeparg_class: str
    decision_vector: pd.DataFrame

    def __init__(self, query: Query, shared_arrays_data: list[tuple]):
        # Retrieve shared memory
        features_memory = shared_memory.SharedMemory(name='features')
        features = np.ndarray(
            shape=shared_arrays_data[0][0], dtype=shared_arrays_data[0][1],
            buffer=features_memory.buf)
        
        clstr_array_memory = shared_memory.SharedMemory(name='clstr')
        clstr_array = np.ndarray(
            shape=shared_arrays_data[1][0], dtype=shared_arrays_data[1][1],
            buffer=clstr_array_memory.buf)
        
        dom_array_memory = shared_memory.SharedMemory(name='dom')
        dom_array = np.ndarray(
            shape=shared_arrays_data[2][0], dtype=shared_arrays_data[2][1],
            buffer=dom_array_memory.buf)
        
        super_array_memory = shared_memory.SharedMemory(name='super')
        super_array = np.ndarray(
            shape=shared_arrays_data[3][0], dtype=shared_arrays_data[3][1],
            buffer=super_array_memory.buf)
        
        amr_idx_array_memory = shared_memory.SharedMemory(name='amr')
        amr_idx_array = np.ndarray(
            shape=shared_arrays_data[4][0], dtype=shared_arrays_data[4][1],
            buffer=amr_idx_array_memory.buf)

        self.deeparg_class = query.get_top_deeparg_classification(False)
        alignments = np.array(query.get_all_alignments())
        get_ref_name = np.vectorize(lambda a: a.get_name())
        get_clstr_amr = np.vectorize(lambda a: f"{a.get_cluster()}|{a.get_classification()}")
        get_dom_amr = np.vectorize(lambda a: f"dom:{a.get_domain_ids()}|{a.get_classification()}")
        get_super_amr = np.vectorize(lambda a: f"super:{a.get_super_ids()}|{a.get_classification()}")
        get_amr = np.vectorize(lambda a: a.get_classification())
        get_bitscore = np.vectorize(lambda a: a.get_bitscore())
        alignments_matrix = pd.DataFrame(data={
            'ref'   : get_ref_name(alignments),
            'clstr' : get_clstr_amr(alignments),
            'dom'   : get_dom_amr(alignments),
            'super' : get_super_amr(alignments),
            'amr'   : get_amr(alignments),
            'bit'   : get_bitscore(alignments)})
        cols = np.concat((features, "final class"), axis=None)
        self.decision_vector = pd.DataFrame(
            data=np.full(shape=(1,len(cols)), fill_value=0.0), columns=cols)
        self.decision_vector.at[0,"final class"] = float(np.extract(
            amr_idx_array[:,0]==self.deeparg_class, amr_idx_array[:,1])[0])
        for row in alignments_matrix.iterrows():
            self.decision_vector.at[0,row[1]['ref']] = row[1]['bit']
            self.decision_vector.at[0,row[1]['clstr']] += (
                row[1]['bit'] / float(np.extract(clstr_array[:,0]==row[1]['clstr'], clstr_array[:,1])[0]))
            self.decision_vector.at[0,row[1]['dom']] += (
                row[1]['bit'] / float(np.extract(dom_array[:,0] == row[1]['dom'][4:], dom_array[:,1])[0]))
            self.decision_vector.at[0,row[1]['super']] += (
                row[1]['bit'] / float(np.extract(super_array[:,0] == row[1]['super'][6:], super_array[:,1])[0]))
        
        # Close shared memory to avoid leaks
        amr_idx_array_memory.close()
        clstr_array_memory.close()
        super_array_memory.close()
        dom_array_memory.close()
        features_memory.close()
            
    def get_decision_vector(self) -> pd.DataFrame :
        return self.decision_vector
    
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