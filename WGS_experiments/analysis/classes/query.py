from classes.alignment import Alignment
from classes.reference import Reference
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
    
    def has_multiple_classes(self) -> bool:
        for a in self.alignments:
            if a.get_classification() != self.get_top_diamond_classification():
                return True
        return False
    
    def get_all_alignments_names(self) -> list[str]:
        return [alignment.get_name() for alignment in self.alignments]
    
    def get_all_alignment_annotations(self) -> list[tuple[str]]:
        return [alignment.get_annotations() for alignment in self.alignments]
    
    def get_all_alignment_bitscores(self) -> list[float] :
        return [alignment.get_bitscore() for alignment in self.alignments]

    def get_top_diamond_alignment(self) -> Alignment :
        return self.alignments[self.top_diamond_alignment]
    
    def get_top_diamond_classification(self) -> str :
        return self.get_top_diamond_alignment().get_classification()

    def get_top_diamond_alignment_annotation(self) -> tuple[str, str, str, str] :
        """Returns clstr, arg, dom, super and amr of best Diamond alignment, respectively"""
        return self.get_top_diamond_alignment().get_annotations()
    
    def get_top_deeparg_hit(self) -> Alignment :
        if not self.deeparg_hit:
            raise RuntimeError("You didn't check to see if query is a deeparg hit")
        
        return self.alignments[self.top_deeparg_hit]
    
    def get_top_deeparg_classification(self, error: bool = True) -> str : 
        if not error and not self.deeparg_hit:
            return "none"
        return self.get_top_deeparg_hit().get_classification()
    
    def passed_cov_threshold(self) -> bool :
        return len(self.alignments) > 0
    
    def create_query_vector(
            self,
            reference_clstr_amr_count: dict[str, int],
            reference_dom_amr_count: dict[str, int],
            reference_super_amr_count: dict[str, int],
            reference_amr_count: dict[str, int]) -> "QueryVector":
        return QueryVector(
            self,
            reference_clstr_amr_count,
            reference_dom_amr_count,
            reference_super_amr_count,
            reference_amr_count)
    
class QueryVector:
    label_count: pd.DataFrame

    deeparg_class: str
    diamond_class: str
    most_freq_class: str

    deeparg_labels: tuple[str,str,str,str,str]
    diamond_labels: tuple[str,str,str,str,str]

    name: str
    
    def __init__(
            self, query: Query, 
            reference_clstr_amr_count: dict[str, int],
            reference_dom_amr_count: dict[str, int],
            reference_super_amr_count: dict[str, int],
            reference_amr_count: dict[str, int]):
        names = np.array(query.get_all_alignments_names())
        annotations = np.array(query.get_all_alignment_annotations())
        bitscores = np.array(query.get_all_alignment_bitscores())

        index = pd.MultiIndex.from_arrays(arrays=[
            annotations[:,0], annotations[:,1], annotations[:,2],
            annotations[:,3], annotations[:,4], names], names=[
                "clstr", "arg", "dom", "super", "amr", "names"])
        
        feature_matrix = pd.DataFrame(
            data={"bitscore": bitscores}, index=index)

        self.label_count = (feature_matrix
            .groupby(level=["clstr", "arg", "dom", "super", "amr"])["bitscore"]
            .count()
            .reset_index()
            .rename(columns={"bitscore":"count"}))
        self.label_count["amr ref count"] = self.label_count.apply(
            lambda x: reference_amr_count[x['amr']], axis=1)
        self.label_count["clstr|amr ref count"] = self.label_count.apply(
            lambda x: reference_clstr_amr_count['|'.join((x['clstr'], x['amr']))], axis=1)
        self.label_count["dom|amr ref count"] = self.label_count.apply(
            lambda x: reference_dom_amr_count['|'.join((x['dom'], x['amr']))], axis=1)
        self.label_count["super|amr ref count"] = self.label_count.apply(
            lambda x: reference_super_amr_count['|'.join((x['super'], x['amr']))], axis=1)
        
        self.name = query.get_query_name()

        self.deeparg_class = query.get_top_deeparg_classification(False)
        self.diamond_class = query.get_top_diamond_classification()
        self.most_freq_class = (
            self.label_count.groupby(by=["amr"])["count"].sum().idxmax()
            if self.label_count.groupby(by=["amr"])["count"].sum().max() > self.label_count.loc[
                self.label_count["amr"]==self.deeparg_class]["count"].sum() else self.deeparg_class)
        self.class_with_big_super = (
            self.label_count[["super", "amr", "super|amr ref count"]].drop_duplicates().groupby(by=["amr"])["super|amr ref count"].max().idxmax()
            if self.label_count[["super", "amr", "super|amr ref count"]].drop_duplicates().groupby(by=["amr"])["super|amr ref count"].max().max() > self.label_count.loc[
                self.label_count["amr"]==self.deeparg_class]["super|amr ref count"].max() else self.deeparg_class)
        self.class_with_big_clstr = (
            self.label_count[["clstr", "amr", "clstr|amr ref count"]].drop_duplicates().groupby(by=["amr"])["clstr|amr ref count"].max().idxmax()
            if self.label_count[["clstr", "amr", "clstr|amr ref count"]].drop_duplicates().groupby(by=["amr"])["clstr|amr ref count"].max().max() > self.label_count.loc[
                self.label_count["amr"]==self.deeparg_class]["clstr|amr ref count"].max() else self.deeparg_class)

        self.diamond_labels = query.get_top_diamond_alignment_annotation()

        self.label_count["Is Diamond Best Hit Label"] = self.label_count.apply(
            lambda x: tuple(x["clstr":"amr"].to_list()) == self.diamond_labels, axis=1)
        self.label_count["Is Diamond Best-Hit Class"] = self.label_count.apply(
            lambda x: x["amr"] == self.diamond_class, axis=1)
        self.label_count["Is Most Frequent Class"] = self.label_count.apply(
            lambda x: x["amr"] == self.most_freq_class, axis=1)
        self.label_count["Is DeepARG Class"] = self.label_count.apply(
            lambda x: x["amr"] == self.deeparg_class, axis=1)
        self.label_count["Is Class with Biggest Superfam"] = self.label_count.apply(
            lambda x: x["amr"] == self.class_with_big_super, axis=1)
        self.label_count["Is Class with Biggest Cluster"] = self.label_count.apply(
            lambda x: x["amr"] == self.class_with_big_clstr, axis=1)
        

        self.label_count["Diamond Class"] = self.diamond_class
        self.label_count["Diamond clstr"] = self.diamond_labels[0]
        self.label_count["Diamond dom"] = self.diamond_labels[2]
        self.label_count["Diamond super"] = self.diamond_labels[3]
        self.label_count["Most Frequent Class"] = self.most_freq_class
        self.label_count["DeepARG Class"] = self.deeparg_class
        self.label_count["Class with Biggest Superfam"] = self.class_with_big_super
        self.label_count["Class with Biggest Cluster"] = self.class_with_big_clstr
        self.label_count["Query"] = self.name
        
    def has_multiple_possible_classes(self) -> bool :
        return len(self.label_count["amr"].drop_duplicates()) > 1
    
    def is_deeparg_unexpected(self) -> bool :
        return (self.deeparg_class != self.diamond_class) or (self.deeparg_class != self.most_freq_class)
    
    def is_deeparg_hit(self) -> bool :
        return self.deeparg_class != "none"
    
    def get_deeparg_class(self) -> str :
        return self.deeparg_class
    
    def get_query_name(self) -> str :
        return self.name
    
    def get_label_counts(self) -> pd.DataFrame :
        return self.label_count
    
    def get_feature_matrix(self) -> pd.DataFrame :
        return self.feature_matrix