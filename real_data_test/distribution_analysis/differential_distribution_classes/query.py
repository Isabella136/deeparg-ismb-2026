from differential_distribution_classes.alignment import Alignment
from differential_distribution_classes.reference import Reference

class Query:
    deeparg_hit: bool
    top_deeparg_hit: int
    above_cov_thresh: bool
    top_diamond_alignment: int
    alignments: list[Alignment]

    def __init__(self, row: list[str], ref_dict: dict[str, Reference], ls_model: bool):
        self.deeparg_hit = False
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
    
    def get_top_deeparg_classification(self) -> str : 
        return self.get_top_deeparg_hit().get_classification()
    
    def get_top_deeparg_hit_domain_identifiers(self) -> tuple[str, str, str, str] :
        """Returns clstr|amr, dom|arg|amr, dom|amr, and super|amr ids of DeepARG hit, respectively"""
        return self.get_top_deeparg_hit().get_domain_identifiers()