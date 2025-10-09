from differential_distribution_classes.alignment import Alignment
from differential_distribution_classes.reference import Reference

class Query:
    def __init__(self, row: list[str], ref_dict: dict[str, Reference]):
        self.alignments = [Alignment(row, ref_dict)]
        self.top_diamond_alignment = 0
        self.deeparg_hit = False

    def add_alignment(self, row: list[str], ref_dict: dict[str, Reference]):
        self.alignments.append(Alignment(row, ref_dict))
        if self.alignments[-1].get_bitscore() > self.get_top_diamond_alignment().get_bitscore():
            self.top_diamond_alignment = len(self.alignments) - 1

    def add_deeparg_hit(self, best_hit: str):
        self.deeparg_hit = True
        for index, alignment in enumerate(self.alignments):
            if alignment.get_name() == best_hit:
                self.top_deeparg_hit = index
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

    def get_top_diamond_alignment_domain_identifiers(self) -> list[tuple[str, str, str, str]] :
        return self.get_top_diamond_alignment().get_domain_identifiers()
    
    def get_top_deeparg_hit(self) -> Alignment :
        if not self.deeparg_hit:
            raise RuntimeError("You didn't check to see if query is a deeparg hit")
        
        return self.alignments[self.top_deeparg_hit]
    
    def get_top_deeparg_classification(self) -> str : 
        return self.get_top_deeparg_hit().get_classification()
    
    def get_top_deeparg_hit_domaind_identifiers(self) -> list[tuple[str, str, str, str]] :
        return self.get_top_deeparg_hit().get_domain_identifiers()