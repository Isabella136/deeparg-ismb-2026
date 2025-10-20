from differential_distribution_classes.domain import DomainContainer, Domain
from differential_distribution_classes.reference import Reference

class Alignment(DomainContainer):
    def __init__(self, row: list[str], ref_dict: dict[str, Reference]):
        self.matching_reference = ref_dict[row[1]]
        self.bitscore = float(row[11])
        self.start = int(row[8])
        self.end = int(row[9])
        self.al_len = int(row[3])

    def get_bitscore(self) -> float:
        return self.bitscore
    
    def get_coverage(self) -> float:
        return (float(self.al_len) / float(self.matching_reference.get_length()))
    
    def get_domains(self) -> list[Domain]:
        aligned_domains = list()
        for domain in self.matching_reference.get_domains():
            if (domain.start < self.end) and (domain.end > self.start):
                aligned_domains.append(domain)            
        return aligned_domains
    
    def get_name(self) -> str:
        return self.matching_reference.get_name()
    
    def get_cluster(self) -> int:
        return self.matching_reference.get_cluster()