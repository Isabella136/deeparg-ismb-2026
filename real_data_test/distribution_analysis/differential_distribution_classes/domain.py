from abc import ABC, abstractmethod

class Domain:
    def __init__(self, row: list[str]) :
        self.cdd_acc = row[7]
        self.super_acc = row[10] if row[1] == "specific" else row[7]
        self.start = int(row[3])
        self.end = int(row[4])

class DomainContainer(ABC):
    @abstractmethod
    def get_domains(self) -> list[Domain] :
        pass
    
    @abstractmethod
    def get_name(self) -> str :
        pass
    
    @abstractmethod
    def get_cluster(self) -> int :
        pass
    
    def get_arg_id(self) -> str :
        return self.get_name().split("|")[-1]
    
    def get_classification(self) -> str :
        return self.get_name().split("|")[-2]
        
    
    def get_domain_identifiers(self) -> tuple[str, str, str, str] :
        cluster_key = "|".join((str(self.get_cluster()), self.get_classification()))
        self_cdd_accs = [domain.cdd_acc for domain in self.get_domains()]
        self_super_accs = [domain.super_acc for domain in self.get_domains()]
        if len(self_cdd_accs) == 0:
            self_cdd_accs.append(self.get_arg_id())
            self_super_accs.append(self.get_arg_id())
        gene_key = "|".join(("$".join(self_cdd_accs), self.get_arg_id(), self.get_classification()))
        cdd_key = "|".join(("$".join(self_cdd_accs), self.get_classification()))
        super_key = "|".join(("$".join(self_super_accs), self.get_classification()))

        return (cluster_key, gene_key, cdd_key, super_key)
        