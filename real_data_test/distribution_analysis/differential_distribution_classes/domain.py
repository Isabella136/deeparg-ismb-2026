from abc import ABC, abstractmethod

class Domain:
    def __init__(self, row: list[str]) :
        self.dom_acc = row[7]
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
    
    def get_groupings(self) -> tuple[str, str, str, str, str] :
        """Returns clstr, arg, dom, super, and amr, respectively"""
        clstr_id = str(self.get_cluster())
        arg_id = self.get_arg_id()
        dom_id = self.get_domain_ids()
        super_id = self.get_super_ids()
        amr_class_id = self.get_classification()
        return (clstr_id, arg_id, dom_id, super_id, amr_class_id)

    def get_domain_ids(self) -> str:
        self_dom_accs = [domain.dom_acc for domain in self.get_domains()]
        if len(self_dom_accs) == 0:
            self_dom_accs.append(self.get_arg_id())
        return "$".join(self_dom_accs)
    
    def get_super_ids(self) -> str:
        self_super_accs = [domain.super_acc for domain in self.get_domains()]
        if len(self_super_accs) == 0:
            self_super_accs.append(self.get_arg_id())
        return "$".join(self_super_accs)
    
    def get_domain_identifiers(self) -> tuple[str, str, str, str] :
        """Returns clstr|amr, dom|arg|amr, dom|amr, and super|amr ids, respectively"""
        clstr_id = "|".join((str(self.get_cluster()), self.get_classification()))
        arg_id = "|".join((self.get_domain_ids(), self.get_arg_id(), self.get_classification()))
        dom_id = "|".join((self.get_domain_ids(), self.get_classification()))
        super_id = "|".join((self.get_super_ids(), self.get_classification()))

        return (clstr_id, arg_id, dom_id, super_id)
        