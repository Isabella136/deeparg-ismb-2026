from differential_distribution_classes.domain import DomainContainer, Domain

class Reference(DomainContainer):
    def __init__(self, arg):
        if isinstance(arg, str):
            self.name = arg
            self.domains = list()
        elif isinstance(arg, list):
            self.name = arg[0].split('>')[-1]
            self.domains = [Domain(arg)]
        else:
            raise RuntimeError("Reference initializer was given the wron argument")

    def add_length_info(self, length: int):
        self.length = length

    def define_cluster(self, cluster: int):
        self.cluster = cluster

    def add_domain(self, row):
        self.domains.append(Domain(row))

    def get_domains(self) -> list[Domain]:
        return self.domains
    
    def get_name(self) -> str:
        return self.name
    
    def get_cluster(self) -> int:
        return self.cluster