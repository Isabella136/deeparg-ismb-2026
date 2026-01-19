from classes.domain import DomainContainer, Domain

class Reference(DomainContainer):
    def __init__(self, arg):
        self.name = arg
        self.domains = list()

    # Add length info when parsing feature fasta file
    def add_length_info(self, length: int):
        self.length = length

    # Add cluster info when parsing 40% database clustering output
    def define_cluster(self, cluster: int):
        self.cluster = cluster

    # Add domain if reference has at least one hit from batch CD search
    def add_domain(self, row):
        self.domains.append(Domain(row))

    def get_length(self) -> int :
        return self.length

    def get_domains(self) -> list[Domain]:
        return self.domains
    
    def get_name(self) -> str:
        return self.name
    
    def get_cluster(self) -> int:
        return self.cluster