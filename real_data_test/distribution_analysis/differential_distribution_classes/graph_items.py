from abc import ABC, abstractmethod

class Vertex(ABC):
    @abstractmethod
    def add_edge(self, state: bool, edge: 'Edge'):
        pass
  
class Edge:
    pair_count = 0
    def __init__(self, vertex_a: Vertex, vertex_b: Vertex):
        self.pair_count += 1
        self.vertex_a = vertex_a
        self.vertex_b = vertex_b

class State:
    edges = list()
    state_count = 0

    def __init__(self, edge: Edge):
        self.edges.append(edge)
        self.state_count += 1

class TrioVertex(Vertex):
    cdd_acc = None
    arg_id = None
    amr_cl = None

    state_a = None
    state_b = None

    ref_count = 0

class ClstrVertex(Vertex):
    clstr = 0
    amr_cl = None

    state_a = None
    state_b = None

    ref_count = 0

class DomainVertex(Vertex):
    cdd_acc = None
    amr_cl = None

    state_a = None
    state_b = None

    ref_count = 0

class SuperVertex(Vertex):
    super_acc = None
    amr_cl = None

    state_a = None
    state_b = None

    ref_count = 0