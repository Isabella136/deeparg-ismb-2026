from abc import ABC, abstractmethod

class Vertex(ABC):
    amr_cl: str
    state_a: "StateA"
    state_b: "StateB"
    ref_count: int

    def get_edge_from_a(self, vertex_b: "Vertex") -> "Edge":
        edge = self.state_a.get_edge_to_b(vertex_b)
        edge.add_to_count()
        return edge

    def add_edge_to_b(self, vertex_a: "Vertex", edge: "Edge"):
        self.state_b.add_edge_from_a(vertex_a, edge)

    def increment_states_but_not_pair(self):
        self.state_a.increment_state_count()
        self.state_b.increment_state_count()

    def get_state_a_to_state_b_info(self) -> str:
        edge_info = self.state_a.get_all_edges_info()
        state_a_count = self.state_a.get_state_count()
        return f"State A count\t{state_a_count}\tSwitches to different state B\t{
            "\t".join(f"{name} ({count})" for (name, count) in edge_info)}"
    
    def get_state_b_from_state_a_info(self) -> str:
        edge_info = self.state_b.get_all_edges_info()
        state_b_count = self.state_b.get_state_count()
        return f"State B count\t{state_b_count}\tSwitches from different state A\t{
            "\t".join(f"{name} ({count})" for (name, count) in edge_info)}"
    
    def get_adjacent_vertices(self) -> list["Vertex"] :
        adjacent_vertices = self.state_a.get_adjacent_vertices()
        adjacent_vertices.extend(self.state_b.get_adjacent_vertices())
        return adjacent_vertices

    @abstractmethod
    def get_name(self) -> str :
        pass
        
  
class Edge:
    pair_count: int
    vertex_a: Vertex
    vertex_b: Vertex

    def __init__(self, vertex_a: Vertex, vertex_b: Vertex):
        self.vertex_a = vertex_a
        self.vertex_b = vertex_b
        self.pair_count = 0

    def add_to_count(self):
        self.pair_count += 1

    def get_count(self) -> int :
        return self.pair_count
    
    def get_vertex_a(self) -> Vertex :
        return self.vertex_a
    
    def get_vertex_b(self) -> Vertex :
        return self.vertex_b

class State(ABC):
    edges: dict[str, Edge]
    state_count: int
    vertex: Vertex

    def __init__(self, vertex: Vertex):
        self.edges = dict()
        self.state_count = 0
        self.vertex = vertex

    def get_state_count(self) -> int :
        return self.state_count
    
    def increment_state_count(self):
        self.state_count += 1

    def get_all_edges_info(self) -> list[tuple[str, int]]:
        edge_info: list[tuple[str, int]] = list()

        for (state_b, edge) in self.edges.items():
            edge_info.append((state_b, edge.get_count()))

        return edge_info
    
    @abstractmethod
    def get_adjacent_vertices(self) -> list[Vertex]:
        pass
        

class StateA(State):
    def __init__(self, vertex: Vertex):
        super().__init__(vertex)
            
    def get_edge_to_b(self, vertex_b: Vertex) -> Edge :
        """Only should be called from function "get_edge_from_a" in Vertex object"""
        self.state_count += 1
        state_b_name = vertex_b.get_name()
        if state_b_name not in self.edges.keys():
            self.edges[state_b_name] = Edge(self.vertex, vertex_b)
        return(self.edges[state_b_name])
    
    def get_adjacent_vertices(self) -> list[Vertex]:
        return [e.get_vertex_b() for e in self.edges.values()]

class StateB(State):
    def __init__(self, vertex: Vertex):
        super().__init__(vertex)

    def add_edge_from_a(self, vertex_a: Vertex, edge: Edge):
        """Only should be called from function "add_edge_to_b" in Vertex object"""
        self.state_count += 1
        state_a_name = vertex_a.get_name()
        if state_a_name not in self.edges.keys():
            self.edges[state_a_name] = edge

    def get_adjacent_vertices(self) -> list[Vertex]:
        return [e.get_vertex_a() for e in self.edges.values()]

class TrioVertex(Vertex):
    dom_acc: str
    arg_id: str

    def __init__(self, id: str, ref_count: int):
        (self.dom_acc, self.arg_id, self.amr_cl) = tuple(id.split("|"))
        self.ref_count = ref_count
        self.state_a = StateA(self)
        self.state_b = StateB(self)

    def get_name(self) -> str:
        return "|".join((self.dom_acc, self.arg_id, self.amr_cl))

class ClstrVertex(Vertex):
    clstr: str

    def __init__(self, id: str, ref_count: int):
        (self.clstr, self.amr_cl) = tuple(id.split("|"))
        self.ref_count = ref_count
        self.state_a = StateA(self)
        self.state_b = StateB(self)

    def get_name(self) -> str:
        return "|".join((self.clstr, self.amr_cl))

class DomainVertex(Vertex):
    dom_acc: str

    def __init__(self, id: str, ref_count: int):
        (self.dom_acc, self.amr_cl) = tuple(id.split("|"))
        self.ref_count = ref_count
        self.state_a = StateA(self)
        self.state_b = StateB(self)

    def get_name(self) -> str:
        return "|".join((self.dom_acc, self.amr_cl))

class SuperVertex(Vertex):
    super_acc: str

    def __init__(self, id: str, ref_count: int):
        (self.super_acc, self.amr_cl) = tuple(id.split("|"))
        self.ref_count = ref_count
        self.state_a = StateA(self)
        self.state_b = StateB(self)

    def get_name(self) -> str:
        return "|".join((self.super_acc, self.amr_cl))