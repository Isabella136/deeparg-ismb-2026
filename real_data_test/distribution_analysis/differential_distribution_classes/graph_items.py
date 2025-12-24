from abc import ABC, abstractmethod
import numpy as np

class Vertex(ABC):
    amr_cl: str
    state_a: "StateA"
    state_b: "StateB"
    ref_count: int

    # Make and get switch pair edge to new label when switching from vertex's label
    def get_edge_from_a(self, vertex_b: "Vertex") -> "Edge":
        edge = self.state_a.get_edge_to_b(vertex_b)
        edge.add_to_count()
        return edge

    # Add switch pair edge from old label when switching to vertex's label
    def add_edge_to_b(self, vertex_a: "Vertex", edge: "Edge"):
        self.state_b.add_edge_from_a(vertex_a, edge)

    # When no switch is happening, but you still want to keep track of counts
    def increment_states_but_not_pair(self):
        self.state_a.increment_state_count()
        self.state_b.increment_state_count()

    # Retrieve info on switch pairs from vertex to new labels
    def get_state_a_to_state_b_info(self) -> str:
        edge_info = self.state_a.get_all_edges_info()
        state_a_count = self.state_a.get_state_count()
        return f"State A count\t{state_a_count}\tSwitches to different state B\t{
            "\t".join(f"{name} ({count})" for (name, count) in edge_info)}"
    
    # Retrieve info on switch pairs to vertex from old labels
    def get_state_b_from_state_a_info(self) -> str:
        edge_info = self.state_b.get_all_edges_info()
        state_b_count = self.state_b.get_state_count()
        return f"State B count\t{state_b_count}\tSwitches from different state A\t{
            "\t".join(f"{name} ({count})" for (name, count) in edge_info)}"
    
    # Get edges representing switch pairs to or from vertex
    def get_edges(self) -> dict[(str, str), "Edge"]:
        edges = self.state_a.get_edges()
        edges.update(self.state_b.get_edges())
        return edges
    
    # Get vertices representing labels that vertex switches to or from
    def get_adjacent_vertices(self) -> dict[str, "Vertex"] :
        adjacent_vertices = self.state_a.get_adjacent_vertices()
        adjacent_vertices.update(self.state_b.get_adjacent_vertices())
        return adjacent_vertices
    
    # Get count when vertex is Diamond best hit
    def get_state_a_count(self) -> np.double :
        return np.double(self.state_a.get_state_count())
    
    # Get count when vertex is DeepARG hit
    def get_state_b_count(self) -> np.double :
        return np.double(self.state_b.get_state_count())
    
    # Get vertex's reference count
    def get_state_i_count(self) -> np.double :
        return np.double(self.ref_count)

    # Defined in children classes, get label/vertex name
    @abstractmethod
    def get_name(self) -> str :
        pass
        
    @abstractmethod
    def has_arg_label(self) -> bool:
        pass
  
class Edge:
    pair_count: int
    vertex_a: Vertex
    vertex_b: Vertex
    has_arg_label: bool

    # Create edge from vertex a to vertex b to represent pair switch
    def __init__(self, vertex_a: Vertex, vertex_b: Vertex):
        self.vertex_a = vertex_a
        self.vertex_b = vertex_b
        self.pair_count = 0
        self.has_arg_label = vertex_a.has_arg_label() or vertex_b.has_arg_label()

    # Increment each time there is a switch from a to b
    def add_to_count(self):
        self.pair_count += 1

    # Get pair switch count
    def get_count(self) -> int :
        return self.pair_count
    
    # Get vertex a, representing Diamond best-hit label
    def get_vertex_a(self) -> Vertex :
        return self.vertex_a
    
    # Get vertex b, representing DeepARG hit label
    def get_vertex_b(self) -> Vertex :
        return self.vertex_b

class State(ABC):
    edges: dict[str, Edge]
    state_count: int
    vertex: Vertex

    # Create state object for given vertex
    def __init__(self, vertex: Vertex):
        self.edges = dict()
        self.state_count = 0
        self.vertex = vertex

    # Get amount of times vertex appear in particular state
    def get_state_count(self) -> int :
        return self.state_count
    
    # Increment each time vertex is in particular state
    def increment_state_count(self):
        self.state_count += 1

    # Get info on all switch pairs vertex is part of under particular state
    def get_all_edges_info(self) -> list[tuple[str, int]]:
        edge_info: list[tuple[str, int]] = list()
        for (other, edge) in self.edges.items():
            edge_info.append((other, edge.get_count()))
        return edge_info
    
    @abstractmethod
    def get_adjacent_vertices(self) -> dict[str, Vertex] :
        pass

    @abstractmethod   
    def get_edges(self) -> dict[(str, str), Edge]:
        pass

class StateA(State):
    def __init__(self, vertex: Vertex):
        super().__init__(vertex)
            
    # Make and get switch pair edge to new label when switching from vertex's label
    def get_edge_to_b(self, vertex_b: Vertex) -> Edge :
        """Only should be called from function "get_edge_from_a" in Vertex object"""
        self.state_count += 1
        state_b_name = vertex_b.get_name()
        if state_b_name not in self.edges.keys():
            self.edges[state_b_name] = Edge(self.vertex, vertex_b)
        return(self.edges[state_b_name])
    
    def get_adjacent_vertices(self) -> dict[str, Vertex]:
        return {e.get_vertex_b().get_name(): e.get_vertex_b() for e in self.edges.values()}
    
    def get_edges(self) -> dict[(str, str), Edge]:
        vertex_a = self.vertex.get_name()
        edges: dict[(str, str), Edge] = dict()
        for edge in self.edges.items(): 
            edges[(vertex_a, edge[0])] = edge[1]
        return edges

class StateB(State):
    def __init__(self, vertex: Vertex):
        super().__init__(vertex)

    # Add switch pair edge from old label when switching to vertex's label
    def add_edge_from_a(self, vertex_a: Vertex, edge: Edge):
        """Only should be called from function "add_edge_to_b" in Vertex object"""
        self.state_count += 1
        state_a_name = vertex_a.get_name()
        if state_a_name not in self.edges.keys():
            self.edges[state_a_name] = edge

    def get_adjacent_vertices(self) -> dict[str, Vertex]:
        return {e.get_vertex_a().get_name(): e.get_vertex_a() for e in self.edges.values()}
    
    def get_edges(self) -> dict[(str, str), Edge]:
        vertex_b = self.vertex.get_name()
        edges: dict[(str, str), Edge] = dict()
        for edge in self.edges.items(): 
            if edge[0] != vertex_b:
                edges[(edge[0], vertex_b)] = edge[1]
        return edges

# Old class, will need to remove
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
    
# Class for amr label
class AmrVertex(Vertex):
    def __init__(self, amr: str, ref_count: int):
        self.amr_cl = amr
        self.ref_count = ref_count
        self.state_a = StateA(self)
        self.state_b = StateB(self)

    def get_name(self) -> str:
        return self.amr_cl
    
    def has_arg_label(self) -> bool:
        return False

# Class for cluster|amr label
class ClstrVertex(Vertex):
    clstr: str

    def __init__(self, id: str, ref_count: int):
        (self.clstr, self.amr_cl) = tuple(id.split("|"))
        self.ref_count = ref_count
        self.state_a = StateA(self)
        self.state_b = StateB(self)

    def get_name(self) -> str:
        return "|".join((self.clstr, self.amr_cl))
    
    def has_arg_label(self) -> bool:
        return False

# Class for domain|amr label
class DomainVertex(Vertex):
    dom_acc: str
    is_arg_label: bool

    def __init__(self, id: str, ref_count: int, is_arg_label: bool):
        (self.dom_acc, self.amr_cl) = tuple(id.split("|"))
        self.ref_count = ref_count
        self.state_a = StateA(self)
        self.state_b = StateB(self)
        self.is_arg_label = is_arg_label

    def get_name(self) -> str:
        return "|".join((self.dom_acc, self.amr_cl))
    
    def has_arg_label(self) -> bool:
        return self.is_arg_label

# Class for superclass|amr label
class SuperVertex(Vertex):
    super_acc: str
    is_arg_label: bool

    def __init__(self, id: str, ref_count: int, is_arg_label: bool):
        (self.super_acc, self.amr_cl) = tuple(id.split("|"))
        self.ref_count = ref_count
        self.state_a = StateA(self)
        self.state_b = StateB(self)
        self.is_arg_label = is_arg_label

    def get_name(self) -> str:
        return "|".join((self.super_acc, self.amr_cl))
    
    def has_arg_label(self) -> bool:
        return self.is_arg_label