from differential_distribution_classes.graph_items import Vertex, Edge
from differential_distribution_classes.math_util import relative_abundance
import pandas as pd
import numpy as np

# SubGraph includes a maximum connected subgraph from graph
class SubGraph:
    subgraph_index: int
    vertices: dict[str, Vertex]
    edges: dict[(str, str), Edge]
    switch_pair_counts: np.array
    reverse_switch_pair_counts: np.array
    dia_ref_abundance: np.array
    dia_a_abundance: np.array
    dia_b_abundance: np.array
    dee_ref_abundance: np.array
    dee_a_abundance: np.array
    dee_b_abundance: np.array
    abundance_i: np.array
    abundance_a: np.array
    abundance_b: np.array
    relative_i: np.array
    relative_a: np.array
    relative_b: np.array
    a_to_b_sign: np.array

    def __init__(self, index: int, vertex_entry: tuple[str, Vertex]):
        self.subgraph_index = index
        self.vertices = dict([vertex_entry])
        self.edges = vertex_entry[1].get_edges()

    def add_vertex(self, vertex_entry: tuple[str, Vertex]):
        self.vertices.update([vertex_entry])
        self.edges.update(vertex_entry[1].get_edges())

    def get_size(self):
        return len(self.vertices)

    # Those are raw abundance of vertex at each state
    def calc_abundance(self):
        self.abundance_i = np.array([
            ver.get_state_i_count() for ver in self.vertices.values()])
        self.abundance_a = np.array([
            ver.get_state_a_count() for ver in self.vertices.values()])
        self.abundance_b = np.array([
            ver.get_state_b_count() for ver in self.vertices.values()])
        self.a_to_b_sign = np.sign([
            ver.get_state_b_count() - ver.get_state_a_count() for ver in self.vertices.values()])
    
    # Those are abundances of vertex at each state relative to rest of connected subgraph
    def calc_relative_abundance(self):
        self.relative_i = relative_abundance(np.array([
            ver.get_state_i_count() for ver in self.vertices.values()]))
        self.relative_a = relative_abundance(np.array([
            ver.get_state_a_count() for ver in self.vertices.values()]))
        self.relative_b = relative_abundance(np.array([
            ver.get_state_b_count() for ver in self.vertices.values()]))
        
    def get_connected_table(self, sample, identity, model) -> pd.DataFrame:
        return (pd.DataFrame({
            "sample" : sample,
            "alignment identity" : identity,
            "model" : model,
            "label" : self.vertices.keys(),
            "subgraph" : self.subgraph_index,
            "reference abundance" : np.int_(self.abundance_i),
            "reference relative abundance" : self.relative_i,
            "state A abundanc" : np.int_(self.abundance_a),
            "state A relative abundance" : self.relative_a,
            "state B abundance" : np.int_(self.abundance_b),
            "state B relative abundance" : self.relative_b,
            "a_to_b_sign": np.int_(self.a_to_b_sign),
            "label is arg": [ver.has_arg_label() for ver in self.vertices.values()]}))
    
    # Those are abundance information for each switch pair
    def calc_pair(self):
        self.switch_pair_counts = np.array([
            edge.get_count() for edge in list(self.edges.values())])
        self.reverse_switch_pair_counts = np.array([
            0 if (dee, dia) not in self.edges else self.edges[(dee, dia)].get_count()
            for dia, dee in list(self.edges.keys())])
        self.dia_ref_abundance = np.array([
            np.int_(self.vertices[ver[0]].get_state_i_count()) for ver in list(self.edges.keys())])
        self.dia_a_abundance = np.array([
            np.int_(self.vertices[ver[0]].get_state_a_count()) for ver in list(self.edges.keys())])
        self.dia_b_abundance = np.array([
            np.int_(self.vertices[ver[0]].get_state_b_count()) for ver in list(self.edges.keys())])
        self.dee_ref_abundance = np.array([
            np.int_(self.vertices[ver[1]].get_state_i_count()) for ver in list(self.edges.keys())])
        self.dee_a_abundance = np.array([
            np.int_(self.vertices[ver[1]].get_state_a_count()) for ver in list(self.edges.keys())])
        self.dee_b_abundance = np.array([
            np.int_(self.vertices[ver[1]].get_state_b_count()) for ver in list(self.edges.keys())])
    
    def get_pair_table(self, sample, identity, model) -> pd.DataFrame:
        return (pd.DataFrame({
            "sample" : sample,
            "alignment identity" : identity,
            "model" : model,
            "diamond best-hit label": [ver[0] for ver in list(self.edges.keys())],
            "diamond best-hit reference abundance" : self.dia_ref_abundance,
            "diamond best-hit state a abundance" : self.dia_a_abundance,
            "diamond best-hit state b abundance" : self.dia_b_abundance,
            "deeparg hit label": [ver[1] for ver in list(self.edges.keys())],
            "deeparg hit reference abundance" : self.dee_ref_abundance,
            "deeparg hit state a abundance" : self.dee_a_abundance,
            "deeparg hit state b abundance" : self.dee_b_abundance,
            "switch pair count": self.switch_pair_counts,
            "reverse switch pair count": self.reverse_switch_pair_counts,
            "pair contains arg label": [edge.has_arg_label for edge in self.edges.values()]}))

class Graph:
    subgraphs: list[SubGraph]
    def __init__(self, all_vertices: dict[str, Vertex]):
        self.subgraphs = list()
        while len(all_vertices) > 0:
            last_vertex = all_vertices.popitem()
            self.subgraphs.append(SubGraph(len(self.subgraphs), last_vertex))
            vertex_queue = list(last_vertex[1].get_adjacent_vertices().items())
            while len(vertex_queue) > 0:
                curr_vertex = vertex_queue.pop(0)
                if all_vertices.get(curr_vertex[0]) is None:
                    continue
                all_vertices.pop(curr_vertex[0])

                self.subgraphs[-1].add_vertex(curr_vertex)
                vertex_queue.extend(curr_vertex[1].get_adjacent_vertices().items())
            self.subgraphs[-1].calc_abundance()
            self.subgraphs[-1].calc_relative_abundance()
            self.subgraphs[-1].calc_pair()
    
    def get_connected_table(self, sample, identity, model) -> pd.DataFrame :
        return(pd.concat([
            subgraph.get_connected_table(sample, identity, model) for subgraph in self.subgraphs 
            if subgraph.get_size() > 1]))

    def get_pair_table(self, sample, identity, model) -> pd.DataFrame :
        return(pd.concat([
            subgraph.get_pair_table(sample, identity, model) for subgraph in self.subgraphs 
            if subgraph.get_size() > 1]))
    