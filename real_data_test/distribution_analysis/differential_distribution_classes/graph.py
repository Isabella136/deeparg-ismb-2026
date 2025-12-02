from differential_distribution_classes.graph_items import Vertex
from differential_distribution_classes.math_util import array_wise_clr_transform, relative_abundance
import pandas as pd
import numpy as np

# SubGraph includes a maximum connected subgraph from graph
class SubGraph:
    subgraph_index: int
    vertices: dict[str, Vertex]
    clr_transform_i: np.array
    clr_transform_a: np.array
    clr_transform_b: np.array
    relative_i: np.array
    relative_a: np.array
    relative_b: np.array
    a_to_b_sign: np.array

    def __init__(self, index: int, vertex_entry: tuple[str, Vertex]):
        self.subgraph_index = index
        self.vertices = dict([vertex_entry])

    def add_vertex(self, vertex_entry: tuple[str, Vertex]):
        self.vertices.update([vertex_entry])

    def calc_clr_transform(self):
        self.clr_transform_i = array_wise_clr_transform(np.array([
            ver.get_state_i_count() for ver in self.vertices.values()]))
        self.clr_transform_a = array_wise_clr_transform(np.array([
            ver.get_state_a_count() for ver in self.vertices.values()]))
        self.clr_transform_b = array_wise_clr_transform(np.array([
            ver.get_state_b_count() for ver in self.vertices.values()]))
        self.a_to_b_sign = np.sign([
            ver.get_state_b_count() - ver.get_state_a_count() for ver in self.vertices.values()])
        
    def calc_relative_abundance(self):
        self.relative_i = relative_abundance(np.array([
            ver.get_state_i_count() for ver in self.vertices.values()]))
        self.relative_a = relative_abundance(np.array([
            ver.get_state_a_count() for ver in self.vertices.values()]))
        self.relative_b = relative_abundance(np.array([
            ver.get_state_b_count() for ver in self.vertices.values()]))
        
    def get_clr_transform(self) -> pd.DataFrame:
        return (pd.DataFrame({
                "subgraph" : np.full(shape=self.clr_transform_i.shape, fill_value=self.subgraph_index),
                "state I clr" : np.sign(self.clr_transform_i) * np.abs(self.clr_transform_i),
                "state A clr" : np.sign(self.clr_transform_a) * np.abs(self.clr_transform_a),
                "state B clr" : np.sign(self.clr_transform_b) * np.abs(self.clr_transform_b), 
                "state B - state A sign" : self.a_to_b_sign},
            index=self.vertices.keys()))
    def get_relative_abundance(self):
        return(pd.DataFrame({
            "subgraph" : np.full(shape=self.clr_transform_i.shape, fill_value=self.subgraph_index),
            "relative abundance I" : self.relative_i,
            "relative abundance A" : self.relative_a,
            "relative abundance B" : self.relative_b},
            index=self.vertices.keys()))

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
            self.subgraphs[-1].calc_clr_transform()
            self.subgraphs[-1].calc_relative_abundance()
    
    def get_connected_clr_transform(self) -> pd.DataFrame :
        return(pd.concat([subgraph.get_clr_transform() for subgraph in self.subgraphs]))
    
    
    def get_relative_abundance(self) -> pd.DataFrame:
        return (pd.concat([subgraph.get_relative_abundance() for subgraph in self.subgraphs]))
    