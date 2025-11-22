from differential_distribution_classes.graph_items import Vertex
from differential_distribution_classes.math_util import clr_transform
import pandas as pd
import numpy as np

class SubGraph:
    subgraph_index: int
    vertices: dict[str, Vertex]
    clr_transform_i: np.array
    clr_transform_a: np.array
    clr_transform_b: np.array
    a_to_b_sign: np.array

    def __init__(self, index: int, vertex_entry: tuple[str, Vertex]):
        self.subgraph_index = index
        self.vertices = dict([vertex_entry])

    def add_vertex(self, vertex_entry: tuple[str, Vertex]):
        self.vertices.update([vertex_entry])

    def calc_clr_transform(self):
        self.clr_transform_i = clr_transform(np.array([
            ver.get_state_i_count() for ver in self.vertices.values()]))
        self.clr_transform_a = clr_transform(np.array([
            ver.get_state_a_count() for ver in self.vertices.values()]))
        self.clr_transform_b = clr_transform(np.array([
            ver.get_state_b_count() for ver in self.vertices.values()]))
        self.a_to_b_sign = np.sign([
            ver.get_state_b_count() - ver.get_state_a_count() for ver in self.vertices.values()])
        
    def get_clr_transform(self) -> pd.DataFrame:
        return (pd.DataFrame({
                "subgraph" : np.full(shape=self.clr_transform_i.shape, fill_value=self.subgraph_index),
                "state I clr" : np.sign(self.clr_transform_i) * np.abs(self.clr_transform_i),
                "state A clr" : np.sign(self.clr_transform_a) * np.abs(self.clr_transform_a),
                "state B clr" : np.sign(self.clr_transform_b) * np.abs(self.clr_transform_b), 
                "state B - state A sign" : self.a_to_b_sign },
            index=self.vertices.keys()))


class Graph:
    subgraphs: list[SubGraph]

    def __init__(self, all_vertices: dict[str, Vertex]):
        self.subgraphs = list()
        while len(all_vertices) > 0:
            last_vertex = all_vertices.popitem()
            self.subgraphs.append(SubGraph(len(self.subgraphs), last_vertex))
            vertex_queue = last_vertex[1].get_adjacent_vertices()
            while len(vertex_queue) > 0:
                curr_vertex = vertex_queue.pop(0)
                curr_name = curr_vertex.get_name()
                if all_vertices.get(curr_name) is None:
                    continue
                all_vertices.pop(curr_name)
                self.subgraphs[-1].add_vertex((curr_name, curr_vertex))
                vertex_queue.extend(curr_vertex.get_adjacent_vertices())
            self.subgraphs[-1].calc_clr_transform()
    
    def get_clr_transform(self) -> pd.DataFrame :
        return(pd.concat([subgraph.get_clr_transform() for subgraph in self.subgraphs]))
    