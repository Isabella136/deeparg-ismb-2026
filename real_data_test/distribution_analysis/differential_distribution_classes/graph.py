from differential_distribution_classes.graph_items import Vertex

class SubGraph:
    vertices: dict[str, Vertex]
    geom_mean_i: float
    geom_mean_a: float
    geom_mean_a_to_b: float

    def __init__(self, vertex_entry: tuple[str, Vertex]):
        self.vertices = dict(vertex_entry)

    def add_vertex(self, vertex_entry: tuple[str, Vertex]):
        self.vertices.update(vertex_entry)


class Graph:
    subgraphs: list[SubGraph]

    def __init__(self, all_vertices: dict[str, Vertex]):
        self.subgraphs = list()
        while all_vertices.len() > 0:
            last_vertex = all_vertices.popitem()
            self.subgraphs.append(SubGraph(last_vertex))
            vertex_queue = last_vertex[1].get_adjacent_vertices()
            while vertex_queue.len() > 0:
                curr_vertex = vertex_queue.pop(0)
                curr_name = curr_vertex.get_name()
                if all_vertices.get(curr_name) is None:
                    continue
                all_vertices.pop(curr_name)
                self.subgraphs[-1].add_vertex((curr_name, curr_vertex))
                vertex_queue.extend(curr_vertex.get_adjacent_vertices)