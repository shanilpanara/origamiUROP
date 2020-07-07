import numpy as np


class Vertex:
    """
    A Vertex is a point in space that is connected to other vertices
    
    Defined simply by: x, y and z (in 3D)
    """

    def __init__(self, position: np.ndarray):
        self.position = position

    def distance(self, another_vertex: Vertex) -> float:
        return np.linalg.norm(self.position - another_vertex.position)


class Edge:
    """ 
    An Edge is a combination of two vertices and stores geometric information about this connection 
    """

    def __init__(self, vertex_1: Vertex, vertex_2: Vertex):
        self.vertices = (vertex_1, vertex_2)

    def __len__(self):
        return self.vertices[0].distance(self.vertices[1])

    @property
    def vector(self) -> np.ndarray:
        return self.vertices[1] - self.vertices[0]


class Polygon:
    """
    A Polygon is a cyclic combination of 3 or more edges

    The first vertex on the starting edge = the last vertex on the final edge
    """

    def __init__(self, edges: np.ndarray):
        self.edges = edges
        self.n_sides = edges.shape[0]  # no of rows
