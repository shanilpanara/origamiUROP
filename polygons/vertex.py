import numpy as np


class Vertex:
    """
    A Vertex is a point in space that is connected to other vertices
    """

    def __init__(self, position: np.ndarray):
        self.position = position

    def distance(self, vertex: Vertex) -> float:
        return np.linalg.norm(self.position - vertex.position)


class Edge:
    """ 
    An Edge connects two vertices and stores geometric information about this connection 
    """

    def __init__(self, vertex_1: Vertex, vertex_2: Vertex):
        self.vertices = (vertex_1, vertex_2)

    def __len__(self):
        return self.vertices[0].distance(self.vertices[1])

    @property
    def vector(self) -> np.ndarray:
        return self.vertices[1] - self.vertices[0]
