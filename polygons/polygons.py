import numpy as np


class Vertex:
    """
    A Vertex is a point in space that is connected to other vertices
    
    Defined simply by: x, y and z (in 3D)
    """

    def __init__(self, position: np.ndarray):
        self.position = position
        self.x = position[:, 0]
        self.y = position[:, 1]

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

    def __init__(self, verticies: np.ndarray):
        """
        4 or more vertex objects given, each value on a given row
        
        What to do:
        - Edges (using the above class)
        - Interior angles
            - use cosine formula with vectors 
        """
        self.n_sides = self.edges.shape[0]  # no. of sides = no of rows
        self.interior_angles = interiorAngles()

        @property
        def interiorAngles(self) -> np.ndarray:  # containing floats
            return

