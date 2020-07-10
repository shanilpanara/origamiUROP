import re
from typing import List

import numpy as np
import meshio
import matplotlib.pyplot as plt
import matplotlib.tri as tri

# class Vertex:
#     """
#     A Vertex is a point in space that is connected to other vertices

#     Defined simply by: x, y and z (in 3D)
#     """

#     def __init__(self, position: np.ndarray):
#         self.position = position
#         self.x = position[:, 0]
#         self.y = position[:, 1]

#     def distance(self, another_vertex: "Vertex") -> float:
#         return np.linalg.norm(self.position - another_vertex.position)

EDGE_TYPE = {0: "boundary", 1: "scaffold", 2: "staple"}

class Edge():
    """ 
    An Edge is a combination of two vertices and stores geometric information about this connection 
    """

    def __init__(self, vertex_1: list, vertex_2: list, edge_type: int = 0):
        self.vertices = np.array([vertex_1, vertex_2])

        # don't use type because that's a protected function in Python
        self.kind = EDGE_TYPE[edge_type]
        
        # these belong to individual nucleotides not the whole edge
        # self.linear_velocity = 0
        # self.angular_velocity = 0

    @property
    def length(self):
        return np.linalg.norm(self.vertices[1] - self.vertices[0])

    @property
    def vector(self) -> np.ndarray:
        return self.vertices[1] - self.vertices[0]

    @property
    def unit_vector(self) -> np.ndarray:
        return self.vector / (self.vector ** 2).sum() ** 0.5


def define_edges(vertices_of_polygon, index_1, index_2, edge_type):
    return Edge(vertices_of_polygon[index_1, :], vertices_of_polygon[index_2, :],)


class BoundaryPolygon():
    """
    A Polygon is a cyclic combination of 3 or more edges

    The first vertex on the starting edge = the last vertex on the final edge
    """

    def __init__(self, vertices: np.ndarray):
        """ 4 or more vertex objects given, each value on a given row """

        assert vertices.shape[1] == 3 or vertices.shape[1] == 2
        self.vertices = vertices

    @property
    def n_edges(self) -> np.array:
        return self.vertices.shape[0]

    @property
    def edges(self) -> List[Edge]:
        """
        Iterate over all vertices and return a list of Edge instances
        """
        edges = []
        for i in range(self.n_edges):
            edges.append(
                Edge(self.vertices[i-1], self.vertices[i], 0)
            )
        return edges

    @property
    def x(self) -> np.ndarray:
        return self.vertices[:, 0]

    @property
    def y(self) -> np.ndarray:
        return self.vertices[:, 1]

    @property
    def z(self) -> np.ndarray:
        try:
            return self.vertices[:, 2]
        except IndexError as err:
            raise err(f'Trying to access Polygon.z but {self} is only 2D!')

    def __repr__(self) -> str:
        return f"<Polygon{self.vertices.shape[1]}D Vertices[{self.vertices.shape[0]}]>"

    def __str__(self) -> str:
        return "This polygon has {} edges".format(self.n_edges)

    def write_STL(self, fout : str):
        triangulation = tri.Triangulation(self.vertices[:, 0], self.vertices[:, 1])
        triangles = triangulation.get_masked_triangles()
        cells = [("triangle", triangles)]
        meshio.write_points_cells(fout, self.vertices, cells)

    def write_PLY(self, fout : str, comments : List[str] = []):
        """
        Writes the BoundaryPolygon to a ASCII PLY File.
        """

        # begin header
        output_string = "PLY\nformat ascii 1.0\n"
        for comment in comments:
            output_string += f'comment {comment}\n'
        output_string += f"element vertex {self.vertices.shape[0]}\n"
        output_string += ''.join(
                [f"property float {i}\n" for i in ['x', 'y', 'z']]
            )
        output_string += \
            "element face 1\nproperty list uchar int vertex_index\nend_header\n"
        # end of header

        # vertex position list
        # regex substitutes out '[' and ']' characters and the replace
        # removes the space at the beginning of each new line
        output_string += \
            re.sub(r"\[|\]", '', self.vertices.__str__()).replace('\n ', '\n')

        # only one face in face list e.g. for a square 4 0 1 2 3
        output_string += f"\n{self.vertices.shape[1]} "
        output_string += ' '.join([str(i) for i in range(self.vertices.shape[1])])
        output_string += '\n'

        # finished generating file so now writing
        with open(fout, 'w') as f:
            f.write(output_string)

    def plot2D(
            self, 
            ax : plt.Axes = None, 
            fout : str = None, 
            show : bool = True,
            **kwargs
        ):
        """
        Assumes that the shape is 2D and lies on the z=0 plane.

        """
        if not ax:
            fig, ax = plt.subplots()
        ax.plot(self.x, self.y, 'k-', **kwargs)
        ax.set_aspect("equal", "datalim")
        if fout:
            fig.savefig(fout)
        if show:
            fig.show()
            
# ---plot---#
def plot_the_vertices(vertices: np.ndarray):
    fig = plt.figure()
    nodes = vertices
    triangulation = tri.Triangulation(nodes[:, 0], nodes[:, 1])
    hello = triangulation.get_masked_triangles()
    plt.triplot(triangulation, "-k")
    plt.axes().set_aspect("equal", "datalim")
    fig.show()

    fig2 = plt.figure()
    new_nodes_x = np.append(nodes[:, 0], nodes[0, 0])
    new_nodes_y = np.append(nodes[:, 1], nodes[0, 1])
    print(new_nodes_x, new_nodes_y)
    plt.plot(new_nodes_x, new_nodes_y)
    plt.axes().set_aspect("equal", "datalim")
    fig2.show()

    print(SHAPE.x)
    return


# ---run the script---#
def make_polygon(polygon_vertices: np.ndarray):

    plot_the_vertices(polygon_vertices)

    i = 1
    # Current Functionalities
    print(f"This shape has {SHAPE.edges.__len__()} edges \n")
    print(
        f"The coordinates of edge number {i+1} are: {str(SHAPE.edges[i].vertices[0])} and {str(SHAPE.edges[i].vertices[1])} \n"
    )
    print(f"It has a size of {round(SHAPE.edges[i].__len__(),4)}\n")
    print(f"This edge is a type of {SHAPE.edges[i].type} edge \n")
    print(f"It's vector is {SHAPE.edges[i].vector}\n")
    print(f"It's unit vector is {SHAPE.edges[i].unit_vector}")
    return
