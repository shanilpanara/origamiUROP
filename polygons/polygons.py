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


class Edge:
    """ 
    An Edge is a combination of two vertices and stores geometric information about this connection 
    """

    def __init__(self, vertex_1: list, vertex_2: list, edge_type: int = 0):
        self.vertices = np.array([vertex_1, vertex_2])
        self.type = EDGE_TYPE[edge_type]
        ## DEFAULT INITIALISED VALUES
        self.dir = np.array([0.0, 0.0, 1.0])  # direction
        self.linear_velocity = 0.0
        self.angular_velocity = 0.0

    def __len__(self):
        return np.linalg.norm(self.vertices[1] - self.vertices[0])

    @property
    def vector(self) -> np.ndarray:
        return self.vertices[1] - self.vertices[0]

    @property
    def unit_vector(self):
        return self.vector / (self.vector ** 2).sum() ** 0.5


def define_edges(vertices_of_polygon, index_1, index_2, edge_type):
    return Edge(vertices_of_polygon[index_1, :], vertices_of_polygon[index_2, :],)


class BoundaryPolygon:
    """
    A Polygon is a cyclic combination of 3 or more edges

    The first vertex on the starting edge = the last vertex on the final edge
    """

    def __init__(self, vertices: np.ndarray):
        """ 4 or more vertex objects given, each value on a given row """

        self.vertices = vertices
        self.x = vertices[:, 0]
        self.y = vertices[:, 1]
        self.z = vertices[:, 2]
        self.n_edges = self.vertices.shape[0]  # no. of sides = no of rows

        self.edges = []

        # makes a index which "wraps around on itself, e.g. 4 edges -> [0,1,2,3,0]"
        counter = np.append(np.arange(0, self.n_edges), 0)
        for i in range(self.n_edges):
            self.edges.append(
                define_edges(self.vertices, counter[i], counter[i + 1], 0)
            )

    def __str__(self):
        return "This polygon has {} edges".format(self.n_edges)

    def write_STL(self, filename: str):
        triangulation = tri.Triangulation(self.vertices[:, 0], self.vertices[:, 1])
        triangles = triangulation.get_masked_triangles()
        cells = [("triangle", triangles)]
        filename = filename + ".stl"
        print(filename)
        meshio.write_points_cells(filename, self.vertices, cells)


# ---plot---#
def plot_the_vertices(vertices: np.ndarray):
    fig = plt.figure()
    nodes = vertices
    triangulation = tri.Triangulation(nodes[:, 0], nodes[:, 1])
    hello = triangulation.get_masked_triangles()
    plt.triplot(triangulation, "-k")
    plt.axes().set_aspect("equal", "datalim")
    plt.show()
    return hello


# ---run the script---#
def make_polygon(polygon_vertices: np.ndarray):
    shape_object = BoundaryPolygon(polygon_vertices)

    hello = plot_the_vertices(polygon_vertices)

    i = 1
    # Current Functionalities
    print(f"This shape has {shape_object.edges.__len__()} edges \n")
    print(
        f"The coordinates of edge number {i+1} are: {str(shape_object.edges[i].vertices[0])} and {str(shape_object.edges[i].vertices[1])} \n"
    )
    print(f"It has a size of {round(shape_object.edges[i].__len__(),4)}\n")
    print(f"This edge is a type of {shape_object.edges[i].type} edge \n")
    print(f"It's vector is {shape_object.edges[i].vector}\n")
    print(f"It's unit vector is {shape_object.edges[i].unit_vector}")
    return shape_object, hello


###---stored arrays of some shapes---###
square = np.array([[0, 0, 0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
trapezium = np.array([[0, 0, 0], [1.0, 10.0, 0.0], [11.0, 10.0, 0.0], [15.0, 0.0, 0.0]])

no = 0.866025
hexagon = np.array(
    [[0, 1, 0], [no, 0.5, 0], [no, -0.5, 0], [0, -1, 0], [-no, -0.5, 0], [-no, 0.5, 0]]
)


###---run the program---###
SHAPE, hello = make_polygon(hexagon)
