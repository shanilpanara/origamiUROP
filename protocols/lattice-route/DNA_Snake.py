from origamiUROP.oxdna import System
from origamiUROP.lattice import LatticeRoute
from origamiUROP.polygons import BoundaryPolygon

import numpy as np


def generate(polygon_vertices: np.ndarray):
    polygon = BoundaryPolygon(polygon_vertices)
    lattice = polygon.lattice()

    nodes = lattice.route()

    route = LatticeRoute(nodes)

    system = System(np.array([50.0, 50.0, 50.0]))
    system.add_strand(route)

    route.plot()
    system.write_oxDNA()


def generate_square():
    square = np.array([[0, 0, 0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
    generate(square * 10)


if __name__ == "__main__":
    generate_square()
