from origamiUROP.oxdna import System
from origamiUROP.lattice import LatticeRoute
from origamiUROP.polygons import BoundaryPolygon

import numpy as np

from os import path

ROOT = "/".join(path.abspath(__file__).split("/")[:-1])


def generate(polygon_vertices: np.ndarray):
    polygon = BoundaryPolygon(polygon_vertices)
    lattice = polygon.lattice()

    nodes = lattice.route()

    route = LatticeRoute(nodes)
    system = route.system()
    system.add_strand(route)

    route.plot()
    system.write_oxDNA(root=ROOT)


def generate_square():
    square = np.array([[0, 0, 0], [10.0, 1.0, 0.0], [10.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
    generate(square)


if __name__ == "__main__":
    generate_square()
