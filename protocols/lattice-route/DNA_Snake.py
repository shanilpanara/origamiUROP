from origamiUROP.oxdna import System
from origamiUROP.lattice import LatticeRoute
from origamiUROP.polygons import BoundaryPolygon
from matplotlib import pyplot as plt

import numpy as np

from os import path

ROOT = "/".join(path.abspath(__file__).split("/")[:-1])


def generate(polygon_vertices: np.ndarray):
    polygon = BoundaryPolygon(polygon_vertices)
    lattice = polygon.lattice()
    route = lattice.route()

    system= route.system()
    system.write_oxDNA(root=ROOT)

    route.plot()

    return


def generate_square():
    square = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 4.0, 0.0], [0.0, 4.0, 0.0]])
    generate(square)

if __name__ == "__main__":
    generate_square()
