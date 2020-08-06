from origamiUROP.oxdna import System
from origamiUROP.lattice import LatticeRoute
from origamiUROP.polygons import BoundaryPolygon

import numpy as np


def test_square():
    square = np.array([[0, 0, 0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
    polygon = BoundaryPolygon(square * 10)
    lattice = polygon.lattice()
    nodes = lattice.route()
    route = LatticeRoute(nodes)
    #route.plot()
    system = System(np.array([50.0, 50.0, 50.0]))
    system.add_strand(route)
    system.write_oxDNA()


if __name__ == "__main__":
    test_square()
