from origamiUROP.oxdna import System
from origamiUROP.lattice import LatticeRoute
from origamiUROP.polygons import BoundaryPolygon
from matplotlib import pyplot as plt

import numpy as np

from os import path

ROOT = "/".join(path.abspath(__file__).split("/")[:-1])


def generate(polygon_vertices: np.ndarray, DNAout: str = None, PLOTout: str = "connected"):
    polygon = BoundaryPolygon(polygon_vertices)
    lattice = polygon.lattice(straightening_factor=5)
    lattice.plot(lattice.final_coords, fout=PLOTout)
    # lattice.plot(lattice.final_array, fout="straightedge")

    # route = lattice.route()
    # route.plot()

    if DNAout:
        system= route.system()
        system.write_oxDNA(prefix = DNAout, root=ROOT)

    return lattice.get_crossovers_coords()

if __name__ == "__main__":
    square = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 4.0, 0.0], [0.0, 4.0, 0.0]])
    trap = np.array([[0.,0.,0.],[1.5,6.,0.],[8.5,6.,0.],[10.,0.,0.]])
    no = 0.85
    hexagon = np.array(
        [[0, 1, 0], [no, 0.5, 0], [no, -0.5, 0], [0, -1, 0], [-no, -0.5, 0], [-no, 0.5, 0]])  
    plus = np.array([
        [1.,0.,0.], [2.,0.,0.], [2.,1.,0.], [3.,1.,0.], [3.,2.,0.], [2.,2.,0.],
        [2.,3.,0.], [1.,3.,0.], [1.,2.,0.], [0.,2.,0.], [0.,1.,0.], [1.,1.,0.]])
    diamond = np.array([[1.,0.,0.],[2.,1.,0.],[1.,2.,0.],[0.,1.,0.]])


    lattice = generate(square, DNAout=None, PLOTout="square_con")
    lattice = generate(trap*2, DNAout=None, PLOTout="trap_con")
    lattice = generate(hexagon*10, DNAout=None, PLOTout="hex_con")
    lattice = generate(plus*10, DNAout=None, PLOTout="plus_con")
    lattice = generate(diamond*[10,5,0], DNAout=None, PLOTout="diam_con")

