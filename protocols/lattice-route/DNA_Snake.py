from origamiUROP.oxdna import System
from origamiUROP.lattice import LatticeRoute
from origamiUROP.polygons import BoundaryPolygon
from matplotlib import pyplot as plt

import numpy as np

from os import path

ROOT = "/".join(path.abspath(__file__).split("/")[:-1])


def generate(polygon_vertices: np.ndarray, DNAout: str = None, PLOTout: str = "connected"):
    polygon = BoundaryPolygon(polygon_vertices)
    lattice = polygon.lattice(straightening_factor=5, start_side="left")
    lattice.plot(lattice.final_coords, fout=PLOTout, poly = False, lattice_points= True)

    route = lattice.route()
    route.plot()

    if DNAout:
        system= route.system()
        system.write_oxDNA(prefix = DNAout, root=ROOT)

    return lattice.get_crossovers_coords()

if __name__ == "__main__":
    square = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 4.0, 0.0], [0.0, 4.0, 0.0]])
    trap = np.array([[0.,0.,0.],[1.5,6.,0.],[8.5,6.,0.],[10.,0.,0.]])
    trapREV = np.array([[0.,10.,0.],[2.5,4.,0.],[7.5,4.,0.],[10.,10.,0.]])
    no = 0.85
    hexagon = np.array(
        [[0, 1, 0], [no, 0.5, 0], [no, -0.5, 0], [0, -1, 0], [-no, -0.5, 0], [-no, 0.5, 0]])  
    plus = np.array([
        [1.,0.,0.], [2.,0.,0.], [2.,1.,0.], [3.,1.,0.], [3.,2.,0.], [2.,2.,0.],
        [2.,3.,0.], [1.,3.,0.], [1.,2.,0.], [0.,2.,0.], [0.,1.,0.], [1.,1.,0.]])
    diamond = np.array([[1.,0.,0.],[2.,1.,0.],[1.,2.,0.],[0.,1.,0.]])
    triangle = np.array([[0,0,0],[5,9,0],[10,0,0]])


    lattice = generate(square, DNAout=None, PLOTout="square_con")
    lattice = generate(trap*2, DNAout=None, PLOTout="trap_con")
    lattice = generate(hexagon*18, DNAout=None, PLOTout="hex_con")
    lattice = generate(plus*6, DNAout=None, PLOTout="plus_con")
    lattice = generate(diamond*20, DNAout=None, PLOTout="diam_con")
    lattice = generate(trapREV*5, DNAout=None,PLOTout="trapREV_con")
    lattice = generate(triangle*4, DNAout=None,PLOTout="triangle")

