from origamiUROP.oxdna import System
from origamiUROP.lattice import LatticeRoute, DNASnake
from origamiUROP.polygons import BoundaryPolygon

import numpy as np

from os import path

ROOT = "/".join(path.abspath(__file__).split("/")[:-1])


def generate(polygon_vertices: np.ndarray, DNAout: str = None, PLOTout: str = None):
    polygon = BoundaryPolygon(polygon_vertices)
    lattice = polygon.dna_snake(straightening_factor=5, start_side="left")
    route = lattice.route()

    if PLOTout:
        plot_list = [lattice.quantised_array, lattice.crossover_coords]
        lattice.plot(plot_list, fout=PLOTout, poly = True)
    route.plot()

    if DNAout:
        system = route.system()
        system.write_oxDNA(prefix = DNAout, root=ROOT)


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


    # generate(square, DNAout="square", PLOTout="squares")
    # generate(trap*8, DNAout="trapezium", PLOTout="trapeziums")
    # generate(hexagon*18, DNAout="hexagon", PLOTout="hexagons")
    # generate(plus*6, DNAout="plus", PLOTout="plus")
    # generate(diamond*[20,30,0], DNAout="diamond", PLOTout="diamond")
    # generate(trapREV*5, DNAout="trapezium_rev",PLOTout="trapezium_rev")
    generate(triangle*6, DNAout="triangle",PLOTout="triangle")

    # generate(plus*6, DNAout="plus", PLOTout="plus")

