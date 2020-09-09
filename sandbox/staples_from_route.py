from origamiUROP.lattice import LatticeRoute, LatticeEdge, LatticeNode
from origamiUROP.oxdna import Strand, System
from origamiUROP.tools import DNANode, DNAEdge
from origamiUROP.lattice.utils import find_crossover_locations
from origamiUROP.polygons import BoundaryPolygon

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from typing import List
from itertools import cycle
from copy import deepcopy

class StapleNode(DNANode):
    def __init__(self, position : np.ndarray):
        super().__init__(position)

class StapleEdge(DNAEdge):
    def __init__(self, vertex_1: StapleNode, vertex_2: StapleNode):
        super().__init__(vertex_1, vertex_2)

class StapleRoute(Strand):
    def __init__(self, scaffold_rows: List[Strand], nodes: List[StapleNode] = []):
        self._nodes = nodes
        self._scaffold_rows = scaffold_rows
        self.update_strand_and_nucleotides()
        super().__init__(self._nucleotides)


    @property
    def nodes(self) -> List[StapleNode]:
        return self._nodes

    @property
    def edges(self):
        _edges = [StapleEdge(node, self.nodes[i+1]) for i, node in enumerate(self.nodes[:-1])]
        return _edges
    
    @property
    def scaffold_rows(self) -> List[Strand]:
        """ Returns scaffold as it's constituent rows"""
        return self._scaffold_rows

    def update_strand_and_nucleotides(self, **kwargs):
        self._nucleotides = []
        for i, edge in enumerate(self.edges):
            if i % 2 == 1:
                continue

            x1 = int(edge.vertices[0][0])
            x2 = int(edge.vertices[1][0])
            row = int(edge.vertices[0][1])
            if len(self.scaffold_rows) - 1 == row:
                break
            nucleotides = []
            scaffold_row = self.scaffold_rows[row]
            
            if x1 > x2:
                x1, x2 = x2, x1
                for nucleotide in scaffold_row.nucleotides[x1:x2+1]:
                    nucleotides.append(nucleotide.make_across())
                # switch the order of the nucleotides back again
                nucleotides = nucleotides[::-1]
            else:
                for nucleotide in scaffold_row.nucleotides[::-1][x1:x2+1]:
                    nucleotides.append(nucleotide.make_across())
 
            for nucleotide in nucleotides:
                self._nucleotides.append(nucleotide)

class StapleCollection:
    """ Probably need to edit this such that
    it becomes easy to join staples on the same row together
    not sure, how to implement that yet """
    def __init__(self, strands: List[Strand] = [], route: LatticeRoute = None):
        self._strands = strands
        self.route = route
                
        # Import scaffold strands as its constituent rows
        self.scaffold_strands = route.get_strand()[1]
        # Import scaffold information class
        self.scaffold_rows = Scaffold(route)
        
        # Two copies of scaffold array. The stapled array will change as staples are added
        self.scaffold_array = self.scaffold_rows.array
        self.stapled_array = deepcopy(self.scaffold_array)
    
    @property
    def staples(self) -> List[Strand]:
        return self._strands

    @property
    def n_staples(self) -> int:
        return len(self._strands)

    def add_staples(self, staple_strand: Strand):
        self._strands.append(staple_strand)

    def generate_staples(self, route : LatticeRoute, staple_width = 16):
        # deal with start case

        # deal with end case

        # deal with cases where next row is bigger

        # deal with cases where next row is smaller

        # deal with the rest :)
        

    #     # get all the information from rows

        
    #     for row, info in scaffold.info.iterrows():
    #         x2 = info["bounds"][0]
    #         if info["start side"] == "left":
    #             x1 = x2 + staple_width
    #         elif info["start side"] == "right":
    #             x1 = x2 - staple_width
            
    #         staple_nodes = [
    #             StapleNode([x1, row, 0.0]),
    #             StapleNode([x2, row, 0.0]),
    #             StapleNode([x2, row+1, 0.0]),
    #             StapleNode([x1, row+1, 0.0])]

    #         staples.append(StapleRoute(scaffold_rows, staple_nodes))

    #     return StapleCollection(staples)
        pass

    def plot_nodes(self, strand: Strand, ax, colour = 'r', width = 0.01, **kwargs):
        nodes = np.array(strand.nodes)
        #plt.grid(True)
        ax.plot(nodes[:, 0], nodes[:, 1], 'k', ms = 0.5, alpha = 0)
        ax.xaxis.set_major_locator(MultipleLocator(2))
        ax.yaxis.set_major_locator(MultipleLocator(5))
        ax.set_xlabel("No. of nucleotides")
        ax.set_ylabel("No. of strands")
        for edge in strand.edges:
            ax.arrow(
                edge.vertices[0][0], 
                edge.vertices[0][1], 
                edge.vector[0], 
                edge.vector[1], 
                width=width,
                color = colour,
                head_width = 0.1,
                length_includes_head=True, **kwargs)
    

    def plot(self, route: LatticeRoute, fout: str = None, colour: str = None):
        fig, ax = plt.subplots()
        if not colour:
            colours = cycle(('r','b'))
        else:
            colours = cycle((colour))
        if route:
            self.plot_nodes(strand = route, ax = ax, colour = 'k', width = 0.05, alpha = 0.05)
        for staple in self.staples:

            self.plot_nodes(strand = staple, ax = ax, colour = next(colours))
        
        plt.gca().set_aspect(5)
        if fout:
            plt.savefig(fout, dpi=500)

        plt.show()

class Scaffold:
    def __init__(self,route: LatticeRoute, staple_widths = [5, 16, 26]):
        self.route = route
        self.edges = route.edges[0::2] #only rows which are horizontal
        self.staple_widths = staple_widths
        self.sizes = self.get_row_sizes()
        self.start_side = self.get_start_side()

        self.startx = self.get_startx()
        self.endx = self.get_endx()
        self.bounds = self.get_bounds()

        self.array = self.get_scaffold_array()


    def get_row_sizes(self) -> list:
        """ Returns length of each row """
        return [edge.nt_length + 1 for edge in self.edges]
    
    def get_startx(self) -> list:
        """ Returns x coord of start point on each row """
        return [int(edge.vertices[0][0]) for edge in self.edges]
    
    def get_endx(self) -> list:
        """ Returns x coord of end point on each row """
        return [int(edge.vertices[1][0]) for edge in self.edges]
    
    def get_bounds(self) -> list:
        """ Returns x coord of start and end point on each row """
        return [[int(edge.vertices[0][0]),int(edge.vertices[1][0])] for edge in self.edges]    
    
    def get_start_side(self) -> list:
        """ Returns left or right as the start side of each row"""
        nodes0 = self.route.nodes[0::2]
        nodes1 = self.route.nodes[1::2]
        return ["left" if node0[0]-node1[0] < 0 else "right" for (node0,node1) in zip(nodes0, nodes1)]

    def get_scaffold_array(self) -> np.ndarray:
        # Find max row width and no. of columns and create empty grid
        scaffold = np.zeros( (len(self.edges), max(self.sizes)) )
        for row in range(len(self.edges)):
            scaffold[row][min(self.bounds[row]):max(self.bounds[row])+1] = 1
        return scaffold


    def n_staples(self, staple_width = None):
        """ Returns no. of staples per row """
        if not staple_width: # default to 16
            staple_width = self.staple_widths[1]
        return [int(i/staple_width) for i in self.sizes]
    
    def unpaired_bases(self, staple_width = None):
        """ Returns no. of unpaired bases per row """
        if not staple_width: # default to 16
            staple_width = self.staple_widths[1]
        return [int(i % staple_width) for i in self.sizes]

    @property
    def info(self, staple_width = None):
        """
        Returns information for each row in the Lattice Route
        Specifically its size and start/end point
        And the number of staples 
        """
        return pd.DataFrame({
                "size": self.sizes,
                "start side": self.start_side,
                "bounds": self.bounds,
                "staples (5)": self.n_staples(self.staple_widths[0]),
                "staples (16)": self.n_staples(self.staple_widths[1]),
                "unpaired bases (5)": self.unpaired_bases(self.staple_widths[0]),
                "unpaired bases (16)": self.unpaired_bases(self.staple_widths[1]),
            })
        

# def side_staples(route : LatticeRoute, staple_width = 16):
#     # import scaffold as formed of its constituent rows
#     scaffold_rows = route.get_strand()[1]
#     scaffold = Scaffold(route)
#     staples = []

#     for row, info in scaffold.info.iterrows():
#         x2 = info["bounds"][0]
#         if info["start side"] == "left":
#             x1 = x2 + staple_width
#         elif info["start side"] == "right":
#             x1 = x2 - staple_width
        
#         staple_nodes = [
#             StapleNode([x1, row, 0.0]),
#             StapleNode([x2, row, 0.0]),
#             StapleNode([x2, row+1, 0.0]),
#             StapleNode([x1, row+1, 0.0])]

#         staples.append(StapleRoute(scaffold_rows, staple_nodes))

#     return StapleCollection(staples)

def test_StapleRoute(route: LatticeRoute):
    x1 = 0
    x2 = 15
    row = 0
    scaffold_rows = route.get_strand()[1]
    staple_nodes = [
        StapleNode([x1, row, 0.0]),
        StapleNode([x2, row, 0.0]),
        StapleNode([x2, row+1, 0.0]),
        StapleNode([x1, row+1, 0.0])]
    
    return StapleRoute(scaffold_rows, staple_nodes)

def generate_route(polygon_vertices: np.ndarray):
    polygon = BoundaryPolygon(polygon_vertices)
    lattice = polygon.dna_snake(straightening_factor=5, start_side="left")
    return lattice.route()

if __name__ == "__main__":
    half_turn_indices   = [4, 15, 25, 36, 46, 56, 67, 77, 88, 98, 109]
    staple_lengths      = [9, 31, 51, 73]

    route_vertices = [
        [0., 0.,0.],
        [56.,0.,0.],
        [56.,1.,0.],
        [0., 1.,0.],
        [0., 2.,0.],
        [56.,2.,0.],
        [56., 3.,0.],
        [0., 3.,0.],
        [0., 4.,0.],
        [56., 4.,0.],
        [56.,5.,0.],
        [0., 5.,0.],
        [0., 6.,0.],
        [56.,6.,0.],
        [56., 7.,0.],
        [0., 7.,0.],
        [0., 8.,0.],
        [88., 8.,0.],
        [88., 9.,0.],
        [0., 9.,0.],
        [0., 10.,0.],
        [88., 10.,0.],
        [88., 11.,0.],
        [0., 11.,0.],
        [0., 12.,0.],
        [56., 12.,0.],
        [56., 13.,0.],
        [0., 13.,0.],
        [0., 14.,0.],
        [56., 14.,0.],
        [56, 15, 0],
        [0,15,0]
    ]
    trapREV = np.array([[0.,10.,0.],[2.5,4.,0.],[7.5,4.,0.],[10.,10.,0.]])
    no = 0.85
    hexagon = np.array(
        [[0, 1, 0], [no, 0.5, 0], [no, -0.5, 0], [0, -1, 0], [-no, -0.5, 0], [-no, 0.5, 0]])  
    plus = np.array([
        [1.,0.,0.], [2.,0.,0.], [2.,1.,0.], [3.,1.,0.], [3.,2.,0.], [2.,2.,0.],
        [2.,3.,0.], [1.,3.,0.], [1.,2.,0.], [0.,2.,0.], [0.,1.,0.], [1.,1.,0.]])
    diamond = np.array([[1.,0.,0.],[2.,1.,0.],[1.,2.,0.],[0.,1.,0.]])
    triangle = np.array([[0,0,0],[5,9,0],[10,0,0]])
    stacked_I = np.array([
    [0.,0.,0.],[3.,0.,0.],[3.,1.,0.],[2.,1.,0.], [2.,2.,0.],[3.,2.,0.],
    [3.,3.,0.],[2.,3.,0.],[2.,4.,0.],[3.,4.,0.],[3.,5.,0.],[0.,5.,0.],[0.,4.,0.],[1.,4.,0.],
    [1.,3.,0.],[0.,3.,0.], [0.,2.,0.],[1.,2.,0.],[1.,1.,0.],[0.,1.,0.]
    ])

    # nodes = [LatticeNode(np.array(i)) for i in route_vertices]
    # route = LatticeRoute(nodes)

    """ Generate the route from polygon """
    route = generate_route(trapREV*2)
    route.plot()

    # """ Run side_staples """
    # collection = side_staples(route)
    # collection.plot(route)
    # system = route.system(collection.staples)
    # system.write_oxDNA("lol")
    

    """ Test Scaffold """
    # test, scaf = test_StapleRoute(route)
    # system = System(np.array([50,50,50]))
    # system.add_strands(scaf)
    # system.write_oxDNA("scaffold")
    scaffold = Scaffold(route)
    # scaffold.describe()