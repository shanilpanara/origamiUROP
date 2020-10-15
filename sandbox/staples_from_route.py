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

from origamiUROP.lattice.utils import find_crossover_locations

class StapleNode(DNANode):
    def __init__(self, position: np.ndarray):
        super().__init__(position)


class StapleEdge(DNAEdge):
    def __init__(self, vertex_1: StapleNode, vertex_2: StapleNode):
        super().__init__(vertex_1, vertex_2)


class StapleRoute(Strand):
    def __init__(self,
                 scaffold_rows: List[Strand],
                 nodes: List[StapleNode] = []):
        self._nodes = nodes
        self._scaffold_rows = scaffold_rows
        self.update_strand_and_nucleotides()
        super().__init__(self._nucleotides)

    @property
    def nodes(self) -> List[StapleNode]:
        return self._nodes

    @property
    def edges(self):
        _edges = [StapleEdge(node, self.nodes[i + 1]) for i, node in enumerate(self.nodes[:-1])]
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
            
            # Find start of staple
            node1 = int(edge.vertices[0][0])
            node2 = int(edge.vertices[1][0])

            x1 = node1
            x2 = node2
            # if len(self.scaffold_rows) - 1 == row:
            #     break
            
            row = int(edge.vertices[0][1]) # == int(edge.vertices[1][0])
            scaffold_row = self.scaffold_rows[row]
            nucleotides = []
            

            if x1 > x2:
                # this doesnt make sense to slice, hence switch ----x2-------x1  to  x1-------x2----
                x1, x2 = x2, x1
                for nucleotide in scaffold_row.nucleotides[x1:x2 + 1]:
                    nucleotides.append(nucleotide.make_across())
                # switch the order of the nucleotides back again
                nucleotides = nucleotides[::-1]
            else:
                for nucleotide in scaffold_row.nucleotides[::-1][x1:x2 + 1]:
                    nucleotides.append(nucleotide.make_across())

            for nucleotide in nucleotides:
                self._nucleotides.append(nucleotide)


class StapleCollection:
    def __init__(self, strands: List[Strand] = [], route: LatticeRoute = None):
        self._staples = strands
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
        return self._staples

    @property
    def n_staples(self) -> int:
        return len(self._staples)

    def add_staples(self, staple_strand: Strand):
        self._staples.append(staple_strand)

    def generate_staple(self, node1: int, node2: int, row_indices: List[int], shift) -> StapleRoute:
        node_cycle = cycle([node1, node2, node2, node1])
        staple_nodes = []

        if shift:
            next(node_cycle)
            next(node_cycle)

        for row in row_indices:
            # Generate nodes
            x_3p = next(node_cycle)
            staple_nodes.append(StapleNode([x_3p, row, 0.0]))
            x_5p = next(node_cycle)
            staple_nodes.append(StapleNode([x_5p, row, 0.0]))

            # Update stapled array
            bounds = [x_3p, x_5p]
            self.stapled_array[row][min(bounds):max(bounds) + 1] = 2

        return StapleRoute(self.scaffold_strands, staple_nodes)

    @staticmethod
    def compute_turning_points(start_site, start_side, staple_width, opposite = False):
        # staple must go in the OPPOSITE direction of strand
        node2 = start_site
        if not opposite:
            node1 = node2 + staple_width - 1 if start_side == "left" else node2 - staple_width + 1
        else:
            node1 = node2 + staple_width - 1 if start_side == "right" else node2 - staple_width + 1
        return node1, node2

    @staticmethod
    def all_same(items) -> bool:
        """Checks if all elements of a list are equal"""
        return all(x == items[0] for x in items)

    def run_algorithm_1(self, staple_width=16):
        # see pandas dataframe in self.scaffold_rows.info to understand values
        df = self.scaffold_rows.info()
        staples = []
        for i in range(len(df) - 1):
            if df.loc[i, "size"] < staple_width or df.loc[i+1, "size"] < staple_width:
                print(df.loc[i, "size"])
                staple_size = 5
            else:
                staple_size = staple_width

            """ COMPUTE TURNING POINTS """
            node1, node2 = self.compute_turning_points(df.loc[i, "start x"], df.loc[i, "start side"], staple_size)
            """ DETERMINE WHICH ROWS TO STAPLE OVER"""
            # if shift is True, staples are generated in the reverse direction in "generate_staple()"
            flip = False

            ### DEAL WITH 3 ROWS CASES ###
            # if i is not 0, the edge lines up on rows i-1 to i+1 and edge of i-1 != i-2 or i = 1 (where i-2 doesnt exist)
            if not i == 0 and len(set((df.loc[[i - 1, i + 1], "end x"].tolist() + [df.loc[i, "start x"]]))) == 1  and (i == 1 or df.loc[i - 1, "end x"] != df.loc[i - 2, "start x"]):
                rows = [i - 1, i, i + 1]
                flip = True
            # if i is not the last one, the edge lines up on rows i to i+1 and edge of i+1 != i+2 or i = N-2 (where i+3 doesnt exist)
            elif not i + 1 == len(df) - 1 and self.all_same(df.loc[[i, i + 2], "start x"].tolist() + [df.loc[i+1, "end x"]]) and (i == len(df) - 3 or df.loc[i + 2, "start x"] != df.loc[i + 3, "end x"]) :
                rows = [i, i + 1, i + 2]

            ### DEAL WITH CASES WHERE TOP IS SMALLER ###
            elif df.loc[i+1, "end x"] < df.loc[i, "start x"]:
                node1, node2 = self.compute_turning_points(df.loc[i, "end x"],df.loc[i, "start side"],
                                                           staple_size, opposite = True)
                rows = [i, i+1]
                flip = True
            
            ### DEAL WITH CASES WHERE TOP IS BIGGER ###
            elif df.loc[i+1, "end x"] > df.loc[i, "start x"]:
                node1, node2 = self.compute_turning_points(df.loc[i, "end x"],df.loc[i, "start side"],
                                                           staple_size, opposite = True)
                rows = [i, i+1]
                flip = True

            ### DEAL WITH THE CASES: ROW i == OCCUPIED ###
            ## SHORTEN
            # if True:
            #     pass
            # ## SKIP
            # else:
            #     continue

            else:
                rows = [i, i + 1]


            ### SHIFT NODE INDEX RELATIVE TO THE FIRST NUCLEOTIDE FROM THE STARTING SIDE
            # if df.loc[i, "start side"] == "right":
            #     node1 += df.loc[i, "r shift"]
            #     node2 += df.loc[i, "r shift"]
            # elif df.loc[i, "start side"] == "left":
            #     node1 += df.loc[i, "l shift"]
            #     node2 += df.loc[i, "l shift"]

            # GENERATE AND ADD STAPLES TO THE LIST
            staples.append(self.generate_staple(node1, node2, rows, flip))

        self._staples += staples

    def plot_nodes(self, strand: Strand, ax, colour='r', width=0.2, **kwargs):
        nodes = np.array(strand.nodes)
        #plt.grid(True)
        ax.plot(nodes[:, 0], nodes[:, 1], 'k', ms=0.5, alpha=0)
        ax.xaxis.set_major_locator(MultipleLocator(20))
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
                    color=colour,
                    head_width=0.1,
                    length_includes_head=True,
                    **kwargs)

    def plot(self, fout: str = None, colour: str = None):
        fig, ax = plt.subplots()
        if not colour:
            colours = cycle(('r', 'b', 'g', 'k'))
        else:
            colours = cycle((colour))
        if route:
            self.plot_nodes(strand = route, ax = ax, colour = 'k', width = 0.05, alpha = 0.4)
        for staple in self.staples:
            self.plot_nodes(strand = staple, ax = ax, colour = next(colours), alpha = 0.8)
        
        plt.gca().set_aspect(5)
        if fout:
            plt.savefig(fout, dpi=500)

        plt.show()


class Scaffold:
    def __init__(self, route: LatticeRoute, staple_widths=[5, 16, 26]):
        self.route = route
        self.edges = route.edges[0::2]  #only rows which are horizontal
        self.staple_widths = staple_widths
        self.sizes = self.get_row_sizes()
        self.start_side = self.get_start_side()

        self.startx = self.get_startx()
        self.endx = self.get_endx()
        self.bounds = self.get_bounds()

        self.array = self.get_scaffold_array()
        self.r_shift = self.get_r_shift()
        self.l_shift = self.get_l_shift()

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
        return [[int(edge.vertices[0][0]),
                 int(edge.vertices[1][0])] for edge in self.edges]

    def get_start_side(self) -> list:
        """ Returns left or right as the start side of each row"""
        nodes0 = self.route.nodes[0::2]
        nodes1 = self.route.nodes[1::2]
        return ["left" if node0[0] - node1[0] < 0 else "right" for (node0, node1) in zip(nodes0, nodes1)]

    def get_scaffold_array(self) -> np.ndarray:
        # Find max row width and no. of columns and create empty grid
        scaffold = np.zeros((len(self.edges), max(self.sizes)))
        for row in range(len(self.edges)):
            scaffold[row][min(self.bounds[row]):max(self.bounds[row]) + 1] = 1
        return scaffold
    
    def get_r_shift(self) -> int:
        max_x = max(self.sizes)
        return [max_x - start_x for start_x in self.startx]

    def get_l_shift(self) -> int:
        min_x = 0
        return [min_x - start_x for start_x in self.startx]

    def n_staples(self, staple_width=None):
        """ Returns no. of staples per row """
        if not staple_width:  # default to 16
            staple_width = self.staple_widths[1]
        return [int(i / staple_width) for i in self.sizes]

    def unpaired_bases(self, staple_width=None):
        """ Returns no. of unpaired bases per row """
        if not staple_width:  # default to 16
            staple_width = self.staple_widths[1]
        return [int(i % staple_width) for i in self.sizes]

    def info(self, staple_width=None):
        """
        Returns information for each row in the Lattice Route
            Size: no. of nt in row
            Bounds: [start, end] points of row
            Startx/Endx: same as above but split up
        And the number of staples 
        """
        return pd.DataFrame({
                "size": self.sizes,
                "start side": self.start_side,
                "start x": self.startx,
                "end x": self.endx,
                "bounds": self.bounds,
                "r shift": self.r_shift,
                "l shift": self.l_shift,
                "staples (5)": self.n_staples(self.staple_widths[0]),
                "staples (16)": self.n_staples(self.staple_widths[1]),
                "unpaired bases (5)": self.unpaired_bases(self.staple_widths[0]),
                "unpaired bases (16)": self.unpaired_bases(self.staple_widths[1]),
            })
        

def side_staples(route: LatticeRoute, staple_width=16):
    # import scaffold as formed of its constituent rows
    scaffold_rows = route.get_strand()[1]
    scaffold = Scaffold(route)
    staples = []

    for row, info in scaffold.info().iterrows():
        x2 = info["bounds"][0]
        if info["start side"] == "left":
            x1 = x2 + staple_width - 1
        elif info["start side"] == "right":
            x1 = x2 - staple_width + 1

        staple_nodes = [
            StapleNode([x1, row, 0.0]),
            StapleNode([x2, row, 0.0]),
            StapleNode([x2, row + 1, 0.0]),
            StapleNode([x1, row + 1, 0.0])
        ]

        staples.append(StapleRoute(scaffold_rows, staple_nodes))

    return staples


def test_StapleRoute(route: LatticeRoute):
    x1 = 0
    x2 = 15
    row = 0
    scaffold_rows = route.get_strand()[1]
    staple_nodes = [
        StapleNode([x1, row, 0.0]),
        StapleNode([x2, row, 0.0]),
        StapleNode([x2, row + 1, 0.0]),
        StapleNode([x1, row + 1, 0.0])
    ]

    return StapleRoute(scaffold_rows, staple_nodes)


def generate_route(polygon_vertices: np.ndarray):
    polygon = BoundaryPolygon(polygon_vertices)
    lattice = polygon.dna_snake(straightening_factor=5, start_side="left", grid_size= [0.34, 2.5])
    return lattice.route()


if __name__ == "__main__":
    half_turn_indices = [4, 15, 25, 36, 46, 56, 67, 77, 88, 98, 109]
    staple_lengths = [9, 31, 51, 73]

    route_vertices = [[0., 0., 0.], [56., 0., 0.], [56., 1., 0.], [0., 1., 0.],
                      [0., 2., 0.], [56., 2., 0.], [56., 3., 0.], [0., 3., 0.],
                      [0., 4., 0.], [56., 4., 0.], [56., 5., 0.], [0., 5., 0.],
                      [0., 6., 0.], [56., 6., 0.], [56., 7., 0.], [0., 7., 0.],
                      [0., 8., 0.], [88., 8., 0.], [88., 9., 0.], [0., 9., 0.],
                      [0., 10., 0.], [88., 10., 0.], [88., 11., 0.],
                      [0., 11., 0.], [0., 12., 0.], [56., 12., 0.],
                      [56., 13., 0.], [0., 13., 0.], [0., 14., 0.],
                      [56., 14., 0.], [56, 15, 0], [0, 15, 0]]
    trapREV = np.array([[0., 10., 0.], [2.5, 4., 0.], [7.5, 4., 0.],
                        [10., 10., 0.]])
    no = 0.85
    hexagon = np.array([[0, 1, 0], [no, 0.5, 0], [no, -0.5, 0], [0, -1, 0],
                        [-no, -0.5, 0], [-no, 0.5, 0]])
    plus = np.array([[1., 0., 0.], [2., 0., 0.], [2., 1., 0.], [3., 1., 0.],
                     [3., 2., 0.], [2., 2., 0.], [2., 3., 0.], [1., 3., 0.],
                     [1., 2., 0.], [0., 2., 0.], [0., 1., 0.], [1., 1., 0.]])
    diamond = np.array([[1., 0., 0.], [2., 1., 0.], [1., 2., 0.], [0., 1.,
                                                                   0.]])
    triangle = np.array([[0, 0, 0], [5, 9, 0], [10, 0, 0]])
    stacked_I = np.array([[0., 0., 0.], [3., 0., 0.], [3., 1.,
                                                       0.], [2., 1., 0.],
                          [2., 2., 0.], [3., 2., 0.], [3., 3.,
                                                       0.], [2., 3., 0.],
                          [2., 4., 0.], [3., 4., 0.], [3., 5.,
                                                       0.], [0., 5., 0.],
                          [0., 4., 0.], [1., 4., 0.], [1., 3.,
                                                       0.], [0., 3., 0.],
                          [0., 2., 0.], [1., 2., 0.], [1., 1., 0.],
                          [0., 1., 0.]])
    square = np.array([[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.]])
    # nodes = [LatticeNode(np.array(i)) for i in route_vertices]
    # route = LatticeRoute(nodes)
    """ Generate the route from polygon """
    route = generate_route(square * 20)
    route.plot()

    # """ Run side_staples """
    # collection = side_staples(route)
    # collection.plot(route)
    # system = route.system(collection)
    # system.write_oxDNA("trap w staples")
    """ Test Scaffold """
    # test, scaf = test_StapleRoute(route)
    # system = System(np.array([50,50,50]))
    # system.add_strands(scaf)
    # system.write_oxDNA("scaffold")
    # scaffold = Scaffold(route)
    # display(scaffold.info())
    # scaffold.describe()
    """ Test Staple Collection with new workflow """
    collection = StapleCollection(route=route)
    collection.run_algorithm_1()
    collection.plot()
    plt.imshow(collection.stapled_array[::-1], aspect = 5)
    system = route.system(collection.staples)
    system.write_oxDNA("stapled_hex")