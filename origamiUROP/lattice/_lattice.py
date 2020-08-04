from typing import List

import numpy as np
from shapely.geometry import MultiPoint
from shapely import geometry

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.pyplot import MultipleLocator

from operator import itemgetter
from copy import deepcopy
import itertools

from origamiUROP.lattice import edge, node, route
from origamiUROP.oxdna.strand import POS_STACK


def find_crossover_locations(
    max_size: int = 30, bp_per_turn: float = 10.45, kind: str = "half",
) -> dict:
    """
    Function which returns an array containing the best locations 
    for half or whole turns of the DNA strand where the strand will
    stay in plane and have the least strain/stress build-up.
    
    This is computed given the number of basepairs per turn
    & locations are given in no. of base pairs

    Arguments:

        bp_per_turn --- (default: 10.45)
        kind --- either half turns or whole turns
        max_size --- max no. of turns to base pairs to the list
    """
    crossover_location = [0]
    no_of_bp = 0

    # half turn
    if kind == "half":
        i = 1  # every odd number
        while no_of_bp < max_size:
            # minus 1 at the end because Python indexing starts at 0
            no_of_bp = int(round(bp_per_turn * 0.5 * i, 0) - 1)
            crossover_location.append(no_of_bp)
            i += 2

    # whole turns
    elif kind == "whole":
        i = 2  # every even number
        while no_of_bp < max_size:
            # minus 1 at the end because Python indexing starts at 0
            no_of_bp = int(round(bp_per_turn * 0.5 * i, 0) - 1)
            crossover_location.append(no_of_bp)
            i += 2

    return crossover_location


def round_to_multiple(n, mo=0.34, decimal_places=2):
    """
    Function rounds to the nearest multiple of value given
    Returns output to (default) 2 decimal places

    Arguments:

        n --- value (integer or float) to round  
        mo --- "multiple_of" is the value (integer or float) which 
        we want to round to a multiple of  
        decimal_places --- no. of decimals to return
    """
    a = (n // mo) * mo  # Smaller multiple
    b = a + mo  # Larger multiple
    closest_multiple = b if n - a > b - n else a  # Return of closest of two
    return round(closest_multiple, decimal_places)


def modify_lattice_row(grid: np.ndarray, difference: np.ndarray):
    """
    Modifies the lattice sites for a given row in the `grid` by 
    adding or removing the number of sites equal to the value which 
    is stored in `difference` for that row

    Arguments:
        grid --- (n, x) numpy array
        difference --- (n, 1) numpy array
    
    e.g. a difference of -5 or +11, remove (1 -> 0) or add (0 -> 1) 
    to the lattice row, respectively. Each removal/addition occurs on 
    the opposite side of the lattice than the previous removal/addition
    

    """
    # For the difference required to reach a certain no. of nucleotides on each row

    for row, value in enumerate(difference):
        # Alternate between left and right
        sides = {"left": "right", "right": "left"}
        side = "left"  # begin at the left side
        while value != 0:
            # find location first/last 1 in the grid
            # i.e. one_right = x-coord (x location of the last 1 from the right)
            one_right = int(np.argwhere(grid[row])[:, 0].max())
            one_left = int(np.argwhere(grid[row])[:, 0].min())

            if np.sum(grid[row]) == 0:
                value = 0

            if value > 0:  #  add 1's to edges
                if side == "left":
                    grid[row, one_left - 1] = 1
                else:
                    grid[row, one_right + 1] = 1
                value -= 1
            else:  # value < 0 i.e. remove 1's from edges
                if side == "left":
                    grid[row, one_left] = 0
                else:
                    grid[row, one_right] = 0
                value += 1
            side = sides[side]  # change left -> right OR right -> left

    return grid


class Lattice:
    def __init__(
        self,
        polygon: np.ndarray,
        x_spacing: float = POS_STACK,
        y_spacing: float = 1.00,
        bp_per_turn: float = 10.45,
    ):
        self.polygon_array = polygon.astype(dtype=np.float64)
        self.polygon = geometry.Polygon(polygon)
        self.x_spacing = x_spacing
        self.y_spacing = y_spacing
        self.bp_per_turn = bp_per_turn
        """
        Lattice class forms a set of points given a polygon where a DNA Origami Scaffold
        can be laid.

        `lattice` refers to a coordinate system
        `grid` refers to an array of 1's and 0's

        Arguments:
        polygon --- a set of vertices given in order, which form a closed shape
        x_spacing --- (default = 0.34) the x distance between each site in the coordinate grid
        y_spacing --- (default = 1.00) the y distance between each site in the coordinate grid
        """

    @property
    def intersect_coords(self):
        # Calculate Bounds of the box around lattice
        x_min, y_min, x_max, y_max = self.polygon.bounds
        x_min = round_to_multiple(x_min, self.x_spacing)
        x_max = round_to_multiple(x_max, self.x_spacing)
        y_min = round_to_multiple(y_min, self.y_spacing)
        y_max = round_to_multiple(y_max, self.y_spacing)

        # Initialise large grid of points spaced using x & y spacing
        # This is to "quantise" our grid
        x_grid = np.linspace(x_min, x_max, int((x_max - x_min) / self.x_spacing) + 2)
        y_grid = np.linspace(y_min, y_max, int((y_max - y_min) / self.y_spacing) + 1)
        bounds_grid = np.transpose(
            [np.tile(x_grid, len(y_grid)), np.repeat(y_grid, len(x_grid))]
        )
        points = MultiPoint(bounds_grid)

        # Shapely intersect grid with polygon and store lattice coordinates in "lattice"
        # (Shapely isn't well documented enough and there might have been a better way
        # to do this)
        result = points.intersection(self.polygon)
        lattice = np.array(result.__geo_interface__["coordinates"])

        # Normalise values to integer values, e.g. 0.34 -> 1
        lattice[:, 0] /= self.x_spacing
        lattice[:, 1] /= self.y_spacing
        lattice = np.around(lattice).astype(int)

        # Move grid to begin at [0,0,(0)]
        x_min = lattice[:, 0].min()
        y_min = lattice[:, 1].min()
        lattice[:, 0] -= x_min
        lattice[:, 1] -= y_min
        return lattice

    @property
    def intersect_array(self) -> np.ndarray:
        """ 
        Make binary array representing the lattice
        where 1 = scaffold site, 0 = not a scaffold site
        """
        lattice = deepcopy(self.intersect_coords)
        y_max = lattice[:, 1].max()
        x_max = lattice[:, 0].max()

        grid = np.zeros((y_max + 1, x_max + 1)).astype(int)
        for point in lattice:
            grid[point[1], point[0]] = 1
        return grid

    @property
    def adjusted_array(self):
        """
        Adjust each row to contain a multiple of 32 scaffold sites
        Values are rounded up OR down to the closest multiple
        """
        # Add a border of 16 0's around the lattice
        grid = deepcopy(self.intersect_array)
        grid = np.pad(grid, pad_width=16, mode="constant", constant_values=0)

        # Find the number of nucleotide sites per row and store in a numpy array
        nt_per_row = []
        for row in range(len(grid)):
            nt_per_row.append(np.sum(grid[row]))
        nt_per_row = np.array(nt_per_row)

        rounded_to_32 = lambda x: round_to_multiple(x, 32, 0)
        nt_per_row_round = np.array(list(map(rounded_to_32, nt_per_row)))
        nt_per_row_diff = nt_per_row_round - nt_per_row

        grid = modify_lattice_row(grid, nt_per_row_diff)

        return grid

    @property
    def straight_edge_array(self):
        """
        Algorithm to straighten out edges of the lattice
        - achieved by shifting each row to begin at a multiple of 16
        """
        grid = deepcopy(self.adjusted_array)

        # For a given number of rows
        #   Find starting_location
        #   Calc Rounded_location
        #   Change_location = Rounded - starting
        #   np.roll(grid[row])
        straightening_factor = 4
        for row in range(np.shape(grid)[0]):
            if np.sum(grid[row]) == 0:  # to account for padding
                continue
            else:
                starting_location = int(np.argwhere(grid[row])[0])
                rounded_location = round_to_multiple(starting_location, 4, 0)
                change_in_location = rounded_location - starting_location
                grid[row] = np.roll(grid[row], change_in_location)
        return grid

    @property
    def final_coords(self):
        """Convert binary array to a set of coordinates"""
        # Remove padding and find coordinates of updated system
        grid = deepcopy(self.straight_edge_array)
        y_min = np.argwhere(grid)[:, 0].min()
        y_max = np.argwhere(grid)[:, 0].max()
        x_min = np.argwhere(grid)[:, 1].min()
        x_max = np.argwhere(grid)[:, 1].max()
        print(
            "Converted binary to coords, the xy bounds are: \n"
            f"x: {x_min}, {x_max}. y: {y_min}, {y_max}."
        )

        grid = grid[y_min : y_max + 1, x_min : x_max + 1]
        lattice = np.argwhere(grid)
        # swap columns, np.argwhere returns x and y the wrong way around
        lattice[:, 0], lattice[:, 1] = lattice[:, 1], lattice[:, 0].copy()
        return lattice

    def crossovers_binary(self):
        grid = deepcopy(self.straight_edge_array)
        # copy the size of grid but fill it with zeros
        crossovers_array = np.zeros(np.shape(grid))

        max_width = np.shape(grid)[1]
        possible_crossovers = find_crossover_locations(max_width)
        print("Possible Crossovers", possible_crossovers)
        for row in range(np.shape(grid)[0]):

            # to account for padding
            lattice_sites = np.sum(grid[row])
            if lattice_sites == 0:
                continue

            # find x of first point
            left_bound = np.argwhere(grid[row])[0]
            for i, bp in enumerate(possible_crossovers):
                if bp > lattice_sites:
                    continue
                crossovers_array[row, left_bound + bp] = 1 * grid[row, left_bound + bp]

        return crossovers_array

    def crossovers_coords(self):
        """Convert binary array to a set of coordinates"""
        crossovers = deepcopy(self.crossovers_binary())
        y_min = np.argwhere(crossovers)[:, 0].min()
        y_max = np.argwhere(crossovers)[:, 0].max()
        x_min = np.argwhere(crossovers)[:, 1].min()
        x_max = np.argwhere(crossovers)[:, 1].max()
        print(
            "Converted crossover binary to coords, the xy bounds are: \n"
            f"x: {x_min}, {x_max}. y: {y_min}, {y_max}."
        )

        crossovers = crossovers[y_min : y_max + 1, x_min : x_max + 1]
        coords = np.argwhere(crossovers)
        # swap columns, np.argwhere returns x and y the wrong way around
        coords[:, 0], coords[:, 1] = coords[:, 1], coords[:, 0].copy()

        return coords

    def plotPolygon(self, ax, nodes: np.ndarray, coords: bool):
        polygon = deepcopy(self.polygon_array)
        x_min = polygon[:, 0].min()
        y_min = polygon[:, 1].min()

        # shift polygon to (0,0)
        polygon[:, 0] -= x_min
        polygon[:, 1] -= y_min

        # Normalise values to integer values, e.g. 0.34 -> 1
        polygon[:, 0] /= self.x_spacing
        polygon[:, 1] /= self.y_spacing

        # add first value to start
        polygon = np.vstack((polygon, polygon[0]))

        if coords:
            # Center polygon
            polygon[:, 0] += (nodes[:, 0].max() - polygon[:, 0].max()) / 2
        else:
            # flip
            polygon[:, 1] = polygon[:, 1].max() - polygon[:, 1]
            # add padding to match
            polygon += [16, 16, 0]

        ax.plot(polygon[:, 0], polygon[:, 1], "r-")

        return ax

    def plotCrossovers(self, ax, coords: bool):

        if coords:
            crossovers = deepcopy(self.crossovers_coords())
            ax.plot(
                crossovers[:, 0], crossovers[:, 1], "bo", ms=1.5,
            )
        else:
            crossovers = deepcopy(self.crossovers_binary())
            crossovers = np.ma.masked_where(crossovers == 0, crossovers)
            cmap = plt.get_cmap("cool")
            ax.imshow(crossovers[::-1], cmap)

    def plot(
        self,
        lattice: np.ndarray,
        ax: plt.Axes = None,
        fout: str = None,
        show: bool = False,
    ):

        nodes = np.array(lattice)
        if not ax:
            fig, ax = plt.subplots()
        # plt.grid(True)
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontsize(4)
        ax.xaxis.set_major_locator(MultipleLocator(32))
        ax.yaxis.set_major_locator(MultipleLocator(10))
        ax.set_xlabel("No. of nucleotides")
        ax.set_ylabel("No. of strands")

        if np.shape(nodes)[1] == 2:  # Plot Coordinates
            ax.plot(nodes[:, 0], nodes[:, 1], "ko", ms=0.5)
            self.plotPolygon(ax, nodes, coords=True)
            self.plotCrossovers(ax, coords=True)
            ax.set_title("Lattice plotted from coords")

        else:  # Plot binary array
            cmap = plt.get_cmap("hot")
            ax.imshow(nodes[::-1], cmap)
            self.plotPolygon(ax, nodes, coords=False)
            self.plotCrossovers(ax, coords=False)
            ax.set_title("Lattice plotted from array")

        plt.gca().set_aspect(3)
        plt.show()
        if fout:
            plt.savefig(f"{fout}.png", dpi=500)


square = np.array([[0, 0, 0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
trapezium = np.array([[0, 0, 0], [1.0, 10.0, 0.0], [12.0, 10.0, 0.0], [13.0, 0.0, 0.0]])
line = np.array([[0, 0, 0], [40, 0, 0], [40, 20, 0], [0, 20, 0]])

no = 0.5
hexagon = np.array(
    [[0, 1, 0], [no, 0.5, 0], [no, -0.5, 0], [0, -1, 0], [-no, -0.5, 0], [-no, 0.5, 0]]
)

star = np.array(
    [
        [0, 300, 0],
        [100, 110, 0],
        [300, 70, 0],
        [160, -90, 0],
        [190, -300, 0],
        [0, -210, 0],
        [-190, -300, 0],
        [-160, -90, 0],
        [-300, 70, 0],
        [-100, 110, 0],
    ]
)

if __name__ == "__main__":
    polygon = trapezium * 5
    lattice = Lattice(polygon)
    # lattice.plot(lattice.adjusted_array, fout="Unstraightened")
    # lattice.plot(lattice.adjusted_coords, fout="polygonCoords")
    # lattice.plot(lattice.final_coords)
    lattice.plot(lattice.final_coords, fout="Crossovers_on_Trapezium_Coords")

    print("completed")
