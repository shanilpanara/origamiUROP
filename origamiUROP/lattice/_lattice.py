from typing import List

import numpy as np
from shapely.geometry import MultiPoint
from shapely import geometry

import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator

from operator import itemgetter
from copy import deepcopy
import itertools

from origamiUROP.lattice import edge, node, route
from origamiUROP.oxdna.strand import POS_STACK


def round_to_multiple(n, mo=0.34, decimal_places=2):
    """
    Function rounds to the nearest multiple of value given
    Returns output to (default) 2 decimal places

    Arguments:

    n --- value (integer or float) to round  
    mo --- "multiple_of" is the value (integer or float) which we want 
    to round to a multiple of  
    decimal_places --- no. of decimals to return
    """
    a = (n // mo) * mo  # Smaller multiple
    b = a + mo  # Larger multiple
    closest_multiple = b if n - a > b - n else a  # Return of closest of two
    return round(closest_multiple, decimal_places)


class Lattice:
    def __init__(
        self,
        polygon: np.ndarray,
        x_spacing: float = POS_STACK,
        y_spacing: float = 1.00,
    ):
        self.polygon = geometry.Polygon(polygon)
        self.x_spacing = x_spacing
        self.y_spacing = y_spacing
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
        x_min, y_min, x_max, y_max = polygon.bounds
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
        result = points.intersection(polygon)
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
    def intersect_array(self):
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

        # For the difference required to reach a multiple of 32 nt on each row
        # i.e. -5 or +11, remove (0) or add (1) to the lattice
        # alternating each removal or addition to either side
        for row, value in enumerate(nt_per_row_diff):
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

    @property
    def second_adjust_array(self):
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

        for row in range(np.shape(grid)[0]):
            if np.sum(grid[row]) == 0:  # to account for padding
                continue
            else:
                starting_location = int(np.argwhere(grid[row])[0])
                rounded_location = round_to_multiple(starting_location, 16, 0)
                change_in_location = rounded_location - starting_location
                grid[row] = np.roll(grid[row], change_in_location)
        return grid

    @property
    def adjusted_coords(self):
        """Convert binary array to a set of coordinates"""
        # Remove padding and find coordinates of updated system
        grid = deepcopy(self.second_adjust_array)
        y_min = np.argwhere(grid)[:, 0].min()
        y_max = np.argwhere(grid)[:, 0].max()
        x_min = np.argwhere(grid)[:, 1].min()
        x_max = np.argwhere(grid)[:, 1].max()
        print(
            "Converted binary to coords, the xy bounds are: \n,"
            f" x: {x_min}, {x_max}. y: {y_min}, {y_max}."
        )

        grid = grid[y_min : y_max + 1, x_min : x_max + 1]
        lattice = np.argwhere(grid)
        # swap columns, np.argwhere returns x and y the wrong way around
        lattice[:, 0], lattice[:, 1] = lattice[:, 1], lattice[:, 0].copy()
        return lattice

    def plotArray(
        self,
        value: np.ndarray,
        ax: plt.Axes = None,
        fout: str = None,
        show: bool = False,
    ):

        nodes = np.array(value)
        if not ax:
            fig, ax = plt.subplots()
        plt.grid(True)
        plt.rc("xtick", labelsize=6)  # fontsize of the tick labels
        plt.rc("ytick", labelsize=6)
        ax.xaxis.set_major_locator(MultipleLocator(32))
        ax.yaxis.set_major_locator(MultipleLocator(10))
        ax.set_xlabel("No. of nucleotides")
        ax.set_ylabel("No. of strands")
        # ax.plot(nodes[:, 0], nodes[:, 1], "ko", ms=0.5)
        ax.imshow(nodes[::-1])
        plt.gca().set_aspect("equal")
        if fout:
            plt.savefig(f"{fout}.png", dpi=500)
        if show:
            plt.show()

    def plotCoords(
        self,
        value: np.ndarray,
        ax: plt.Axes = None,
        fout: str = None,
        show: bool = False,
    ):

        nodes = np.array(value)
        if not ax:
            fig, ax = plt.subplots()
        plt.grid(True)
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontsize(4)
        ax.xaxis.set_major_locator(MultipleLocator(32))
        ax.yaxis.set_major_locator(MultipleLocator(10))
        ax.set_xlabel("No. of nucleotides")
        ax.set_ylabel("No. of strands")
        ax.plot(nodes[:, 0], nodes[:, 1], "ko", ms=0.5)
        plt.gca().set_aspect("equal")
        if fout:
            plt.savefig(f"{fout}.png", dpi=500)
        if show:
            plt.show()


square = np.array([[0, 0, 0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
trapezium = np.array([[0, 0, 0], [1.0, 10.0, 0.0], [2.0, 10.0, 0.0], [3.0, 0.0, 0.0]])
line = np.array([[0, 0, 0], [40, 0, 0], [40, 20, 0], [0, 20, 0]])

no = 0.866025
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
    polygon = geometry.Polygon(trapezium * 30)
    lattice = Lattice(polygon)
    # lattice.plotArray(lattice.adjusted_array, fout="First")
    # lattice.plotArray(lattice.second_adjust_array, fout="Second")
    lattice.plotCoords(lattice.adjusted_coords, fout="Straightened")
