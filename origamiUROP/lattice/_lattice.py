from typing import List

import numpy as np
from shapely.geometry import MultiPoint
from shapely import geometry

import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator

# from operator import itemgetter
from copy import deepcopy
import itertools

from .node import LatticeNode
from .route import LatticeRoute
from ..oxdna.strand import POS_STACK, FENE_LENGTH

from bisect import bisect_left, bisect_right


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


def find_closest(myList, myNumber):
    """
    Credit: https://stackoverflow.com/questions/12141150
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return after + 1
    else:
        return before + 1


def modify_lattice_row(grid: np.ndarray, difference: np.ndarray or List, change_side: str = "yes"):
    """
    Modifies the lattice sites for a given row in the `grid` by 
    adding or removing the number of sites equal to the value which 
    is stored in `difference` for that row

    Arguments:
        grid --- (n, x) numpy array of 1's and 0's
        difference --- (n, 1) numpy array
        change_side --- ("Yes") otherwise "onlyleft" or "onlyright"
        arguments can be given, where adjustments will only be made
        on the prescribed side
    
    e.g. a difference of -5 or +11, remove (1 -> 0) or add (0 -> 1) 
    to the lattice row, respectively. Each removal/addition occurs on 
    the opposite side of the lattice than the previous removal/addition
    

    """
    # Alternate between left and right
    sides = {"left": "right", "right": "left","onlyleft":"left", "onlyright":"right"}
    change_side = change_side.lower()
    message = "change_side can only have the following values: \n'yes', \n'onlyleft', \n'onlyright'"
    assert change_side in ("yes", "onlyleft", "onlyright"), message

    if not isinstance(difference, np.ndarray):
        try:
            difference = np.array([difference])
        except TypeError:
            raise TypeError(f"Difference needs to be given as a list, not {type(difference)}")
    
    if len(np.shape(grid)) != 2:
        try:
            grid = np.array([grid])
            assert len(np.shape(grid)) == 2
        except TypeError:
            raise TypeError("The value hasn't been given as a list or a np.ndarray")

    # For the difference required to reach a certain no. of nucleotides on each row
    for row, value in enumerate(difference):
        if change_side == "yes":
            side = "left"  # begin at the left side
        else:
            side = sides[change_side]

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

            if change_side == "yes":
                side = sides[side]  # change left -> right OR right -> left

    return grid

def side(row_no: int, R: bool = None):
    """
    Returns either 'onlyleft' or 'onlyright' depending on 
    which row we are modifying, i.e. row0 will always be modified
    at its end so if R = True, then right is the end of row0,
    row1 will always be modified at the opposite side to row0
    """
    dict = {(0,None):"onlyleft", 
            (0, True):"onlyright", 
            (1,None):"onlyright",
            (1,True):"onlyleft"}
    return dict[(row_no, R)]

def all_same(items):
    return all(x == items[0] for x in items)

class Lattice:
    def __init__(
        self,
        polygon_vertices: np.ndarray,
        grid_size: List[float] = [POS_STACK, 1.00],
        bp_per_turn: float = 10.45,
        straightening_factor: int = 8,
        start_side = "left"
    ):
        """
        Lattice class forms a set of points given a polygon where a DNA Origami Scaffold
        can be laid.

        `lattice` refers to a coordinate system & `grid` refers to an array of 1's and 0's
        Both `lattice` and `grid` represent sites where scaffold can be laid

        Arguments:
        polygon --- a set of vertices given in order, which form a closed shape
        x_spacing --- (default = 0.34) the x distance between each site in the coordinate grid
        y_spacing --- (default = 1.00) the y distance between each site in the coordinate grid
        """
        self.polygon_array = polygon_vertices.astype(dtype=np.float64)
        self.polygon = geometry.Polygon(polygon_vertices)
        self.grid_size = grid_size
        self.x_spacing = grid_size[0]
        self.y_spacing = grid_size[1]
        self.bp_per_turn = bp_per_turn
        self.straightening_factor = straightening_factor
        self.padding = 40
        self.start_side = start_side.lower()
        
        assert self.start_side in ["left","right"], "start_side must be: 'left' or 'right'"

    @property
    def intersect_polygon(self):  # Returns coordinates
        """
        Creates grid with dimensions contained in self.gridsize;
        Creates an intersection with the polygon;
        
        Returned value: set of coordinates detailing only the grid
        points which were overlaid on the polygon
        """
        # Calculate bounds of the grid around the polygon
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
    def intersect_array(self):  # Returns array
        """ 
        Make binary array representing the lattice sites
        where 1 = scaffold site, 0 = not a scaffold site
        """
        lattice = deepcopy(self.intersect_polygon)
        y_max = lattice[:, 1].max()
        x_max = lattice[:, 0].max()

        grid = np.zeros((y_max + 1, x_max + 1)).astype(int)
        for point in lattice:
            grid[point[1], point[0]] = 1
        return grid

    @property
    def quantise_array_rows(self):  # Returns array
        """
        Adjust each row to contain scaffold sites where the first and last
        site correlate to a position where the crossover occurs.
        Values are rounded up OR down to the closest half turn location.

        i.e. a row of length 9 will round down to 5 lattice sites
        or 43 will round up to 46
        """
        grid = deepcopy(self.intersect_array)

        # Find possible crossover locations
        max_width = np.shape(grid)[1]
        self.poss_cross = find_crossover_locations(max_width+self.padding)

        # Add a border of 16 0's around the lattice
        grid = np.pad(grid, pad_width=self.padding, mode="constant", constant_values=0)

        # Find the number of nucleotide sites per row and store in a numpy array
        nt_per_row = []
        for row in range(len(grid)):
            nt_per_row.append(np.sum(grid[row]))
        nt_per_row = np.array(nt_per_row)

        # Find closest crossover and modify rows ensure first/last sites are crossover sites

        # values 1-4 round up
        nt_per_row_round_1 = [5 if 0 < i < 5 else i for i in nt_per_row]
        # other values round to closest half turn 
        closest_crossover = lambda x: find_closest(self.poss_cross, x)
        nt_per_row_round_2 = np.array(list(map(closest_crossover, nt_per_row_round_1)))
        nt_per_row_diff = nt_per_row_round_2 - nt_per_row

        grid = modify_lattice_row(grid, nt_per_row_diff)

        return grid

    @property
    def straight_edge_array(self):
        """
        Algorithm to straighten out edges of the lattice
        - achieved by shifting each row to begin at a multiple of "sf" straightening_factor
        """
        grid = deepcopy(self.quantise_array_rows)

        # For a given number of rows
        #   Find starting_location
        #   Calc Rounded_location
        #   Change_location = Rounded - starting
        #   np.roll(grid[row])

        sf = self.straightening_factor
        for row in range(np.shape(grid)[0]):
            if np.sum(grid[row]) == 0:  # to account for padding
                continue
            else:
                starting_location = int(np.argwhere(grid[row])[0])
                rounded_location = round_to_multiple(starting_location, sf, 0)
                change_in_location = rounded_location - starting_location
                grid[row] = np.roll(grid[row], change_in_location)
        return grid

    @property
    def align_rows_array(self):
        """
        Algorithm to ensure every row connects to the next
        """            

        grid = deepcopy(self.straight_edge_array)
        # half turn nt sizes [5,16,26,37,47, etc] instead of half turn indexes [0,4,15, etc]
        poss_cross = np.add(self.poss_cross[1:],1)
        width_tracker = []

        for index in range(np.shape(grid)[0]-3):

            row0 = grid[index] 
            row1 = grid[index+1]
            row2 = grid[index+2]
            row3 = grid[index+3]
            row_width = [np.sum(row0),np.sum(row1), np.sum(row2), np.sum(row3)]
            width_tracker.append(int(row_width[0]))
            if 0 in row_width[0:2]: continue # just to bext index loop
            
            # R represents the side where the first crossover occurs
            # Hence, if self.start_side = "left"
            # for 0,2,4 etc -> R = True, otherwise R = None
            
            if self.start_side == "left":
                R = True if index % 2 == 0 else None
            else:
                R = None if index % 2 == 0 else True

            # calculate current difference between end/start points of row0 & row1
            if R:
                right0 = int(np.argwhere(row0)[-1])
                right1 = int(np.argwhere(row1)[-1])
                if right0-right1 == 0: continue
                assert right0-right1 != 0
            else:
                left0 = int(np.argwhere(row0)[0])
                left1 = int(np.argwhere(row1)[0])
                if left0-left1 == 0: continue
                assert left0 - left1 != 0
            
            ### Initialise some useful values
            # find index in poss_cross correlating to no. of possible crossovers in that row
            cross_idx_bottom = bisect_left(poss_cross, row_width[0])
            cross_idx_top = bisect_left(poss_cross, row_width[1])
            TminusB = row_width[1]-row_width[0]
            top_bigger = True if TminusB > 0 else False
            extra_turns = cross_idx_bottom - cross_idx_top # in whole row
            row0_diff = 0

            # Same size rows, just shift them left or right
            if TminusB == 0:
                if R:
                    row1_roll = right0 - right1
                else:
                    row1_roll = left0 - left1    
                row1 = np.roll(row1, row1_roll)          
            
            elif abs(extra_turns) > 1:         
                if R:
                    extra_turns_right = int(round(abs(right1-right0) / self.bp_per_turn, 0))
                    if not top_bigger:
                        # shorten the end of row 0
                        if row_width[0] <= poss_cross[3] and not all_same(width_tracker[-3:-1]):
                            row0_diff = poss_cross[cross_idx_bottom - 1] - row_width[0]
                        else:
                            row0_diff = poss_cross[cross_idx_bottom - extra_turns_right] - row_width[0]
                    else:
                        # lengthen the end of row 0
                        if row_width[0] <= poss_cross[3] and not all_same(width_tracker[-3:-1]):
                            row0_diff = poss_cross[cross_idx_bottom + 1] - row_width[0]
                        else:
                            row0_diff = poss_cross[cross_idx_bottom + extra_turns_right] - row_width[0] 
                    

                elif not R:
                    extra_turns_left = int(round(abs(left1-left0) / self.bp_per_turn, 0))
                    if not top_bigger:
                        # shorten the end of row 0
                        if row_width[0] <= poss_cross[3] and not all_same(width_tracker[-3:-1]):
                            row0_diff = poss_cross[cross_idx_bottom - 1] - row_width[0]
                        else:
                            row0_diff = poss_cross[cross_idx_bottom - extra_turns_left] - row_width[0] 
                    else:
                        # lengthen the end of row 0
                        if row_width[0] <= poss_cross[3] and not all_same(width_tracker[-3:-1]):
                            row0_diff = poss_cross[cross_idx_bottom + 1] - row_width[0]
                        else:
                            row0_diff = poss_cross[cross_idx_bottom + extra_turns_left] - row_width[0] 
                
                row0 = modify_lattice_row(row0, row0_diff, side(0, R))

                if R:
                    row1_roll = right0 - right1 + row0_diff
                    row1 = np.roll(row1, row1_roll)
                elif not R:
                    row1_roll = left0 - left1 - row0_diff
                    row1 = np.roll(row1, row1_roll)

            elif abs(extra_turns) == 1:
                if cross_idx_bottom >= 1: #16 or bigger            
                    if top_bigger:
                        if not row_width[1] >= max(width_tracker[-4:-1]) and row_width[1] == row_width[2] == row_width[3]:
                            row1_diff = poss_cross[cross_idx_top - 1] - row_width[1]
                            row1 = modify_lattice_row(row1, row1_diff, side(1,R)) 
                        # remove 1 to row1 -> so bigger shapes are less skewed
                        # only if row width ISNT the bigger and the next few rows arent equal to each other
                        elif cross_idx_bottom > 4 and not row_width[1] > max(width_tracker[-3:-1]): 
                            row1_diff = poss_cross[cross_idx_top - 1] - row_width[1]
                            row1 = modify_lattice_row(row1, row1_diff, side(1,R)) 
                            # pass
                        # if the next 2 rows are equal and the past few rows were bigger than row0
                        elif row_width[2] in [row_width[1], 0]:
                            if max(width_tracker[-3:-1]) >= row_width[0]:
                                # shorten row1 by 1 crossover
                                row1_diff = poss_cross[cross_idx_top - 1] - row_width[1]
                                row1 = modify_lattice_row(row1, row1_diff, side(1,R))
                    
                    elif not top_bigger:
                        if width_tracker[-2] == row_width[1] == row_width[2] and not row_width[0] > max(width_tracker[-6:-1]):
                            # shorten ends of row 0/1 by 1 crossover
                            print(index - 40, row_width[0:3])
                            row0_diff = poss_cross[cross_idx_bottom - 1] - row_width[0]
                            row0 = modify_lattice_row(row0, row0_diff, side(0,R)) 
                        elif min(width_tracker[-3:-1]) in [max(width_tracker[-3:-1]),0]:
                            if row_width[3] == row_width[2] == row_width[1] or row_width[3] == row_width[2] == 0:  
                                if cross_idx_bottom > 4:
                                    # shorten ends of row 0/1 by 1 crossover
                                    row0_diff = poss_cross[cross_idx_bottom - 1] - row_width[0]
                                    row0 = modify_lattice_row(row0, row0_diff, side(0,R)) 
                                    row1_diff = poss_cross[cross_idx_top - 1] - row_width[1]
                                    row1 = modify_lattice_row(row1, row1_diff, side(1,R))
                        elif cross_idx_bottom == 3 and row_width[1]==row_width[2]:
                            row1_diff = poss_cross[cross_idx_top + 1 ] - row_width[1]
                            row1 = modify_lattice_row(row1, row1_diff, side(1,R)) 
                        # where the last 2 rows have been equal and the next 3 rows are equal (or equal to 0)
                        # add 1 from row1 -> so bigger shapes are less skewed
                        elif cross_idx_bottom > 4 and not row_width[1] > max(width_tracker[-3:-1]): 
                            row1_diff = poss_cross[cross_idx_top + 1 ] - row_width[1]
                            row1 = modify_lattice_row(row1, row1_diff, side(1,R)) 
                            # pass



                elif cross_idx_bottom == 0: #5
                    if top_bigger:
                        if row_width[2] < row_width[1]:
                            # row 1: 16 --> 5
                            row1_diff = poss_cross[cross_idx_top - 1] - row_width[1]
                            row1 = modify_lattice_row(row1, row1_diff, side(1,R))
                    elif not top_bigger:
                        #just shift
                        pass


                if R:
                    row1_roll = right0 - right1 + row0_diff
                    row1 = np.roll(row1, row1_roll)
                elif not R:
                    row1_roll = left0 - left1 - row0_diff
                    row1 = np.roll(row1, row1_roll)
            
            # UPDATE ROWS IN GRID
            grid[index] = row0
            grid[index+1] = row1
            width_tracker[index] = np.sum(row0)

        return grid

    @property
    def final_array(self):
        """
        Returns lattice points as an array of 1's and 0's

        1's represent lattice sites where scaffold can be placed
        """
        return self.align_rows_array

    @property
    def final_coords(self):
        """Returns lattice points as a set of coordinates"""
        # Remove padding and find coordinates of updated system
        grid = deepcopy(self.straight_edge_array)
        y_min = np.argwhere(grid)[:, 0].min()
        y_max = np.argwhere(grid)[:, 0].max()
        x_min = np.argwhere(grid)[:, 1].min()
        x_max = np.argwhere(grid)[:, 1].max()

        grid = grid[y_min : y_max + 1, x_min : x_max + 1]
        lattice = np.argwhere(grid)
        # swap columns, np.argwhere returns x and y the wrong way around
        lattice[:, 0], lattice[:, 1] = lattice[:, 1], lattice[:, 0].copy()
        return lattice

    def get_crossovers(self):
        """ Returns a binary array, where 1 represents a potential crossover location"""
        grid = deepcopy(self.final_array)
        # copy the size of grid but fill it with zeros
        crossovers_array = np.zeros(np.shape(grid))

        sides = {"left": "right", "right": "left"}
        side = self.start_side  # begin at the left / rightside
        for row in range(np.shape(grid)[0]):
            # to account for padding, calc no. of lattice sites on row
            lattice_sites = np.sum(grid[row])
            if lattice_sites == 0:  # i.e. when row is just padding
                continue

            # find x of first/last lattice site in the row
            left_bound = np.argwhere(grid[row])[0]
            right_bound = np.argwhere(grid[row])[-1]

            for bp in self.poss_cross:
                if bp > lattice_sites:
                    continue
                if side == "left":
                    crossovers_array[row, left_bound + bp] = (
                        1 * grid[row, left_bound + bp]
                    )
                else:  # if side == right
                    crossovers_array[row, right_bound - bp] = (
                        1 * grid[row, right_bound - bp]
                    )

            side = sides[side]

        return crossovers_array

    def get_crossovers_coords(self):
        """Convert binary array to a set of coordinates"""
        crossovers = deepcopy(self.get_crossovers())
        y_min = np.argwhere(crossovers)[:, 0].min()
        y_max = np.argwhere(crossovers)[:, 0].max()
        x_min = np.argwhere(crossovers)[:, 1].min()
        x_max = np.argwhere(crossovers)[:, 1].max()

        crossovers = crossovers[y_min : y_max + 1, x_min : x_max + 1]
        coords = np.argwhere(crossovers)
        # swap columns, np.argwhere returns x and y the wrong way around
        coords[:, 0], coords[:, 1] = coords[:, 1], coords[:, 0].copy()

        return coords

    def route(self, *args, **kwargs) -> List[LatticeNode]:
        """Generate Scaffold Route, returns list of LatticeNode objects"""
        coords = deepcopy(self.get_crossovers_coords())

        # add 3rd column (z = 0)
        shape = np.shape(coords)
        if shape[1] != 3:
            coords = np.c_[coords, np.zeros(shape[0])]

        crossovers_per_row = 2
        lattice_rows = int(coords[:, 1].max() + 1)  # coords counts from 0
        vertices_in_route = int(crossovers_per_row * lattice_rows)
        vertex_list = np.zeros((vertices_in_route, 3))

        # Find final crossover from left or right and make it a node
        for row in range(0, lattice_rows):
            vertex_index_L = bisect_left(coords[:, 1], row)
            vertex_index_R = bisect_right(coords[:, 1], row) - 1
            if row % 2 == 0:  # if even
                vertex_list[row * 2] = coords[vertex_index_L]
                vertex_list[row * 2 + 1] = coords[vertex_index_R]
            else:  # if odd
                vertex_list[row * 2] = coords[vertex_index_R]
                vertex_list[row * 2 + 1] = coords[vertex_index_L]

        # print(vertex_list)

        node_list = [LatticeNode(i) for i in vertex_list]
        return LatticeRoute(node_list, *args, **kwargs)

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
            polygon += [self.padding, self.padding, 0]

        ax.plot(polygon[:, 0], polygon[:, 1], "r-")

        return ax

    def plotCrossovers(self, ax, coords: bool):

        if coords:
            crossovers = deepcopy(self.get_crossovers_coords())
            ax.plot(
                crossovers[:, 0], crossovers[:, 1], "bx", ms=2.5,
            )
        else:
            crossovers = deepcopy(self.get_crossovers())
            crossovers = np.ma.masked_where(crossovers == 0, crossovers)
            cmap = plt.get_cmap("cool")
            ax.imshow(crossovers[::-1], cmap)

    def plot(self, lattice: np.ndarray, ax: plt.Axes = None, fout: str = None, poly = True, cross = True, lattice_points = True):

        nodes = np.array(lattice)
        if not ax:
            fig, ax = plt.subplots()
        # plt.grid(True)
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontsize(4)
        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.yaxis.set_major_locator(MultipleLocator(10))
        ax.set_xlabel("No. of nucleotides")
        ax.set_ylabel("No. of strands")

        if np.shape(nodes)[1] == 2:
            print("Plotting from coords")
            if lattice_points:
                ax.plot(nodes[:, 0], nodes[:, 1], "ko", ms=0.5)
            if poly:
                self.plotPolygon(ax, nodes, coords=True)
            if cross:
                self.plotCrossovers(ax, coords=True)
            ax.set_title("Lattice plotted from coords")

        else:
            cmap = plt.get_cmap("hot")
            if lattice_points:
                ax.imshow(nodes[::-1], cmap)
            if poly:
                self.plotPolygon(ax, nodes, coords=False)
            if cross:
                self.plotCrossovers(ax, coords=False)
            ax.set_title("Lattice plotted from array")

        plt.gca().set_aspect(5)
        if fout:
            plt.savefig(f"{fout}.png", dpi=500)
        
        if not ax:
            plt.show()
    

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
    polygon = square * 10
    lattice = Lattice(polygon, straightening_factor=5)
    # lattice.plot(lattice.adjusted_array, fout="Unstraightened")
    # lattice.plot(lattice.adjusted_coords, fout="polygonCoords")
    # lattice.plot(lattice.final_coords)
    # lattice.plot(lattice.final_coords, fout="Crossovers_on_Trapezium_Coords", show=True)

    lattice.plot(lattice.final_coords, fout="Coords")
    lattice.plot(lattice.final_array, fout="Array")

    print("completed")
