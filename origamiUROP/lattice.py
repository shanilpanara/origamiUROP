import numpy as np
from shapely.geometry import MultiPoint
from shapely import geometry
import matplotlib.pyplot as plt

from operator import itemgetter
from typing import List
import pprint
from copy import deepcopy

# CONSTANTS - in OXDNA units
x_spacing = 0.34  # POS_STACK
y_spacing = 1  # default value for now

# class lattice:
#     def __init__(self, polygon: geometry.Polygon):
#         self.polygon = polygon

#     @property
#     def polygon_bounds(self): ->


# @property
# def ...
square = np.array([[0, 0, 0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
trapezium = np.array([[0, 0, 0], [1.0, 10.0, 0.0], [2.0, 10.0, 0.0], [3.0, 0.0, 0.0]])

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


def plot_lattice(lattice, title: str):
    fig, ax = plt.subplots()
    plt.imshow(lattice)
    plt.gca().invert_yaxis()
    plt.xlabel("Width of rows (bp)")
    plt.ylabel("Rows")
    plt.title(title)


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


polygon = geometry.Polygon(trapezium * 3)
xmin, ymin, xmax, ymax = polygon.bounds
xmin = round_to_multiple(xmin, x_spacing)
xmax = round_to_multiple(xmax, x_spacing)
ymin = round_to_multiple(ymin, y_spacing)
ymax = round_to_multiple(ymax, y_spacing)


# create a big grid of points spaced with "(x/y)_spacing"
x_grid = np.linspace(xmin, xmax, int((xmax - xmin) / x_spacing) + 2)
y_grid = np.linspace(ymin, ymax, int((ymax - ymin) / y_spacing) + 1)
bounds_grid = np.transpose(
    [np.tile(x_grid, len(y_grid)), np.repeat(y_grid, len(x_grid))]
)
points = MultiPoint(bounds_grid)

# intersect grid with polygon
result = points.intersection(polygon)
# x, y = [], []
# for point in result.geoms:
#     x.append(round(point.x / x_spacing, 2))
#     y.append(round(point.y / y_spacing, 2))
# x = np.array(x)
# y = np.array(y)

# print(x.max().astype(int))

# lattice = np.zeros((+1, +1)).astype(int)


# plt.plot(x, y, "ob")
# plt.gca().set_aspect("equal")


xy = np.array(result.__geo_interface__["coordinates"])

# Normalise values
xy[:, 0] /= x_spacing
xy[:, 1] /= y_spacing
xy = np.around(xy).astype(np.uint64)

# Set the lowest values to 0
x_min = xy[:, 0].min()
y_min = xy[:, 1].min()
xy[:, 0] -= x_min
xy[:, 1] -= y_min

# Make numpy grid of 0's
x_max = int(xy[:, 0].max())
y_max = int(xy[:, 1].max())
lattice = np.zeros((y_max + 1, x_max + 1)).astype(int)

# Sort values by row (y axis)
# xy = sorted(xy, key=itemgetter(1))

# Make quantised lattice
for point in xy:
    lattice[point[1], point[0]] = 1
lattice = np.pad(lattice, pad_width=5, mode="constant", constant_values=0)

plot_lattice(lattice, "Initial Quantised Polygon")

# No. of nt per row
nt_per_row = []
for row in range(len(lattice)):
    nt_per_row.append(np.sum(lattice[row]))
nt_per_row = np.array(nt_per_row)

rounded_to_32 = lambda x: round_to_multiple(x, 32, 0)
nt_per_row_round = np.array(list(map(rounded_to_32, nt_per_row)))
nt_per_row_diff = nt_per_row_round - nt_per_row

lattice_edit = deepcopy(lattice)
# for row in range(len(lattice_edit)):
i = len(nt_per_row_diff)
j = 0
for row, value in enumerate(nt_per_row_diff):
    if value != 0:
        # find location of the first 1 from the right in lattice
        # i.e. one_right = x-coord (x location of one from the right)
        if np.sum(lattice_edit[row]) == 0:
            value = 0
            print("breaking code")
        one_right = int(np.argwhere(lattice_edit[row])[-1])

        # one_left = int(np.argwhere(lattice_edit[40])[0])
        if value > 0:  # we need to add 1's
            lattice_edit[row, one_right + 1] = 1
            value -= 1
            # assign a 1 to the value on the right of it
        elif value < 0:
            lattice_edit[row, one_right] = 0
            value += 1
    else:
        pass
    j += 1
    print(f"Complete: {round((i-j)/i*100, 2)}% left")

plot_lattice(lattice_edit, "Second Quantised Polygon")

plot_lattice(lattice_edit - lattice, "Difference")


# fig = plt.figure()
# ax = fig.add_subplot(111)

# plt.matshow(lattice)
# plt.gca().invert_yaxis()
# plt.xlabel('Strand length (bp/nt)')


# This is so badly documented, will need to sift through these datapoints in a logical way somehow...
