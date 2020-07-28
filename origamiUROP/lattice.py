import numpy as np
from shapely.geometry import MultiPoint
from shapely import geometry

# CONSTANTS
x_spacing = 0.34  # POS_STACK
y_spacing = 1.00  # default value for now

# class lattice:
#     def __init__(self, lattice_box):
#         self.box = lattice_box

#     @property
#     def ...


def round_to_multiple(n, multiple_of=0.34, decimal_places=2):
    """
    rounds to nearest multiple of value given, 
    """
    mo = multiple_of

    # Smaller multiple
    a = (n // mo) * mo

    # Larger multiple
    b = a + mo

    # Return of closest of two
    closest_multiple = b if n - a > b - n else a
    return round(closest_multiple, decimal_places)


trapezium = np.array([[0, 0, 0], [1.0, 10.0, 0.0], [11.0, 10.0, 0.0], [15.0, 0.0, 0.0]])
p = geometry.Polygon(trapezium)
xmin, ymin, xmax, ymax = p.bounds


xmin = round_to_multiple(xmin, 0.34)
xmax = round_to_multiple(xmax, 0.34)
ymin = round_to_multiple(ymin, 1.0)
ymax = round_to_multiple(ymax, 1.0)

print(int((xmax - xmin) / x_spacing))
print(int((ymax - ymin) / y_spacing))
x = np.linspace(xmin, xmax, int((xmax - xmin) / x_spacing) + 1)
y = np.linspace(ymin, ymax, int((ymax - ymin) / y_spacing) + 1)
bounds_grid = np.transpose([np.tile(x, len(y)), np.repeat(y, len(x))])
# create a big grid of points
points = MultiPoint(bounds_grid)
result = points.intersection(p)

xy = np.array(result.__geo_interface__["coordinates"])
# This is so badly documented, will need to sift through these datapoints in a logical way somehow...
