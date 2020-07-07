from stl import mesh
import math
import numpy

# Create 3 faces of a cube
data = numpy.zeros(6, dtype=mesh.Mesh.dtype)

# Top of the cube
data["vectors"][0] = numpy.array([[0, 1, 1], [1, 0, 1], [0, 0, 1]])
data["vectors"][1] = numpy.array([[1, 0, 1], [0, 1, 1], [1, 1, 1]])
# Front face
# data["vectors"][2] = numpy.array([[1, 0, 0], [1, 0, 1], [1, 1, 0]])
# data["vectors"][3] = numpy.array([[1, 1, 1], [1, 0, 1], [1, 1, 0]])
# Left face
# data["vectors"][4] = numpy.array([[0, 0, 0], [1, 0, 0], [1, 0, 1]])
# data["vectors"][5] = numpy.array([[0, 0, 0], [0, 0, 1], [1, 0, 1]])

# Since the cube faces are from 0 to 1 we can move it to the middle by
# substracting .5
data["vectors"] -= 0.5

meshes = mesh.Mesh(data.copy())

# Rotate 90 degrees over the Y axis
meshes.rotate([0.0, 0.5, 0.0], math.radians(90))

# Optionally render the rotated cube faces
from matplotlib import pyplot
from mpl_toolkits import mplot3d

# Create a new plot
figure = pyplot.figure()
axes = mplot3d.Axes3D(figure)

# Render the cube faces
axes.add_collection3d(mplot3d.art3d.Poly3DCollection(meshes.vectors))

# Auto scale to the mesh size
scale = numpy.concatenate(meshes).flatten("A")
axes.auto_scale_xyz(scale, scale, scale)

# Show the plot to the screen
pyplot.show()
