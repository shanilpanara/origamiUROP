import matplotlib.pyplot as plt
import matplotlib.tri as tri

# 1D list containing x-coordinate of each node
nodes_x = [0.000, 1.000, 2.000, 3, 1.000, 1.750, 1.000]
# 1D list containing y-coordinate of each node
nodes_y = [0.000, 0.000, 0.500, 1.000, 1.000, 1.300, 1.700]
scalars = [1.000, 2.000, 1.000, 2.000, 7.000, 4.000, 5.000]
triangles = [
    [4, 3, 0],
]

triangulation = tri.Triangulation(nodes_x, nodes_y, triangles)
plt.triplot(triangulation, "-k")
plt.tricontourf(triangulation, scalars)
plt.colorbar()
plt.show()
