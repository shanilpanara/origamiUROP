import numpy as np
import matplotlib.pyplot as plt

from origamiUROP.polygons import Polygon

def square(L : int) -> np.ndarray:
    return np.array([
        [0., 0.],
        [1., 0.],
        [1., 1.],
        [0., 1.],
    ]) * L

def pentagon(L : int) -> np.ndarray:
    return (np.array([
        [0.550, 0.450],
        [0.455, 0.519],
        [0.491, 0.631],
        [0.609, 0.631],
        [0.645, 0.519]
    ]) - 0.5) * L

class PointInPolygon2D:
    def __init__(self, polygon_array, point, gridsize=0.1):
        self.polygon = np.concatenate([polygon_array, np.array([polygon_array[0, :]])], axis=0)
        self.point = point
        self.obj = Polygon(polygon_array)

    @property
    def box(self):
        return np.array([
            self.polygon[:, 0].min(),
            self.polygon[:, 0].max(),
            self.polygon[:, 1].min(),
            self.polygon[:, 1].max()
        ])

    @property
    def in_box(self):
        box = self.box
        return (
            self.point[0] > box[0] and \
            self.point[0] < box[1] and \
            self.point[1] > box[2] and \
            self.point[1] < box[3]
        )

    @property
    def in_polygon(self):
        if not self.in_box:
            return False
        # do ray casting
        intersections = 0

        # choose a point outside the bounding box - keep it simple,
        # use the y-axis of the point
        start = np.array([
            self.box[0], # start from LHS of box
            self.point[1] # stay on y=A axis
        ])

        # create a vector equation of the line
        # connecting the points
        vector = np.array([
            self.box[0] - self.point[0], 
            0.
        ])

        # equation = start + lambda * vector

        # use linear algebra to check how many intersections
        # there are
        return True

    def plot(self):
        fig, ax = plt.subplots()
        colour = 'g' if self.in_polygon else 'r'
        box = np.array([
            [self.box[0], self.box[2]],
            [self.box[1], self.box[2]],
            [self.box[1], self.box[3]],
            [self.box[0], self.box[3]],
            [self.box[0], self.box[2]]
        ])
        ax.plot(self.polygon[:, 0], self.polygon[:, 1], 'k-')
        ax.plot(box[:, 0], box[:, 1], 'k:')
        ax.plot(self.point[0], self.point[1], f'{colour}o')
        plt.show()

def main(**kwargs):
    obj = PointInPolygon2D(pentagon(5), np.array([0.5, 0.5]))
    print(obj.box)
    print(obj.in_box)
    obj.plot()
    return 

if __name__ == '__main__':
    main()