#!/usr/bin/env python
"""
pytest module for origamiUROP.polygons
"""
from os import path

import numpy as np

from origamiUROP.polygons import Edge, BoundaryPolygon

ROOT = "/".join(path.abspath(__file__).split("/")[:-1])


def test_Edge():
    edge = Edge(np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, -1.0]))
    assert edge.length == 2
    assert edge.kind == "boundary"
    assert np.linalg.norm(edge.vector - np.array([0, 0, -2])) == 0.0
    assert np.linalg.norm(edge.unit_vector - np.array([0, 0, -1])) == 0.0
    return edge


def test_BoundaryPolygon():
    square = BoundaryPolygon(np.array([[1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 0],]))

    assert square.n_edges == len(square.edges)
    assert square.x.sum() == 2
    assert square.y.sum() == 2
    assert square.z.sum() == 0
    assert isinstance(square.edges[0], Edge)

    # square.plot2D(show=False)
    square.write_STL(f"{ROOT}/square.stl")
    square.write_PLY(f"{ROOT}/square.ply")


if __name__ == "__main__":
    test_Edge()
    test_BoundaryPolygon()
