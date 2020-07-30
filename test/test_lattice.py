import numpy as np

from origamiUROP.lattice import LatticeNode, LatticeEdge, LatticeRoute

def test_node():
    node = LatticeNode(np.array([0., 1., 0.]))
    assert node.angle == None
    node.vector_3p = np.array([0., 1., 0.])
    node.vector_5p = np.array([1., 0., 0.])
    assert node.angle == np.pi / 2
    return

def test_edge():
    nodes = [
        LatticeNode(np.array([0., 1., 0.])),
        LatticeNode(np.array([0., 3., 0.]))
    ]
    edge = LatticeEdge(*nodes)
    assert edge.length == 2.0
    return

def test_route():
    return

if __name__=='__main__':
    test_node()
    test_edge()
    test_route()