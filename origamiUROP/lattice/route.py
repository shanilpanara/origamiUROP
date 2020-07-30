import numpy as np
import matplotlib.pyplot as plt

from ..tools import DNANode, DNAEdge
from .edge import LatticeEdge
from .node import LatticeNode
from ..oxdna import Strand, System

class LatticeRoute(Strand):
    """
    Strand which follows the path defined by a set of
    LatticeNodes which are joined by LatticeEdges
    """
    def __init__(self):
        # we want to generate the Edges dynamically not the nodes
        # the nodes will be accessible as the list they are
        # contained in
        self.nodes = []
    
    @property
    def edges(self):
        _edges = [LatticeEdge(node, self.nodes[i+1]) for i, node in enumerate(self.nodes[:-1])]
        return _edges

    def plot(self):
        nodes = np.array(self.nodes)
        plt.plot(nodes[:, 0], nodes[:, 1], 'ko')
        for edge in self.edges:
            plt.arrow(edge.vertices[0], edge.vertices[1], edge.vector[0], edge.vector[1])
        plt.show()

    def system(self, **kwargs):
        _system = System(kwargs.get('box', np.array([50., 50., 50.])))
        _system.add_strand(self)
        return _system