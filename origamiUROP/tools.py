from typing import List

import numpy as np

from .oxdna.strand import Strand, generate_helix
from .oxdna.nucleotide import Nucleotide

def get_rotation_matrix(axis, anglest):
    """
    Copied from https://github.com/rgatkinson/oxdna/blob/master/UTILS/utils.py 
    The argument anglest can be either an angle in radiants
    (accepted types are float, int or np.float64 or np.float64)
    or a tuple [angle, units] where angle a number and
    units is a string. It tells the routine whether to use degrees,
    radiants (the default) or base pairs turns
    axis --- Which axis to rotate about
        Ex: [0,0,1]
    anglest -- rotation in radians OR [angle, units]
        Accepted Units:
            "bp"
            "degrees"
            "radiants"
        Ex: [np.pi/2] == [np.pi/2, "radians"]
        Ex: [1, "bp"]
    """
    if not isinstance(anglest, (np.float64, np.float32, float, int)):
        if len(anglest) > 1:
            if anglest[1] in ["degrees", "deg", "o"]:
                angle = (np.pi / 180.0) * anglest[0]
                # angle = np.deg2rad (anglest[0])
            elif anglest[1] in ["bp"]:
                # Allow partial bp turns
                angle = float(anglest[0]) * (np.pi / 180.0) * 35.9
                # angle = int(anglest[0]) * (np.pi / 180.) * 35.9
                # Older versions of numpy don't implement deg2rad()
                # angle = int(anglest[0]) * np.deg2rad(35.9)
            else:
                angle = float(anglest[0])
        else:
            angle = float(anglest[0])
    else:
        angle = float(anglest)  # in degrees, I think

    axis = np.array(axis)
    axis /= np.sqrt(np.dot(axis, axis))
    ct = np.cos(angle)
    st = np.sin(angle)
    olc = 1.0 - ct
    x, y, z = axis

    return np.array(
        [
            [olc * x * x + ct, olc * x * y - st * z, olc * x * z + st * y],
            [olc * x * y + st * z, olc * y * y + ct, olc * y * z - st * x],
            [olc * x * z - st * y, olc * y * z + st * x, olc * z * z + ct],
        ]
    )

class DNAObject:
    def __init__(self):
        self.units = 'oxdna'

class DNANode(np.ndarray):
    """
    Abstract class for use with DNAEdge that helps to determine
    angles and vectors needed to generate stable structures.
    """

    def __new__(cls, *args, **kwargs):
        return np.ndarray.__new__(cls, (3)) * 0.0

    def __init__(self, position: np.ndarray):
        self[:] = position[:]
        self._vector_3p = None
        self._vector_5p = None
        self._pos_3p = None
        self._pos_5p = None
        self._a1_3p = None
        self._a1_5p = None
        self._a3_3p = None
        self._a3_5p = None
        self._angle_3p = None
        self._angle_5p = None

    @property
    def angle(self) -> float:
        """
        Returns the angle between the two vectors
        leaving the Node in radians
        """
        if self._vector_3p is None:
            return None
        if self.vector_5p is None:
            return None

        return np.arccos(np.dot(-self._vector_3p, self._vector_5p))

    @property
    def vector_3p(self) -> np.ndarray:
        """
        Vector entering/leaving the node from the 3' direction

        vector_3p and vector_5p follow on from each other which
        means that they cannot start/finish at the same point
        """
        return self._vector_3p

    @vector_3p.setter
    def vector_3p(self, new_vector: np.ndarray):
        self._vector_3p = new_vector

    @property
    def vector_5p(self) -> np.ndarray:
        """
        Vector leaving/entering the node from the 5' direction

        vector_3p and vector_5p follow on from each other which
        means that they cannot start/finish at the same point
        """
        return self._vector_5p

    @vector_5p.setter
    def vector_5p(self, new_vector: np.ndarray):
        self._vector_5p = new_vector

    @property
    def a1_3p(self):
        return self._a1_3p

    @a1_3p.setter
    def a1_3p(self, new_vector : np.ndarray):
        self._a1_3p = new_vector / np.linalg.norm(new_vector)

    @property
    def a1_5p(self):
        return self._a1_3p

    @a1_5p.setter
    def a1_5p(self, new_vector : np.ndarray):
        self._a1_5p = new_vector / np.linalg.norm(new_vector)

    @property
    def a3_3p(self):
        return self._a1_3p

    @a3_3p.setter
    def a3_3p(self, new_vector : np.ndarray):
        self._a3_5p = new_vector / np.linalg.norm(new_vector)

    @property
    def a3_5p(self):
        return self._a1_3p

    @a3_5p.setter
    def a3_5p(self, new_vector : np.ndarray):
        self._a3_5p = new_vector / np.linalg.norm(new_vector)

    @property
    def pos_3p(self):
        return self._pos_3p

    @pos_3p.setter
    def pos_3p(self, new_vector : np.ndarray):
        self._pos_3p = new_vector

    @property
    def pos_5p(self):
        return self._pos_5p

    @pos_5p.setter
    def pos_5p(self, new_vector : np.ndarray):
        self._pos_5p = new_vector

class DNAEdge:
    """
    Abstract class that allows subclasses to generate oxDNA Strand
    instances along a vector.
    """

    def __init__(
        self, vertex_1: (np.ndarray or DNANode), vertex_2: (np.ndarray or DNANode)
    ):
        if not isinstance(vertex_1, DNANode):
            vertex_1 = DNANode(vertex_1)
        if not isinstance(vertex_2, DNANode):
            vertex_2 = DNANode(vertex_2)
        self.vertices = (vertex_1, vertex_2)
        self.vertices[0].vector_5p = self.vector
        self.vertices[0].a3_5p = self.unit_vector
        
        self.vertices[1].vector_3p = -self.vector
        self.vertices[1].a3_3p = self.unit_vector
        self.vertices[1].a1_3p = None

    def strand(self, sequence: str = None, **kwargs) -> List[Strand]:

        if not sequence:
            # in future version, this will not be so straightforward
            no_of_nucleotides_in_edge = self.nt_length
        else:
            no_of_nucleotides_in_edge = len(sequence)
            if len(sequence) >= self.nt_length:
                print(
                    f"FYI: The Length of `sequence` is longer than the max no. of nucleotides "
                    f"that can be contained within this edge, i.e. {self.nt_length} nucleotides"
                )

        if self.vertices[0].a1_5p:
            a1 = self.vertices[0].a1_5p
        else:
            a1 = np.array([1., 0., 0.])

        strands = generate_helix(
            bp=no_of_nucleotides_in_edge,
            sequence=sequence,
            start_pos=self.vertices[0],
            back_orient_a1=a1,
            base_orient_a3=self.unit_vector,
            **kwargs,
        )
        self.vertices[0].update_from_nucleotide(strands[0].nucleotides[0], '3p')
        self.vertices[1].update_from_nucleotide(strands[0].nucleotides[-1], '5p')
        return strands

    def segments(self) -> float:
        return

    def node(self, node_3p: DNANode = None, node_5p: DNANode = None) -> DNANode:
        """
        Returns a DNANode for the opposite end of the DNANode provided in 
        parameters
        """
        if not (node_3p or node_5p) or (node_3p and node_5p):
            raise TypeError(
                "Only give one node which is at the" " 3' or 5' end of the Edge"
            )
        if node_3p:
            pass
        elif node_5p:
            pass
        else:
            raise TypeError("Shouldn't get to this point")

    @property
    def length(self):
        """The length of the edge in oxdna units (i think)"""
        return np.linalg.norm(self.vector)

    @property
    def nt_length(self):
        """The length of the edge in units of nucleotides"""
        return int(self.length * 2.45)

    @property
    def vector(self) -> np.ndarray:
        return self.vertices[1] - self.vertices[0]

    @property
    def unit_vector(self) -> np.ndarray:
        return self.vector / (self.vector ** 2).sum() ** 0.5

    @property
    def perp_vector(self) -> np.ndarray:
        """Perpendicular vector which lies in the xy plane"""
        return np.cross(self.unit_vector, np.array([0, 0, 1]))
