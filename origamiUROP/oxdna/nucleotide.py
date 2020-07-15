import numpy as np


class Nucleotide:
    def __init__(
        self,
        cm_pos: np.ndarray,
        a1: np.ndarray,
        a3: np.ndarray,
        v: np.ndarray = np.array([0.0, 0.0, 0.0]),
        L: np.ndarray = np.array([0.0, 0.0, 0.0]),
    ):
        self.cm_pos = cm_pos
        self._a1 = a1
        self._a3 = a3
        self._v = v
        self._L = L

    @property
    def pos_base(self):
        """
        Returns the position of the base centroid.
        Note that cm_pos is the centroid of the backbone and base.
        """
        return self.cm_pos + self._a1 * POS_BASE

    @property
    def pos_stack(self):
        return self.cm_pos + self._a1 * POS_STACK

    @property
    def pos_back(self):
        """
        Returns the position of the backbone centroid.
        Note that cm_pos is the centroid of the backbone and base.
        """
        return self.cm_pos + self._a1 * POS_BACK

    @property  # although this wasn't a property before, not sure why?
    def pos_back_rel(self):
        """
        Returns the position of the backbone centroid relative to the centre of mass
        i.e. it will be a vector pointing from the c.o.m. to the backbone
        """
        return self.pos_back - self.cm_pos

    @property
    def _a2(self):
        return np.cross(self._a3, self._a1)
