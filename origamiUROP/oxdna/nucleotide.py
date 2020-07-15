import numpy as np
import pandas as pd

# Emperical oxDNA constants
POS_BACK = -0.4
POS_STACK = 0.34
POS_BASE = 0.4

class Nucleotide:
    """
    Nucleotides compose Strands

    Parameters:
        base - 'A', 'T', 'C' or 'G'
        cm_pos - Center of mass position vector
        a1 - Unit vector indicating orientation of backbone with respect to base
        a3 - Unit vector indicating orientation (tilting) of base with respect to backbone
        v - Linear velocity vector
        L - Angular velocity vector

    Attributes/Properties:



    """

    def __init__(
        self,
        base : str,
        pos_com : np.ndarray,
        a1 : np.ndarray,
        a3 : np.ndarray,
        v : np.ndarray = np.array([0.0, 0.0, 0.0]),
        L : np.ndarray = np.array([0.0, 0.0, 0.0]),
    ):
        self.pos_com = pos_com
        self._a1 = a1
        self._a3 = a3
        self._v = v
        self._L = L
        self._base = base

        # these are accessed when the nucleotide is added
        # to an oxdna.Strand._nucleotides object
        self._strand_index = -1
        self._before = -1
        self._after = -1

    @property
    def pos_base(self):
        """
        Returns the position of the base site
        """
        return self.pos_com + self._a1 * POS_BASE

    @property
    def pos_stack(self):
        """
        Returns the position the stacking site
        """
        return self.pos_com + self._a1 * POS_STACK

    @property
    def pos_back(self):
        """
        Returns the position of the backbone site
        """
        return self.pos_com + self._a1 * POS_BACK

    @property  # although this wasn't a property before, not sure why?
    def pos_back_rel(self):
        """
        Returns the position of the backbone centroid relative to the centre of mass
        i.e. it will be a vector pointing from the c.o.m. to the backbone
        """
        return self.pos_back - self.pos_com

    @property
    def _a2(self):
        return np.cross(self._a3, self._a1)

    @property
    def series(self) -> pd.Series:
        """
        Writes a pd.Series object containing the information
        needed for writing a row in pd.DataFrame that will
        be used for writing to file.
        """
        return pd.Series({
            'position' : self.pos_com,
            'a1' : self._a1,
            'a3' : self._a3,
            'v' : self._v,
            'L' : self._L,
            'base' : self._base,
            'strand' : self._strand,
            '3p' : self._before,
            '5p' : self._after
        })

