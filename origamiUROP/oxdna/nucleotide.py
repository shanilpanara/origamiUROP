import numpy as np

# Emperical oxDNA constants
POS_BACK = -0.4
POS_STACK = 0.34
POS_BASE = 0.4


class Nucleotide:
    """
    Nucleotides compose Strands
    cm_pos --- Center of mass position
        Ex: [0, 0, 0]
    a1 --- Unit vector indicating orientation of backbone with respect to base
        Ex: [1, 0, 0]
    a3 --- Unit vector indicating orientation (tilting )of base with respect to backbone
        Ex: [0, 0, 1]
    base --- Identity of base, which must be designated with either numbers or
        letters (this is called type in the c++ code). Confusingly enough, this
        is similar to Particle.btype in oxDNA.
        
        Number: {0,1,2,3} or any number in between (-inf,-7) and (10, inf)
        To use specific sequences (or an alphabet large than four) one should
        start from the complementary pair 10 and -7. Complementary pairs are
        such that base_1 + base_2 = 3;
        
        Letter: {A,G,T,C} (these will be translated to {0, 1, 3, 2}).
        
        These are set in the dictionaries: number_to_base, base_to_number
    """

    def __init__(
        self,
        cm_pos: np.ndarray,  # position of centre of mass (likely inbetween base/backbone)
        a1: np.ndarray,  # direction of base
        a3: np.ndarray,  # direction of stacking, normal to base
        base: str,  # either 'A' 'C' 'T' or 'G'
        v: np.ndarray = np.array([0.0, 0.0, 0.0]),
        L: np.ndarray = np.array([0.0, 0.0, 0.0]),
    ):
        self.cm_pos = cm_pos
        self._a1 = a1
        self._a3 = a3
        self._v = v
        self._L = L
        self._base = base

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
