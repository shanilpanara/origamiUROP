import numpy as np

class Nucleotide():
    def __init__(
        self, 
        cm_pos : np.ndarray,
        a1 : np.ndarray, 
        a3 : np.ndarray,
        v : np.ndarray = np.array([0., 0., 0.]),
        L : np.ndarray = np.array([0., 0., 0.])
    ):
        self.cm_pos = cm_pos
        self._a1 = a1
        self._a3 = a3
        self._v = v
        self._L = L