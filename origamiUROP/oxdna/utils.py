import numpy as np
import scipy.linalg as la

def get_rotation_matrix(axis : np.ndarray, theta : float) -> np.ndarray:
    """Returns a rotation matrix for rotating an angle theta
    (rad) around an axis vector
    """
    return la.expm(np.cross(np.eye(3), axis/la.norm(axis)*theta))
    