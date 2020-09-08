import numpy as np
import scipy.linalg as la

# Debesh's oxDNA constants
SHIFT_BASE = 0.5
SHIFT_ACROSS = 0.56 - (SHIFT_BASE * 0.9)
SHIFT_ROUND = -0.105

def get_rotation_matrix(axis : np.ndarray, theta : float) -> np.ndarray:
    """Returns a rotation matrix for rotating an angle theta
    (rad) around an axis vector
    """
    return la.expm(np.cross(np.eye(3), axis/la.norm(axis)*theta))
    
def next_5p(
    position: np.ndarray,
    a1: np.ndarray, 
    a3: np.ndarray,
    angle: float,
    rise: float
):
    """Returns the pos_com for a new nucleotide in the 5' direction
    """
    # shift round in a1 direction
    rotation_matrix = get_rotation_matrix(a3, angle)
    a1 = np.dot(
            rotation_matrix,
            a1,
        )
    
    # shift up in a3 direction
    new_base = position + rise * a3
    new_pos = new_base - a1 * SHIFT_BASE
    new_pos += SHIFT_ROUND * np.cross(a3, a1)
    return new_pos

def get_box() -> np.ndarray:
    pass

def round_to_multiple(n, mo=0.34, decimal_places=2):
    """
    Function rounds to the nearest multiple of value given
    Returns output to (default) 2 decimal places

    Arguments:

        n --- value (integer or float) to round  
        mo --- "multiple_of" is the value (integer or float) which 
        we want to round to a multiple of  
        decimal_places --- no. of decimals to return
    """
    a = (n // mo) * mo  # Smaller multiple
    b = a + mo  # Larger multiple
    closest_multiple = b if n - a > b - n else a  # Return of closest of two
    return round(closest_multiple, decimal_places)
