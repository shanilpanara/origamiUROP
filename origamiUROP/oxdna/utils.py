import numpy as np
import scipy.linalg as la

def get_rotation_matrix(axis : np.ndarray, theta : float) -> np.ndarray:
    """Returns a rotation matrix for rotating an angle theta
    (rad) around an axis vector
    """
    return la.expm(np.cross(np.eye(3), axis/la.norm(axis)*theta))
    
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