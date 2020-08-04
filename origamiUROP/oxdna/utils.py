import numpy as np

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