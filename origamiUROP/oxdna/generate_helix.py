import numpy as np
from origamiUROP.oxdna import Nucleotide, Strand, System

# Constants
PI = np.pi

# Emperical oxDNA constants
POS_BACK = -0.4
POS_STACK = 0.34
POS_BASE = 0.4

FENE_R0_OXDNA = 0.7525
FENE_R0_OXDNA2 = 0.7564
FENE_EPS = 2.0


# Center of the double strand
CM_CENTER_DS = POS_BASE + 0.2

# Ideal distance of base sites of 2 nts (to be base paired in a duplex)
BASE_BASE = 0.3897628551303122

number_to_base = {0: "A", 1: "G", 2: "C", 3: "T"}


def get_rotation_matrix(axis, anglest):
    """
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


def generateHelix(
    BP: int,
    sequence=None,
    start_pos: np.ndarray = np.array([0.0, 0.0, 0.0]),
    back_orient_a1: np.ndarray = np.array([1.0, 0.0, 0.0]),
    base_orient_a3: np.ndarray = np.array([0.0, 0.0, 1.0]),
    initial_rot: float = 0.0,  # radians
    BP_PER_TURN: float = 10.34,
    ds_start=None,
    ds_end=None,
    double=None,
):
    # Must ensure they aren't given as integers, e.g. [1,0,0] needs to be [1.0,0.0,0.0]
    dir = base_orient_a3.astype("float64")
    perp = back_orient_a1.astype("float64")
    ds_start = 0
    ds_end = BP

    v1 = perp
    # Setup initial parameters
    new_strand_1 = Strand([])
    # and we need to generate a rotational matrix
    R0 = get_rotation_matrix(dir, initial_rot)
    # R = get_rotation_matrix(dir, np.deg2rad(35.9))
    R = get_rotation_matrix(dir, [1, BP])
    a1 = v1
    a1 = np.dot(R0, a1)
    rb = np.array(start_pos)
    a3 = dir

    # Set Sequence
    if sequence == None:
        sequence_numbers = np.random.randint(0, 4, BP)
    elif len(sequence) != BP:
        n = BP - len(sequence)
        sequence_numbers += np.random.randint(0, 4, n)

    sequence_base = [number_to_base[i] for i in sequence_numbers]

    # Add nucleotides in canonical double helix
    for i in range(BP):
        new_strand_1.add_nucleotide(
            Nucleotide(sequence_base[i], rb - CM_CENTER_DS * a1, a1, a3,)
        )
        if i != BP - 1:
            a1 = np.dot(R, a1)
            rb += a3 * BASE_BASE

    if double == True:
        new_strand_2 = Strand()
        for i in reversed(range(ds_start, ds_end)):
            # Note that the complement strand is built in reverse order
            nt = new_strand_1._nucleotides[i]
            a1 = -nt._a1
            a3 = -nt._a3
            nt2_cm_pos = -(FENE_EPS + 2 * POS_BACK) * a1 + nt.pos_com
            reverse_seq_number = [3 - sequence_numbers[j] for j in sequence_numbers]
            reverse_seq = [number_to_base[k] for k in reverse_seq_number]
            new_strand_2.add_nucleotide(Nucleotide(reverse_seq[i], nt2_cm_pos, a1, a3))
        return [new_strand_1, new_strand_2]
    else:
        return new_strand_1
