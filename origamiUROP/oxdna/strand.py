"""
Storage and Management of oxDNA strands via the Strand
class
"""
import numpy as np
import pandas as pd

from origamiUROP.oxdna import Nucleotide
import re
from copy import deepcopy

# Constants
PI = np.pi

# Emperical oxDNA constants
POS_BACK = -0.4
POS_STACK = 0.34
POS_BASE = 0.4

FENE_EPS = 2.0

# Center of the double strand
CM_CENTER_DS = POS_BASE + 0.2

# Ideal distance of base sites of 2 nts (to be base paired in a duplex)
BASE_BASE = 0.3897628551303122

number_to_base = {0: "A", 1: "G", 2: "C", 3: "T"}
base_to_number = {"A": 0, "G": 1, "C": 2, "T": 3}


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


class Strand:
    """
    Collection of nucleotides in the 3' -> 5' direction
    Can be added to the system using SystemObject.add_strand

    Parameters:
        nucleotides ([]) - list of nucleotides to form the strand
    
    Attributes/Properties:
        index - strand_id
        sequence - string sequence of bases
        nucleotides - list of Nucleotides
        dataframe - table of values for top and conf files
        copy - get a copy of a Strand instance
    """

    def __init__(self, nucleotides: list = []):
        self._nucleotides = nucleotides

        self.index = 1
        self._nucleotide_shift = 0

    def __repr__(self):
        return f"3>[{self.index}]>{self.sequence}>>5"

    @property
    def sequence(self) -> str:
        return "".join([i._base for i in self._nucleotides])

    @property
    def nucleotides(self) -> list:
        """
        Returns the list of nucleotides where each nucleotide
        is 
        """
        for i, nucleotide in enumerate(self._nucleotides):
            if i == 0:
                nucleotide._before = -1
            else:
                nucleotide._before = i - 1 + self._nucleotide_shift
            if i == len(self._nucleotides) - 1:
                nucleotide._after = -1
            else:
                nucleotide._after = i + 1 + self._nucleotide_shift
            nucleotide._strand_index = self.index
        return self._nucleotides

    def add_nucleotide(self, Nucleotide: Nucleotide):
        self._nucleotides.append(Nucleotide)
        # print(f"Added nucleotide to the strand, now of length {len(self._nucleotides)}")

    @property
    def dataframe(self) -> pd.DataFrame:
        """
        Returns a pd.Dataframe which can be used to write
        the information for both the configuration and
        topology files.
        """
        result = pd.DataFrame([i.series for i in self.nucleotides])
        return result

    def __len__(self) -> int:
        return len(self._nucleotides)

    @sequence.setter
    def sequence(self, seq: str):
        seq = seq.upper()
        if len(seq) > len(self._nucleotides):
            raise ValueError(
                f"Sequence length must = strand length {len(self._nucleotides)}"
            )
        if re.findall("[^ACGT]", seq):
            raise ValueError("Sequence can only contain A, G, C or T")

        for i, base in enumerate(seq):
            self._nucleotides[i]._base = base

    # When a strand is added using add_strand in system.py
    # We create a new object, so the strand is not accessible
    # from anywhere other than the system
    def copy(self):
        return deepcopy(Strand(self._nucleotides))


def generateHelix(
    bp: int,
    sequence: str = None,
    start_pos: np.ndarray = np.array([0.0, 0.0, 0.0]),
    back_orient_a1: np.ndarray = np.array([1.0, 0.0, 0.0]),
    base_orient_a3: np.ndarray = np.array([0.0, 1.0, 0.0]),
    initial_rot: float = 0.0,  # radians
    BP_PER_TURN: float = 10.34,
    ds_start: int = None,
    ds_end: int = None,
    double: bool = False,
) -> list:
    """
    Generate a strand of DNA around a centerline (a3)
        - ssDNA (default) or dsDNA (double = True)
        - by default, strand propagates in the xy plane

        Arguments:
        bp --- Integer number of bp/nt (required)
        sequence --- String. Should be same length as bp (default None)
            Default (None) generates a random sequence. If sequence shorter
            than bp given, a random sequence is assigned for the remaining
            nucleotides
        start_pos --- Location to begin building the strand 
            (default np.array([0.0, 0.0, 0.0]))
        back_orient_a1 --- Sets a1 unit vector, indicating orientation 
            (tilting) of base with respect to backbone 
            (default np.array([1.0, 0.0, 0.0]))
        base_orient_a3 --- a3 unit vector, indicating orientation of 
            backbone with respect to base
            (default np.array([0.0, 1.0, 0.0]))
        initial_rot --- Rotation of first bp in radians 
            (default 0.0)
        BP_PER_TURN --- Base pairs per complete 2*pi helix turn.
            (default 10.34)
        ds_start --- Index (from 0) to begin double stranded region
            (default None)
        ds_end --- Index (from 0) to end double stranded region 
            (default None)
        double --- Generate dsDNA by forming complementary strand in the
            reverse direction (default False)
        
    """
    # Set Sequence
    if sequence is not None:
        try:
            assert type(sequence) == str
        except TypeError:
            raise TypeError("Sequence must be given as a string")
        if len(sequence) > bp:
            n = len(sequence) - bp
            # print(f"Final {n} bases will not be assigned, seq. too long")
            sequence_base = sequence
        elif len(sequence) < bp:
            n = bp - len(sequence)
            extra_seq_to_add = np.random.randint(0, 4, n)
            extra_seq_to_add = "".join(str(number_to_base[x]) for x in extra_seq_to_add)
            sequence_base = sequence + extra_seq_to_add
        else:
            sequence_base = sequence
    else:
        sequence_numbers = np.random.randint(0, 4, bp)
        sequence_base = [number_to_base[i] for i in sequence_numbers]

    # Ensure vectors are in [1.0,0.0,0.0] format (not [1,0,0])
    dir = base_orient_a3.astype("float64")
    perp = back_orient_a1.astype("float64")
    updated_pos = np.array(start_pos.astype("float64"))

    # Setup
    R0 = get_rotation_matrix(dir, initial_rot)
    R = get_rotation_matrix(dir, [1, bp])
    a1 = np.dot(R0, perp)
    a3 = dir

    ds_start = 0
    ds_end = bp

    # Add nucleotides in canonical double helix
    new_strand_1 = Strand([])
    for i in range(bp):
        new_nt = Nucleotide(sequence_base[i], updated_pos - CM_CENTER_DS * a1, a1, a3,)
        new_strand_1.add_nucleotide(new_nt)
        if i != bp - 1:
            a1 = np.dot(R, a1)
            updated_pos += a3 * BASE_BASE

    if double == True:
        new_strand_2 = Strand()
        sequence_numbers = [base_to_number[i] for i in sequence_base]
        reverse_seq_number = [3 - j for j in sequence_numbers]
        reverse_seq = [number_to_base[k] for k in reverse_seq_number]
        for i in reversed(range(ds_start, ds_end)):
            # Note that the complement strand is built in reverse order
            nt1 = new_strand_1._nucleotides[i]
            a1 = -nt1._a1
            a3 = -nt1._a3
            nt2_pos_com = -(FENE_EPS + 2 * POS_BACK) * a1 + nt1.pos_com
            new_strand_2.add_nucleotide(Nucleotide(reverse_seq[i], nt2_pos_com, a1, a3))
        return [new_strand_1, new_strand_2]
    else:
        return [new_strand_1]
