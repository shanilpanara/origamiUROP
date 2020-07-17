import numpy as np
import pandas as pd
from origamiUROP.oxdna import Nucleotide
import re
from copy import deepcopy


class Strand:
    """
    Collection of nucleotides in the 3' -> 5' direction
    Can be added to the system using SystemObject.add_strand
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
