import numpy as np
import pandas as pd

class Strand:
    """
    Collection of nucleotides in the 3' -> 5' direction
    Can be added to the system using SystemObject.add_strand
    """

    def __init__(self, nucleotides: list = []):
        self._nucleotides = nucleotides

        self.index = -1

    def __repr__(self):
        return f'3>[{self.index}]>{self.sequence}>>>5'

    @property
    def sequence(self) -> str:
        return ''.join([i._base for i in self._nucleotides])

    @property
    def nucleotides(self) -> list:
        """
        Returns the list of nucleotides where each nucleotide
        is 
        """
        for i, nucleotide in self._nucleotides:
            nucleotide._before = i-1
            nucleotide._after = i+1
            nucleotide._strand_index = self.index
        return self._nucleotides

    @property
    def dataframe(self) -> pd.DataFrame:
        """
        Returns a pd.Dataframe which can be used to write
        the information for both the configuration and
        topology files.
        """
        result = pd.DataFrame(
            [i.series for i in self.nucleotides]
        )
        return result

    def __len__(self) -> int:
        return len(self._nucleotides)


    ## this is useful but we will probably change this (also we will use a property setter decorator)
    def set_sequence(self, seq: str):  # enforcing we don't use numbers, just letters :)
        if len(seq) != len(self._nucleotides):
            print("error, seq not the same length as no. of nucleotides")
        self._sequence = seq

    ## this isn't needed
    @property
    def list_sequence(self) -> list:
        return self._sequence

    ## this has be implemented as @property\\def sequence(self)
    @property
    def str_sequence(self) -> str:
        return "".join(self._sequence)
