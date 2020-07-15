import re

import numpy as np
import pandas as pd

CONFIGURATION_COLUMNS = ['position', 'a1', 'a3', 'v', 'L']
TOPOLOGY_COLUMNS = ['base', 'strand', '3p', '5p']

def oxDNA_string(dataframe : pd.DataFrame) -> str:
    output = dataframe.to_string(
        header=False, 
        index=False, 
        justify='left',
        formatters= {
            'position' : lambda x: [f"{i:.4f}" for i in x],
            'a1' : lambda x: [f"{i:.4f}" for i in x],
            'a3' : lambda x: [f"{i:.4f}" for i in x],
            'v' : lambda x: [f"{i:.4f}" for i in x],
            'L' : lambda x: [f"{i:.4f}" for i in x],
        }
    )
    output = re.sub(r"\[|\]|\'|\`|\,", "", output)
    output = output.strip()
    output = output.replace('\n ','\n')
    output += '\n'
    output = output.replace('  ', ' ')
    return output.replace('  ', ' ')

class System:
    """
    Object representing an oxDNA system
    Contains strands
    Arguments:
    box -- the box size of the system
        Ex: box = [50, 50, 50]
    time --- Time of the system
    E_pot --- Potential energy
    E_kin --- Kinetic energy
    """

    def __init__(
        self,
        box : np.ndarray, 
        time : int = 0, 
        E_pot : float = 0.0, 
        E_kin : float = 0.0
    ):

        self.box = box
        self.time = time
        self.E_pot = E_pot
        self.E_kin = E_kin

        self._strands = []

    def __repr__(self) -> str:
        return f"oxDNASystem[strands: {len(self._strands)}, nucleotides: {len(self.nucleotides)}]"

    @property
    def E_tot(self):
        return self.E_pot + self.E_kin

    @property
    def strands(self) -> list:

        # used to set Strand._nucleotide_shift
        shift = 0 
        for i, strand in enumerate(self._strands):
            strand._nucleotide_shift = shift
            strand.index = i
            shift += len(strand)
        return self._strands

    @property
    def nucleotides(self) -> list:
        result = []
        for strand in self.strands:
            result += strand._nucleotides
        return result

    @property
    def dataframe(self) -> pd.DataFrame:
        return pd.concat([i.dataframe for i in self.strands])

    @property
    def configuration(self) -> pd.DataFrame:
        return self.dataframe[CONFIGURATION_COLUMNS]

    @property
    def topology(self) -> pd.DataFrame:
        return self.dataframe[TOPOLOGY_COLUMNS]

    def write_oxDNA(self, prefix: str = "out"):
        """
        Writes two files *.conf and *.top for the
        configuration file and topology file required
        to run a simulation using oxDNA
        """
        with open(f"{prefix}.conf", "w") as f:
            f.write(f"t = {self.time}\n")
            f.write(f"b = {self.box[0]} {self.box[1]} {self.box[2]}\n")
            f.write(f"E = {self.E_pot} {self.E_kin} {self.E_tot}\n")
            f.write(oxDNA_string(self.configuration))

        with open(f'{prefix}.top', 'w') as f:
            f.write(f'{len(self.nucleotides)} {len(self.strands)}\n')
            f.write(oxDNA_string(self.topology))