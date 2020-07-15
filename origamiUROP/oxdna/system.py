import numpy as np


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
        self, box: np.ndarray, time: int = 0, E_pot: float = 0.0, E_kin: float = 0.0
    ):

        self.box = box
        self.time = time
        self.E_pot = E_pot
        self.E_kin = E_kin

        self._strands = []

    @property
    def E_tot(self):
        return self.E_pot + self.E_kin

    @property
    def strands(self) -> list:
        for i, strand in enumerate(self._strands):
            strand.index = i
        return self._strands

    @property
    def nucleotides(self) -> list:
        result = []
        for strand in self.strands:
            result += strand._nucleotides
        return result

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

        with open(f'{prefix}.top', 'w') as f:
            f.write(f'{len(self.nucleotides)} {len(self.strands)}\n')

