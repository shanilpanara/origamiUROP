class Strand:
    """
    Collection of nucleotides in the 3' -> 5' direction
    Can be added to the system using SystemObject.add_strand
    """

    def __init__(self, nucleotides: list = []):
        self._nucleotides = nucleotides
        self._sequence = []

    @property
    def length(self) -> int:
        return len(self._nucleotides)

    def set_sequence(self, seq: str):  # enforcing we don't use numbers, just letters :)
        if len(seq) != len(self._nucleotides):
            print("error, seq not the same length as no. of nucleotides")
        self._sequence = seq

    @property
    def list_sequence(self) -> list:
        return self._sequence

    @property
    def str_sequence(self) -> str:
        return "".join(self._sequence)

    def __repr__(self):
        return f"No. of Nucleotides: {self.length} \nSequence: {self.str_sequence}"
