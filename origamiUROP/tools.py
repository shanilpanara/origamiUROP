from .oxdna.strand import Strand, generate_helix

class DNANode:
    """
    Abstract class for use with DNAEdge that helps to determine
    angles and vectors needed to generate stable structures.
    """
    def __init__(self):
        pass

class DNAEdge:
    """
    Abstract class that allows subclasses to generate oxDNA Strand
    instances along a vector.
    """
    def __init__(self, vertex_1 : np.ndarray, vertex_2 : np.ndarray):
        self.vertices = np.array([vertex_1, vertex_2])

    def strand(self, sequence: str = None, **kwargs) -> List[Strand]:

        if not sequence:
            # in future version, this will not be so
            # straightforward
            no_of_nucleotides_in_edge = int(self.length)

        else:
            no_of_nucleotides_in_edge = len(sequence)

        strands = generate_helix(
            bp=no_of_nucleotides_in_edge,
            sequence=sequence,
            start_pos=self.vertices[0],
            back_orient_a1=self.perp_vector,
            base_orient_a3=self.unit_vector,
            **kwargs,
        )
        return strands

    def segments(self) -> float:
        return

    def node(node_3p : DNANode = None, node_5p : DNANode = None) -> DNANode:
        """
        Returns a DNANode for the opposite end of the DNANode provided in 
        parameters
        """
        if not (node_3p or node_5p) or (node_3p and node_5p):
            raise TypeError(
                "Only give one node which is at the"
                " 3' or 5' end of the Edge"
            )
        if node_3p:
            pass
        elif node_5p:
            pass
        else:
            raise TypeError("Shouldn't get to this point")

    @property
    def length(self):
        return np.linalg.norm(self.vertices[1] - self.vertices[0])

    @property
    def vector(self) -> np.ndarray:
        return self.vertices[1] - self.vertices[0]

    @property
    def unit_vector(self) -> np.ndarray:
        return self.vector / (self.vector ** 2).sum() ** 0.5

    @property
    def perp_vector(self) -> np.ndarray:
        """Perpendicular vector which lies in the xy plane"""
        return np.cross(self.unit_vector, np.array([0, 0, 1]))