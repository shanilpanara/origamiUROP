import numpy as np

from origamiUROP.oxdna.strand import generate_helix, POS_BACK, Strand
from origamiUROP.oxdna import System

FENE_LENGTH = 0.76

def main():
    # generate a strand
    strands = []
    strand = generate_helix(n=20)[0]
    strands.append(strand.copy())

    last_nuc = strands[0].nucleotides[-1]
    direction = -last_nuc._a3
    a1 = -last_nuc._a1

    # ensure the backbone position is FENE_LENGTH away from
    # the backbone position of the previous nucleotide
    start = last_nuc.pos_back + (FENE_LENGTH - POS_BACK) * a1

    # generate strand above that's going in opposite direction
    strand = generate_helix(
        n=20,
        start_position=start,
        direction=direction,
        a1=a1
    )[0]
    strands.append(strand)

    # using the two previously created strands create a new strand that
    # we will add to the system
    strand = Strand(nucleotides=strands[0].nucleotides+strands[1].nucleotides)

    # create system and add the final completed strand
    system = System(np.array([50., 50., 50.]))
    system.add_strand(strand)
    system.write_oxDNA('turns')

if __name__ == '__main__':
    main()