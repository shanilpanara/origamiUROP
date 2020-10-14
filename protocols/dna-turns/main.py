import numpy as np

from origamiUROP.oxdna.strand import generate_helix, POS_BACK, Strand
from origamiUROP.oxdna import System

FENE_LENGTH = 0.76

def main(length=16, n_strands=10):
    # generate a strand
    strands = []
    strand = generate_helix(n=length)[0]
    strands.append(strand.copy())

    for i in range(n_strands-1):

        last_nuc = strands[-1].nucleotides[-1]
        direction = -last_nuc._a3
        a1 = -last_nuc._a1

        # ensure the backbone position is FENE_LENGTH away from
        # the backbone position of the previous nucleotide
        start = last_nuc.pos_back + (FENE_LENGTH - POS_BACK) * a1

        # generate strand above that's going in opposite direction
        strand = generate_helix(
            n=length,
            start_position=start,
            direction=direction,
            a1=a1
        )[0]
        strands.append(strand)

    # using the two previously created strands create a new strand that
    # we will add to the system
    nucleotides = []
    for strand in strands:
        nucleotides += strand.nucleotides
    strand = Strand(nucleotides=nucleotides)

    # create system and add the final completed strand
    system = System(np.array([50., 50., 50.]))
    system.add_strand(strand)
    system.write_oxDNA('turns')

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-l', '--length', default=16, required=False)
    parser.add_argument('-n', '--n-strands', default=10, required=False)
    main(**vars(parser.parse_args()))