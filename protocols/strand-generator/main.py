import numpy as np

from origamiUROP.oxdna import Strand, Nucleotide
from origamiUROP.tools import get_rotation_matrix

def get_pos_com(nuc : Nucleotide) -> np.ndarray:
    # shift round in a2 direction
    a2 = nuc._a2
    # shift up in a3 direction
    a3 = nuc._a3
    return

def main():
    nucleotides = []
    print("Creating a nucleotide:\n")
    nucleotides.append(
        Nucleotide(
            'A',
            np.array([0., 0., 0.]),
            a1 = np.array([1., 0., 0.]),
            a3 = np.array([0., 1., 0.])
        )
    )
    print(f"Nucleotide #0: {nucleotides[0]}")
    print("Creating more nucleotides...\n")
    nucleotides.append(
        Nucleotide(
            'T',
            get_pos_com(nucleotides[-1]),
            a1 = get_rotational_matrix(nucleotides[-1]._a1, [1, 'bp'])
            a3 = np.array([0., 1., 0.])
        )
    )
    return

if __name__ == '__main__':
    print('Running Strand Generator...\n')
    main()