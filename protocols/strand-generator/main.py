import numpy as np

from origamiUROP.oxdna import Strand, Nucleotide, System
from origamiUROP.oxdna.strand import BASE_BASE, POS_BASE, FENE_EPS, POS_BACK, POS_STACK, generate_helix
from origamiUROP.tools import get_rotation_matrix

ACROSS = {
    'A' : 'T',
    'C' : 'G',
    'G' : 'C',
    'T' : 'A',
    'a' : 't',
    't' : 'a',
    'c' : 'g',
    'g' : 'c'
}

def get_3p(nuc : Nucleotide, base : str = 'A') -> Nucleotide:
    return

def get_5p(nuc : Nucleotide, base : str = 'T') -> Nucleotide:
    # shift round in a1 direction
    angle = (np.pi / 180.0) * 34.3
    rotation_matrix = get_rotation_matrix(nuc._a3, angle)
    a1 = np.dot(
            nuc._a1,
            rotation_matrix,
        )
    # shift up in a3 direction
    new_base = nuc.pos_base + BASE_BASE * nuc._a3
    new_pos = new_base - a1 * POS_BASE * 1.5
    return Nucleotide(base, new_pos, a1, nuc._a3.copy())

def get_across(nuc : Nucleotide) -> Nucleotide:
    """
    Returns Nucleotide o
    """
    a1 = -nuc._a1
    a3 = -nuc._a3
    pos_com = nuc.pos_com - a1 * (1.85 + 2 * POS_BACK)
    return Nucleotide(ACROSS[nuc._base], pos_com, a1, a3)

def main():
    nucleotides = []
    print("Creating a nucleotide:")
    nucleotides.append(
        Nucleotide(
            'A',
            np.array([0., -10., 0.]),
            a1 = np.array([1., 0., 0.]),
            a3 = np.array([0., 1., 0.])
        )
    )
    print(f"Nucleotide #0: {nucleotides[0]}")
    print("Creating more nucleotides...\n")
    for i in range(40):
        nucleotides.append(get_5p(nucleotides[-1]))
    strand = Strand(nucleotides=nucleotides)
    print(f"Strand: {strand}")
    system = System(np.array([20., 20., 20.]))
    system.add_strand(strand)
    nucleotides = []
    for nuc in strand.nucleotides[::-1]:
        nucleotides.append(get_across(nuc))
    system.add_strand(Strand(nucleotides))
    system.add_strands(
        generate_helix(40, start_pos=np.array([-10., -10., 0.]), double=True)
    )
    system.write_oxDNA('generator')
    return

if __name__ == '__main__':
    print('Running Strand Generator...\n')
    main()