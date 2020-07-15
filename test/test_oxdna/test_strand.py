from origamiUROP.oxdna import Strand, Nucleotide
import numpy as np
import pandas as pd

def test_Strand():
    strand = Strand([
        Nucleotide('A', np.array([1., 0., 0.]), np.array([1., 0., 0.]), np.array([0, 0., 1.])),
        Nucleotide('A', np.array([2., 0., 0.]), np.array([1., 0., 0.]), np.array([0, 0., 1.])),
        Nucleotide('A', np.array([3., 0., 0.]), np.array([1., 0., 0.]), np.array([0, 0., 1.])),
        Nucleotide('A', np.array([4., 0., 0.]), np.array([1., 0., 0.]), np.array([0, 0., 1.])),
        Nucleotide('A', np.array([5., 0., 0.]), np.array([1., 0., 0.]), np.array([0, 0., 1.])),
        Nucleotide('A', np.array([6., 0., 0.]), np.array([1., 0., 0.]), np.array([0, 0., 1.])),
    ])
    print(strand)
    print(strand.dataframe)
    return

if __name__ == '__main__':
    test_Strand()