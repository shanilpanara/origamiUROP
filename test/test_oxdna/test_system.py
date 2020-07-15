from origamiUROP.oxdna import System, Strand, Nucleotide
import numpy as np
import pandas as pd

def test_System():
    system = System(np.array([50., 50., 50.]))
    strand_1 = Strand([
        Nucleotide('A', np.array([1., 0., 0.]), np.array([1., 0., 0.]), np.array([0, 0., 1.])),
        Nucleotide('A', np.array([2., 0., 0.]), np.array([1., 0., 0.]), np.array([0, 0., 1.])),
        Nucleotide('A', np.array([3., 0., 0.]), np.array([1., 0., 0.]), np.array([0, 0., 1.])),
        Nucleotide('A', np.array([4., 0., 0.]), np.array([1., 0., 0.]), np.array([0, 0., 1.])),
        Nucleotide('A', np.array([5., 0., 0.]), np.array([1., 0., 0.]), np.array([0, 0., 1.])),
        Nucleotide('A', np.array([6., 0., 0.]), np.array([1., 0., 0.]), np.array([0, 0., 1.])),
    ])
    strand_2 = Strand([
        Nucleotide('T', np.array([1., 2., 0.]), np.array([1., 0., 0.]), np.array([0, 0., 1.])),
        Nucleotide('T', np.array([2., 2., 0.]), np.array([1., 0., 0.]), np.array([0, 0., 1.])),
        Nucleotide('T', np.array([3., 2., 0.]), np.array([1., 0., 0.]), np.array([0, 0., 1.])),
        Nucleotide('T', np.array([4., 2., 0.]), np.array([1., 0., 0.]), np.array([0, 0., 1.])),
        Nucleotide('T', np.array([5., 2., 0.]), np.array([1., 0., 0.]), np.array([0, 0., 1.])),
        Nucleotide('T', np.array([6., 2., 0.]), np.array([1., 0., 0.]), np.array([0, 0., 1.])),
    ])
    system._strands += [strand_1, strand_2]
    print(system)
    print(system.dataframe)
    system.write_oxDNA()
    return

if __name__ == '__main__':
    test_System()
