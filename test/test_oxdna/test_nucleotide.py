from origamiUROP.oxdna import Nucleotide
import numpy as np

def test_Nucleotide():
    nucleotide = Nucleotide(
        'A', 
        np.array([1., 0., 0.]),
        np.array([1., 0., 0.]),
        np.array([0., 0., 1.])
    )
    print(nucleotide)
    print(nucleotide.series)

if __name__ == '__main__':
    test_Nucleotide()