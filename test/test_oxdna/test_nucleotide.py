from origamiUROP.oxdna import Nucleotide
import numpy as np


def test_Nucleotide():
    nucleotide = Nucleotide(
        "A",
        np.array([1.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 0.0, 1.0]),
    )

    assert len(nucleotide.series) == 9
    assert (nucleotide.pos_back == [0.6, 0.0, 0.0]).all()
    assert (nucleotide.pos_back_rel == [-0.4, 0.0, 0.0]).all()
    assert (nucleotide.pos_base == [1.4, 0.0, 0.0]).all()
    assert (nucleotide.pos_stack == [1.34, 0.0, 0.0]).all()
    assert (nucleotide._a2 == [0.0, 1.0, 0.0]).all()

    print(nucleotide)
    print(nucleotide.series)
    return nucleotide


if __name__ == "__main__":
    nucleotide = test_Nucleotide()
