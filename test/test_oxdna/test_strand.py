from origamiUROP.oxdna.strand import Strand, generateHelix
from origamiUROP.oxdna.nucleotide import Nucleotide
import numpy as np
import pandas as pd


def test_Strand():
    strand = Strand(
        [
            Nucleotide(
                "A",
                np.array([1.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "A",
                np.array([2.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "A",
                np.array([3.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "A",
                np.array([4.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "A",
                np.array([5.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "A",
                np.array([6.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
        ]
    )

    assert strand.sequence == "AAAAAA"
    assert len(strand.nucleotides) == 6
    assert len(strand.dataframe.columns.values) == 9

    strand.sequence = "cccaaa"
    strand_1 = strand.copy()
    strand_1.sequence = "tTtgGg"

    assert strand.sequence == "CCCAAA"
    assert strand_1.sequence == "TTTGGG"
    assert len(strand_1.nucleotides) == 6
    assert len(strand_1.dataframe.columns.values) == 9

    print("Basic Strand: \n", strand)
    # print(strand.dataframe)
    return


def test_generateHelix():
    ssDNA_helix = generateHelix(40, double=False)

    assert type(ssDNA_helix) == list
    assert len(ssDNA_helix) == 1
    assert len(ssDNA_helix[0]) == 40
    assert len(ssDNA_helix[0].nucleotides) == 40
    assert len(ssDNA_helix[0].dataframe.columns.values) == 9

    assert ssDNA_helix[0].index != 0  # required for .top oxDNA file
    assert ssDNA_helix[0].index == 1  # put here for explanatory purposes
    print("ssDNA: \n", ssDNA_helix[0])
    # print(ssDNA_helix[0].dataframe)

    dsDNA_helix = generateHelix(40, double=True)

    assert type(dsDNA_helix) == list
    assert len(dsDNA_helix) == 2
    assert len(dsDNA_helix[0]) == len(dsDNA_helix[1])

    print("dsDNA: \n", dsDNA_helix[0], "\n", dsDNA_helix[1])
    # print(dsDNA_helix[0].dataframe)

    return


if __name__ == "__main__":

    test_Strand()
    test_generateHelix()
