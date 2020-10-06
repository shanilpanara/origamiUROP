from os import path

from origamiUROP.readers import Reader
from origamiUROP.oxdna import System

ROOT = "/".join(path.abspath(__file__).split("/")[:-1])

def test_oxDNA_reader():
    return

def test_LAMMPS_reader():
    return

if __name__ == "__main__":
    test_oxDNA_reader()
    test_LAMMPS_reader()
