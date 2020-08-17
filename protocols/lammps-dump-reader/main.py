#!/usr/bin/env python
"""Manual input of LAMMPS dump and data files"""
import numpy as np
import pandas as pd
from origamiUROP.oxdna import System, Strand, Nucleotide
from origamiUROP.oxdna.system import (
    CONFIGURATION_COLUMNS, 
    TOPOLOGY_COLUMNS, 
    LMP_COL_ATOMS, 
    LMP_COL_ELLIPSOIDS, 
    LMP_COL_VELOCITIES
)

class LMPNucleotide(Nucleotide):
    def __init__(
        self, 
        _id,
        _type,
        _pos,
        _mol,
        _v,
        _L,
        _shape,
        _quaternion
    ):
        _lmp_orientation = self._lmp_orientation(_quaternion)
        super().__init__(
            self._lmp_base(_type),
            self._lmp_pos_com(_pos),
            _lmp_orientation['a1'],
            _lmp_orientation['a3'],
            self._lmp_v(_v),
            self._lmp_L(_L),
        )
        self._strand_index = _mol - 1
        self.index = _id - 1

    @staticmethod
    def _lmp_orientation(quaternion) -> dict:
        return {
            'a1' : np.array([])
            'a2' : np.array([])
            'a3' : np.array([]),
        }

    @staticmethod
    def _lmp_base(lmp_type) -> str:
        return

    @staticmethod
    def _lmp_pos(lmp_position) -> np.ndarray:
        return

    @staticmethod
    def _lmp_v(lmp_velocity) -> np.ndarray:
        return

    @staticmethod
    def _lmp_L(lmp_angular) -> np.ndarray:
        return


def import_LAMMPS_data(fname: str) -> pd.DataFrame:
    n_atoms = int()
    n_bonds = int()

    box = np.zeros((3, 2))

    # import atom table
    atoms = pd.read_csv(fname, nrows=n_atoms)

    # import velocities
    velocities = pd.read_csv(fname, nrows=n_atoms)

    # import ellipsoid table
    ellipsoids = pd.read_csv(fname, nrows=n_atoms)

    # import bonds
    bonds = pd.read_csv(fname, nrows=n_bonds)

    output = pd.DataFrame()
    assert set(output.columns) == set(LMP_COL_ATOMS + LMP_COL_ELLIPSOIDS + LMP_COL_VELOCITIES)
    return output

def detect_filetype(fname: str) -> str:
    filetype = ""
    return filetype

def main(fname: str, kind: str = None):
    if kind == None:
        kind = detect_filetype(fname)
    return

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("fname", type=str)
    parser.add_argument("")
    main(**vars(parser.parse_args()))