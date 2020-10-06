from os import path

import numpy as np
import pandas as pd

from origamiUROP.readers import Reader, OXDNAReader, LAMMPSDataReader, LAMMPSDumpReader
from origamiUROP.oxdna import System

ROOT = "/".join(path.abspath(__file__).split("/")[:-1])

def test_Reader():
    table = pd.DataFrame({
        'base': ['A', 'T', 'C', 'G'],
        'x': [1.0, 2.0, 3.0, 4.0],
        'y': [2.0, 2.0, 2.0, 2.0],
        'z': [1.0, 1.0, 1.0, 1.0],
        'a1x': [1.0, 1.0, 1.0, 1.0],
        'a1y': [0.0, 0.0, 0.0, 0.0],
        'a1z': [0.0, 0.0, 0.0, 0.0],
        'a3x': [0.0, 0.0, 0.0, 0.0],
        'a3y': [1.0, 1.0, 1.0, 1.0],
        'a3z': [0.0, 0.0, 0.0, 0.0],
        'vx': [0.0, 0.0, 0.0, 0.0],
        'vy': [0.0, 0.0, 0.0, 0.0],
        'vz': [0.0, 0.0, 0.0, 0.0],
        'Lx': [0.0, 0.0, 0.0, 0.0],
        'Ly': [0.0, 0.0, 0.0, 0.0],
        'Lz': [0.0, 0.0, 0.0, 0.0],
        'before': [-1, 0, -1, 2],
        'after': [1, -1, 3, -1],
        'strand': [1, 1, 2, 2],
    })
    metadata = {
        'box': np.array([50., 50., 50.])
    }
    reader = Reader(table, metadata)
    system = reader.system

def test_OXDNAReader():
    reader = OXDNAReader([f'{ROOT}/oxdna.test.conf', f'{ROOT}/oxdna.test.top'])
    reader.system.dataframe
    return

def test_LAMMPSDataReader():
    return

def test_LAMMPSDumpReader():
    return

if __name__ == "__main__":
    test_Reader()
    test_OXDNAReader()
    test_LAMMPSDataReader()
    test_LAMMPSDumpReader()
