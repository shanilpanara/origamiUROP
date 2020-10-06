from typing import List

import numpy as np
import pandas as pd

from .oxdna import Nucleotide, Strand, System

class Reader:

    columns = [
        'base',
        'x'   , 'y'   , 'z'  ,
        'a1x' , 'a1y' , 'a1z',
        'a3x' , 'a3y' , 'a3z',
        'vx'  , 'vy'  , 'vz' ,
        'Lx'  , 'Ly'  , 'Lz' ,
        'before', 'after', 'strand',
    ]

    """
    Low-level class that takes a table of nucleotides and returns
    a `origamiUROP.oxdna.System` instance. Designed to be subclassed,
    a derived class must pass a dataframe (that can be checked using
    static methods) of nucleotides when initialising, preferably using
    `super()` methods.
    """
    def __init__(
        self, 
        nucleotides: pd.DataFrame, 
        metadata: dict,
    ):

        Reader._check_dataframe(nucleotides)
        Reader._check_metadata(metadata)

        self._nucleotides = nucleotides
        self._metadata = metadata
        
        self._validate_strands()

    @staticmethod
    def _check_metadata(metadata: dict):
        if metadata == None:
            return
        assert True

    @staticmethod
    def _check_dataframe(dataframe: pd.DataFrame):
        assert set(dataframe.columns) == set(Reader.columns)

    def _validate_strands(self):
        assert True

    def _nucleotide_from_series(self, series: pd.Series) -> Nucleotide:
        nuc = Nucleotide(
            series['base'],
            np.array([series['x'], series['y'], series['z']]),
            np.array([series['a1x'], series['a1y'], series['a1z']]),
            np.array([series['a3x'], series['a3y'], series['a1z']]),
            v=np.array([series['vx'], series['vy'], series['vz']]),
            L=np.array([series['Lx'], series['Ly'], series['Lz']]),
        )
        nuc._before = series['before']
        nuc._after = series['after']
        return nuc

    @property
    def strands(self) -> List[Strand]:
        strand_list = []
        strand_counter = True
        i = 0
        while strand_counter:
            i += 1
            temp = self._nucleotides[self._nucleotides['strand']==i]
            strand_counter = temp.any().any()
            if not strand_counter:
                break
            nucleotide_list = []
            for index in temp.index:
                nucleotide_list.append(
                    self._nucleotide_from_series(temp.loc[index])
                )
            # TODO: sort nucleotide_list to ensure 
            # before and after are correct
            strand_list.append(Strand(nucleotides=nucleotide_list))
        print(f"Strand List:\n{strand_list}")
        return strand_list
    
    @property
    def system(self) -> System:
        _system = System(
            self._metadata['box'],
            time=self._metadata.get('time', 0),
            E_pot=self._metadata.get('PE', 0.0),
            E_kin=self._metadata.get('KE', 0.0)
        )
        # reverse list is necessary but cannot remember why
        _system.add_strands(self.strands[::-1])
        return _system

class LAMMPSDataReader(Reader):
    def __init__(self, fname: str):
        super().__init__(self.dataframe, metadata=self.metadata)

    @property
    def dataframe(self) -> pd.DataFrame:
        result = pd.DataFrame()
        return result

    @property
    def metadata(self) -> dict:
        result = {}
        return result

class LAMMPSDumpReader(Reader):
    def __init__(self, fnames: List[str]):
        data, dump = LAMMPSDumpReader.detect_filetypes(fnames)
        super().__init__(self.dataframe, metadata=self.metadata)

    @staticmethod
    def detect_filetypes(fnames: List[str]):
        return data, dump

    @property
    def dataframe(self) -> pd.DataFrame:
        result = pd.DataFrame()
        return result

    @property
    def metadata(self) -> dict:
        result = {}
        return result

class OXDNAReader(Reader):
    def __init__(self, fnames: List[str]):
        conf, top = OXDNAReader.detect_filetypes(fnames)
        self._metadata = {}
        self._configuration = self.read_configuration()
        self._topology = self.read_topology()
        super().__init__(self.dataframe, metadata=self.metadata)

    @staticmethod
    def detect_filetypes(fnames: List[str]):
        conf = fnames[0]
        top = fnames[1]
        return conf, top

    def read_configuration(self, fname: str) -> pd.DataFrame:
        return

    def read_topology(self, fname: str) -> pd.DataFrame:
        return

    @property
    def dataframe(self) -> pd.DataFrame:
        result = pd.DataFrame()
        return result

    @property
    def metadata(self) -> dict:
        result = {}
        return result
