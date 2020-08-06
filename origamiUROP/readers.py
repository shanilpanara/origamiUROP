from origamiUROP.oxdna import Nucleotide, Strand, System
import numpy as np
from os import path

allowed_bases = ["A", "C", "G", "T", "a", "c", "g", "t"]


class ReaderTopConf:
    """
    Reads in oxDNA .top and .conf files. Converts them to a origamiUROP oxDNA system
    
    CREDIT:
    Copied, edited and reorganised from:
    --- https://github.com/lorenzo-rovigatti/tacoxDNA/blob/master/src/libs/readers.py 

    """

    def __init__(self, topology_file, configuration_file):
        self._conf = False
        self.conf_file = configuration_file
        self.top_file = topology_file

        if not path.isfile(configuration_file):
            print("Configuration file '%s' is not readable" % configuration_file)

        if not path.isfile(topology_file):
            print("Topology file '%s' is not readable" % topology_file)

    def return_system(self):
        # --- Topology File --- #
        t = open(self.top_file, "r")
        self.N, self.N_strands = [int(x) for x in t.readline().split()]
        self._top_lines = t.readlines()

        if len(self._top_lines) != self.N:
            raise Exception(
                "The number of nucleotides specified in the topology file header (%d) is different from the number of nucleotide lines found in the same file (%d)"
                % (self.N, len(self._top_lines))
            )
        # --- Configuration File --- #
        c = open(self.conf_file, "r")
        try:
            timeline = c.readline()
            time = float(timeline.split()[2])
            box = np.array([float(x) for x in c.readline().split()[2:]])
            E_tot, E_pot, E_kin = [float(x) for x in c.readline().split()[2:]]

        except Exception as e:
            raise Exception(
                "The header lines of the configuration file are invalid (caught a '%s' exception)"
                % e
            )

        # --- Define System --- #
        system = System(box, time, E_pot, E_kin)

        # --- Main reading loop --- #
        strand = False
        current_strand_id = 0
        nuc_list = []
        for line_index, top_line in enumerate(self._top_lines):
            top_line_split = top_line.split()
            strand_id = int(top_line_split[0])  # strand index
            nuc_3 = int(top_line_split[2])
            nuc_5 = int(top_line_split[3])

            if top_line_split[1] in allowed_bases:
                base = top_line_split[1]
            else:
                raise ValueError(
                    "The line n. %d in the topology file contains an incorrect base"
                    % line_index
                )

            if strand_id != current_strand_id:
                # First time around don't do this, wait for nuc_list to fill

                if len(nuc_list) != 0:
                    strand = Strand(nuc_list)
                    system.add_strand(strand)
                    nuc_list = []

                current_strand_id = strand_id

            conf_line = c.readline().split()

            if len(conf_line) == 0:
                raise Exception(
                    "The %d-th nucleotide line in the configuration file is empty"
                    % line_index
                )
            elif len(conf_line) != 15:
                raise Exception(
                    "The %d-th nucleotide line in the configuration file is invalid"
                    % line_index
                )

            pos_com = np.array([float(x) for x in conf_line[0:3]])
            a1 = np.array([float(x) for x in conf_line[3:6]])
            a3 = np.array([float(x) for x in conf_line[6:9]])
            v = np.array([float(x) for x in conf_line[9:12]])
            L = np.array([float(x) for x in conf_line[12:15]])

            nuc_list.append(Nucleotide(base, pos_com, a1, a3, v, L))

        strand = Strand(nuc_list)
        system.add_strand(strand)

        c.close()
        t.close()
        return system
