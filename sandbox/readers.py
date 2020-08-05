from origamiUROP.oxdna import Nucleotide, Strand, System
import numpy as np
import os.path

allowed_bases = ["A", "C", "G", "T", "a", "c", "g", "t"]


class ReaderTopConf:
    """
    Reads in oxDNA .top and .conf files. Converts them to a origamiUROP oxDNA system
    
    CREDIT:
    Copied, reworded and reorganised from:
    --- https://github.com/lorenzo-rovigatti/tacoxDNA/blob/master/src/libs/readers.py 

    """

    def __init__(self, topology_file, configuration_file):
        if not os.path.isfile(configuration_file):
            print("Configuration file '%s' is not readable" % configuration_file)

        if not os.path.isfile(topology_file):
            print("Topology file '%s' is not readable" % topology_file)

        # --- Topology File --- #
        t = open(topology_file, "r")
        self.N, self.N_strands = [int(x) for x in t.readline().split()]
        self._top_lines = t.readlines()

        if len(self._top_lines) != self.N:
            raise Exception(
                "The number of nucleotides specified in the topology file header (%d) is different from the number of nucleotide lines found in the same file (%d)"
                % (self.N, len(self._top_lines))
            )

        # --- Configuration File --- #
        self._conf = None
        c = open(configuration_file, "r")
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
        self._conf = c
        # --- Define System --- #
        self.system = System(box, time, E_pot, E_kin)

    def __del__(self):
        if self._conf:
            self._conf.close()

    def return_system(self):
        strand = False
        current_strand_id = 0
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
                # Add to systemdefine new strand
                if strand:
                    self.system.add_strand(strand)
                strand = Strand()
                current_strand_id = strand_id

            conf_line = self._conf.readline().split()
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
            pos_com = [float(x) for x in conf_line[0:3]]
            a1 = [float(x) for x in conf_line[3:6]]
            a3 = [float(x) for x in conf_line[6:9]]
            v = [float(x) for x in conf_line[9:12]]
            L = [float(x) for x in conf_line[12:15]]

            strand.add_nucleotide(Nucleotide(base, pos_com, a1, a3, v, L))

        self.system.add_strand(strand)

        return self.system


if __name__ == "__main__":
    conf = "./oxdna.out.conf"
    top = "./oxdna.out.top"
    system = ReaderTopConf(top, conf).return_system()
    system.write_oxDNA("reader")
