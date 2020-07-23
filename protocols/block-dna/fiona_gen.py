from origamiUROP.oxdna import System, Strand, Nucleotide
from argparse import ArgumentParser

def generate_system(
        box : float,
        n : int,
        fraction_ds : float
    ) -> System:
    """
    Generates an oxDNA system containing a single piece of DNA
    which has blocks of equal size of double-stranded portions
    and single-stranded portions.

    Parameters:
        n - number of nucleotides e.g. 100
        fraction_ds - fraction of double-strand DNA e.g. 0.5
    
    Returns:
        system - oxDNA system
    """
    system = System()
    return system

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-b', '--box', type=float)
    parser.add_argument('-n', '--number', type=int)
    parser.add_argument('-ds', '--double-stranded', type=float)