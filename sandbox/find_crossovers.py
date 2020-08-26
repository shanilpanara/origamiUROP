import itertools
from origamiUROP.lattice import Lattice

def find_crossover_locations(
    max_size: int = 30,
    bp_per_turn: float = 10.45,
    kind: str = "half",
    origin: int = 0
) -> dict:
    """
    Function which returns an array containing the best locations 
    for half or whole turns of the DNA strand where the strand will
    stay in plane and have the least strain/stress build-up.
    
    This is computed given the number of basepairs per turn
    & locations are given in no. of base pairs

    Arguments:
    bp_per_turn - (default: 10.45)
    kind - either half turns or whole turns
    max_size - max no. of turns to base pairs to the list
    origin - what number crossovers start at, either 0 (python indexing) or 1 
    """

    no_of_bp = 0
    if origin == 0: # i.e. Python indexing
        crossover_location = [0]
        shift = 1
    else:
        crossover_location = [1]
        shift = 0


    # half turn
    if kind == "half":
        i = 1  # every odd number
        while no_of_bp < max_size:
            # minus 1 at the end because Python indexing starts at 0
            no_of_bp = int(round(bp_per_turn * 0.5 * i, 0) - shift)
            crossover_location.append(no_of_bp)
            i += 2

    # whole turns
    elif kind == "whole":
        i = 2  # every even number
        while no_of_bp < max_size:
            # minus 1 at the end because Python indexing starts at 0
            no_of_bp = int(round(bp_per_turn * 0.5 * i, 0) - shift)
            crossover_location.append(no_of_bp)
            i += 2
    return crossover_location


def calculate_row_size(*half_turn_sizes: int):
    """
    Function to calculate the row size given a number of integer arguments.
    Returns: (sum - 1) and (list of half turn sizes)

    Arguments:
    half_turn_sizes - integers representing no. of nucleotides in a half turn
    """
    row_size = 0
    for half_turn_size in half_turn_sizes:
        row_size += half_turn_size
    
    # print("Half Turns at:",list(half_turn_sizes))
    row_size -= 1
    return row_size, list(half_turn_sizes)

def get_row_sizes(max_row_size, columns, bp_per_turn: float = 10.45):
    """
    Returns list and dictionary of row sizes:
        list - all possible sizes of rows as a sum of row sizes for the given number of columns
        dictionary - keys formed of values in list & values are half turn sizes which sum up to value of key

    Arguments:
        max_row_size - maximum number of lattice sites in a single row
        columns - number of pairs of crossovers on each row

    Note all key's are given as Python index values, i.e. they start from 0

    """
    half_turn_sizes = find_crossover_locations(max_row_size, bp_per_turn, origin = 1, )[1:]
    crossover_combinations = itertools.combinations_with_replacement(half_turn_sizes, columns)
    row_sizes_and_turns = list(itertools.starmap(calculate_row_size, crossover_combinations))
    
    # Generate dictionary data structure
    rowsize_dict = dict()
    for item in row_sizes_and_turns:
        key = item[0]
        value = item[1]
        rowsize_dict.setdefault(key,[]).append(value)
    
    # Generate list of only row sizes + remove duplicates + order
    rowsize_list = sorted(set([item[0] for item in row_sizes_and_turns]))
    
    return rowsize_list, rowsize_dict

rowsize_list, rowsize_dict = get_row_sizes(100, 2)
print(rowsize_list)
print(rowsize_dict)

class DNAScaffold(Lattice):
    pass