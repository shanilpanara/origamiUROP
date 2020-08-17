#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 17:43:05 2020

@author: louiserosset
"""
from argparse import ArgumentParser
import numpy as np

from fiona_gen import generate_system
#from origamiUROP.oxdna import System
#import sys
#sys.path.append('/Users/louiserosset/Documents/Github/origamiUROP/origamiUROP/oxdna')
#from sys import System
#from origamiUROP.oxdna import System
from system import System
sim_nb = 1

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-n', '--number', type=int, nargs=3, required=True)
    parser.add_argument('-ds', '--double-stranded', type=int, nargs=3, required=True)
    parser.add_argument('-ss', '--single-stranded', type=int, nargs=3,required=True)
    parser.add_argument('-f', '--output-prefix')
    args = parser.parse_args()
    system = generate_system(
        args.number[0],
        args.double_stranded[0],
        args.single_stranded[0]
    )



    if args.output_prefix:
        system.write_oxDNA_folder(sim_nb,args.output_prefix)