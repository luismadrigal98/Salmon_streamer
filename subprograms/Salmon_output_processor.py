#!/usr/bin/env python3

"""
Main program to run the Salmon pipeline.

This program will take the output directory of the Salmon pipeline and combine all the quantification results into one table.

@Author: Luis Javier Madrigal-Roca & John K. Kelly

@Date: 2024-12-15

"""

import sys
import os
import argparse

# Add the parent directory to sys.path to allow imports from sibling directories
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.postprocessing_utilities import combine_results

def main(args):

    combine_results(args.output, args.result_name, args.mode, args.includes_alternative_genome)

if __name__ == '__main__':
    main()