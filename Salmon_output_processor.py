#!/usr/bin/env python3

"""
Main program to run the Salmon pipeline.

This program will take the output directory of the Salmon pipeline and combine all the quantification results into one table.

@Author: Luis Javier Madrigal-Roca & John K. Kelly

@Date: 2024-12-15

"""

import os
from src.postprocessing_utilities import combine_results

def main():
    parser = argparse.ArgumentParser(description='Salmon postprocessing pipeline. This pipeline allows the combination of the quantification results from the Salmon pipeline into one table.')

    global_group = parser.add_argument_group(
        'Global arguments', 
        description = 'Directory and path settings'
    )

    global_group.add_argument(
        '-o', '--output', 
        required=True, 
        help='Directory with the output of the Salmon pipeline'
    )

    args = parser.parse_args()

    combine_results(args.output)

if __name__ == '__main__':
    main()