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

    global_group.add_argument(
        '--result_name', '--rn',
        default='table.txt',
        help='Name of the output file to store the combined results'
    )

    global_group.add_argument(
        '--mode', '-m',
        default='cmd',
        help='Mode to run the commands. Options are "cmd" for command line and "python" for Python code.'
    )

    args = parser.parse_args()

    combine_results(args.output, args.result_name, args.mode)

if __name__ == '__main__':
    main()