"""
This program will preprocess the raw ouptut from Salmon to be used in voom.

As a result, plots for over and under-expression of genes will be generated, along
other intermediate files useful for upstream analysis.



"""

import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description='Preprocess Salmon output for voom and execute voom-based analysis in R')

if __name__ == '__main__':
    main()