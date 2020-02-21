#!/usr/bin/env python

'''Filters KofamScan results for best hits and converts output.

This script takes KofamScan output file and returns a reformatted tab
separated file (to match the other annotation formats of the pipeline)
with the best matches based on the score column.

This tool takes the following input parameters:

    * input file - kofamscan output file

This script returns the following files:

    * output file - reformatted best hit filter results

This script requires the following packages:

    * argparse
    * re

This file can also be imported as a module and contains the follwing 
functions:

    * besthit_filter_kofamscan - kofamscan reformatter / best hit filter
    * main - the main function of the script

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Friday, August, 9th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import re


def besthit_filter_kofamscan(infile, outfile):
    """Filters kofamscan output for best results based on score param.
    and reformats output to a tab separated file.

    Parameters
    ----------
    infile : str
        The file location of the input file
    outfile : str
        Name/location for the output file

    Returns
    -------
    file
        Best hit filtered file
    """

    d = {}

    with open(infile, 'r') as f:
        # Skip first two lines and define new header
        _ = f.readline()
        _ = f.readline()
        new_header = 'Gene_name\tKO\tthrshld\tscore\tE-value\tKO_def\n'

        for l in f:
            # replace spaces with underscore
            rex = re.compile(r'\s+')
            result = rex.sub('\t', l)

            # split line into columns and define variables
            X = result.rstrip().split('\t')
            hq = X[0] 
            nm = X[1]
            KO = X[2]
            thrshld = X[3]
            scr = float(X[4])
            ev = X[5]
            an = ' '.join(X[6:])
            if hq == '*': an = '* ' + an

            # define what the new line looks like
            new_line = f'{nm}\t{KO}\t{thrshld}\t{scr}\t{ev}\t{an}\n'

            # Add new line to dictionary if it has the highest score
            if nm in d:
                old_scr = float(d[nm].split('\t')[3])
                if scr > old_scr:
                    d[nm] = new_line
            else:
                d[nm] = new_line

    # write new file
    with open(outfile, 'w') as o:
        o.write(new_header)
        for k, v in d.items():
            o.write(v)


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--in_file',
        help='Please specify the input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--out_file',
        help='What do you want to name the output file?',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('Running Script...')
    besthit_filter_kofamscan(args['in_file'], args['out_file'])


if __name__ == "__main__":
    main()
