#!/usr/bin/env python

'''Best Hit Filter for Tabular Blast Output.

This script filters tabular blast output for best hit based on bitscore,
as well as user defined percent match length, and percent identity.

This tool takes the following input parameters:

    * tabular blast input file that includes query and subject lengths
      in columns 13 and 14 (or 12 & 13 for 0-index)
    * percent_match_length to filter for as float (ex: 50.0 for 50%)
    * percent_identity to filter for as float (ex: 40.0 for 40.0%)

This script returns the following files:

    * input_file.filtered_best_hits.blst

This script requires the following packages:

    * argparse

This file can also be imported as a module and contains the follwing 
functions:

    * filter_blast - This function coordinates the filtering.
    * main - the main function of the script

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Thursday, July 23th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse


def best_hits(query, bitscore, d, line, dups):
    """ Filters the besthit based on bitscore """

    if query in d:
        dups += 1
        old_bitscore = float(d[query].split('\t')[11])

        if bitscore > old_bitscore:
            d[query] = line

    else:
        d[query] = line

    return d, dups


def filter_blast(infile, pml, pid):
    """Gets and prints the spreadsheet's header columns

    Parameters
    ----------
    infile : str
        The file location of the input file
    pml : float
        Percent match length of the query read to filter above. (ex: 0.9)
        alignment length / query length
    rl : int
        Length of the query read to filter above. (ex: 70)

    Returns
    -------
    file
        writes out tabular blast lines passing the filter to 
        a new file infile.filtered_best_hits.blst
    """

    d = {} # initialize dictionary for bitscore besthits
    dups = 0 # counter for number of duplicate matches for reads
    fails = 0 # counter for number of matches failing filters
    passes = 0 # counter for number of matches passing filters
    total = 0 # counter for total blast entries in file

    with open(infile, 'r') as f:
        for l in f:
            total += 1
            X = l.rstrip().split('\t')
            query = X[0] # read identifier
            bitscore = float(X[11]) # bitscore
            pIdent = float(X[2]) # percent sequence identity
            aLen = int(X[3]) # read alignment length
            qLen = int(X[12]) # full length of read
            pMatch = aLen / qLen # percent match length of read length

            if pMatch >= (pml/100) and pIdent >= pid:
                d, dups = best_hits(query, bitscore, d, l, dups)
                passes += 1
            else:
                fails += 1

    outfile = infile.split('.')[0] + '.filtered_best_hits.blst'
    with open(outfile, 'w') as o:
        for k,v in d.items():
            o.write(v)

    print('Total number of entries in blast file:', total)
    print('Number of reads failing the filters:', fails)
    print('Number of reads passing the filters:', passes)
    print('Number of duplicate blast matches passing filter:', dups)
    print('Number of best hit entries written to new file:', len(d))


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--in_file',
        help='Please specify the tabular blast input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-pml', '--percent_match_length',
        help='Percent match length to filter for (ex: 50).',
        metavar='',
        type=float,
        required=True
        )
    parser.add_argument(
        '-pid', '--percent_identity',
        help='Percent Identity to filter for (ex: 40).',
        metavar='',
        type=float,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('Running Script...')
    filter_blast(
        args['in_file'],
        args['percent_match_length'],
        args['percent_identity']
        )


if __name__ == "__main__":
    main()
