#!/usr/bin/env python

'''Transform Annotations Results.

This script takes the combined annotations results from Step 04a that 
have had the un-annotated genes added as Hypothetical Genes and 
transforms the combined tsv table in preparation to use with Pandas in
Python for plotting a summary of string matched Gene Types in Step 06.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: January 23rd, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from collections import defaultdict


def gather_ano(file):
    """ processes and filters ano data into a dict. returns dict."""

    # {sequence_name: [U-ID,TrEMBL,T-KO,U-ID,SwissProt,S-KO,KEGG,K-KO]}
    d = defaultdict(list)

    with open(file, 'r') as f:
        _ = f.readline()
        for l in f:
            X = l.rstrip().split('\t')
            gnm = X[0]
            db = X[1]
            uid = X[2]
            if uid == 'n/a': uid = 'No_Match'
            annotation = X[7]
            ko = X[9]
            if ko == 'n/a': ko = 'No_Match'
            if db == 'TrEMBL' or db == 'SwissProt':
                d[gnm].extend([uid, annotation, ko])
            elif db == 'KEGG':
                d[gnm].extend([annotation, ko])

    return d


def combine_results(ano, dbused, out):
    """ preocesses files and outputs combined results to tsv """

    # Parse the input file
    ano_dict = gather_ano(ano)

    # Define header segments for each database.
    header_dict = {
                'TrEMBL': '\tU-ID\tTrEMBL\tT-KO',
                'SwissProt': '\tU-ID\tSwissProt\tS-KO',
                'KEGG': '\tKEGG\tK-KO'
                }
    # Initialize the header string
    header = 'RepSeqName'
    # Build the header string
    for db in dbused:
        header = header + header_dict[db]

    # Write new file
    with open(out, 'w') as o:

        o.write(header + '\n')

        for RepSeqName, annotation in ano_dict.items():
            anos = '\t'.join(annotation)
            newLine = f"{RepSeqName}\t{anos}"
            o.write(newLine + '\n')

    # Function and Script End.

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-a', '--Annotation_file',
        help='The output file from Step 04a through 04b.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-d', '--databases_used',
        help='Please specify which databases you used: KEGG TrEMBL SwissProt',
        metavar='',
        type=str,
        nargs='+',
        required=True
        )
    parser.add_argument(
        '-o', '--out_file',
        help='What do you want to name the output .tsv file?',
        metavar='',
        type=str,
        #required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('Running Script...')
    combine_results(
                        args['Annotation_file'],
                        args['databases_used'],
                        args['out_file']
                        )


if __name__ == "__main__":
    main()
