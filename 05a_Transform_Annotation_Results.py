#!/usr/bin/env python

'''Transform Annotations Results.

This script takes the combined annotations results including the NoMatch
from the 3 databases TrEMBL, SwissProt, and KEGG and transform the tsv
table in preparation for use with Python Pandas and plotting.

This tool takes the following input files:

    * 03_*_ClstrRepSeq_Annotations_NoMatch.tsv

This script returns the following files:

    * tsv file with transformed tsv table.

This script requires the following packages:

    * argparse 
    * collections.defaultdict

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
            db = X[0]
            gnm = X[1]
            uid = X[2]
            if uid == 'n/a': uid = 'No_Match'
            annotation = X[7]
            ko = X[9]
            if ko == 'n/a': ko = 'No_Match'
            if db == 'TrEMBL' or db == 'SwissProt':
                d[gnm].extend([uid, annotation, ko])
            elif db == 'KofamScan':
                d[gnm].extend([annotation, ko])

    return d


def combine_results(ano, out):
    """ preocesses files and outputs combined results to tsv """

    # Parse the input file
    # {RepSeqName: [U-ID,TrEMBL,T-KO,U-ID,SwissProt,S-KO,KEGG,K-KO]}
    ano_dict = gather_ano(ano)

    # Set the output header
    header = (
        'RepSeqName\tU-ID\tTrEMBL\tT-KO\tKEGG\tK-KO\tU-ID\tSwissProt\tS-KO'
         )

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
        help='The 03_*_ClstrRepSeq_Annotations_NoMatch.tsv file.',
        metavar='',
        type=str,
        #required=True
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
                        args['out_file']
                        )


if __name__ == "__main__":
    main()
