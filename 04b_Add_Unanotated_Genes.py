#!/usr/bin/env python

'''Add UnAnnotated Genes to UniProt_Annotated.tsv files

Not all genes find a match to the UniProt DBs but its good to keep
track of them. This script adds no matches to the UniProt annotations.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: June 4th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from collections import defaultdict

def read_fasta(fp):
    ''' This function parses fasta files and yields name, seq'''
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def Add_NoMatch_to_UniProt_Annotations(ta, qf, outfile):
    ''' This function adds unannotated genes to the output tsv'''

    d = defaultdict(list)

    with open(ta, 'r') as f:
            header = f.readline()
            for l in f:
                    X = l.split('\t')
                    id = X[1]
                    d[id].append(l)

    with open(qf, 'r') as f, open(outfile, 'w') as o:
        o.write(header)
        for name, seq in read_fasta(f):
            id = name[1:].split(' ')[0]
            if id in d:
                for i in d[id]:
                    o.write(i)
            elif id not in d:
                o.write(
                    f'NoMatch\t{id}\tn/a\tn/a\tn/a\t{len(seq)}\t'
                    f'n/a\tHypothetical Gene\tn/a\tn/a\tn/a\tn/a\t'
                    f'n/a\tn/a\tn/a\tn/a\n'
                    )

def main():

    # Configure Argument Parser
    import argparse
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-a', '--uniprot_annotated_tsv',
        help='Please specify the UniProt_Annotated.tsv file!',
        required=True,
        metavar='',
        type=str,
        )
    parser.add_argument(
        '-q', '--representative_protein_fasta',
        help='Please specify the representative protein fasta file!',
        required=True,
        metavar='',
        type=str,
        )
    parser.add_argument(
        '-o', '--out_file',
        help='What would you like to call the new Output file? (.tsv)',
        required=True,
        metavar='',
        type=str,
        )
    args=vars(parser.parse_args())

    # Run this scripts main function
    print('Running Script...')
    Add_NoMatch_to_UniProt_Annotations(
                            args['uniprot_annotated_tsv'],
                            args['representative_protein_fasta'],
                            args['out_file']
                            )

if __name__ == "__main__":
        main()
