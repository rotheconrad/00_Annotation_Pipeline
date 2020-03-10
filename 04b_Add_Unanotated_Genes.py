#!/usr/bin/env python

'''Add Un-Annotated Genes to 04a_Combined or Step 03 annotation files

Not all genes find a match to all databases that pass the filters.
This script adds genes without annotations as "Hypothetical Genes" to
the annotation files so they can be counted in Steps 06 or 07.

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

def Add_Hypotheticals_to_Annotations(annotations, repgenes, dbused, outfile):
    ''' This function adds unannotated genes to the output tsv'''

    d = defaultdict(list)

    with open(annotations, 'r') as f:
            header = f.readline()
            # if using only a single database from Step 03
            if len(dbused) == 1:
                for l in f:
                    X = l.split('\t')
                    gid = X[0]
                    db = dbused[0]
                    d[f'{db}_{gid}'].append(l)
            # if using a combined database from Step 04a
            if len(dbused) > 1:
                for l in f:
                        X = l.split('\t')
                        gid = X[0]
                        db = X[1]
                        d[f'{db}_{gid}'].append(l)

    with open(repgenes, 'r') as f, open(outfile, 'w') as o:
        o.write(header)
        for name, seq in read_fasta(f):
            # blast splits fasta sequence names at first white space
            # Does Kofam scan have same behavior?
            gid = name[1:].split(' ')[0]
            for db in dbused:
                cid = f'{db}_{gid}'
                if cid in d:
                    for i in d[cid]:
                        o.write(i)
                elif cid not in d:
                    # if using only a single database from Step 03
                    if len(dbused) == 1:
                        if dbused[0] == 'KEGG':
                            o.write(
                                f'{gid}\tn/a\tn/a\tn/a\tn/a\t'
                                f'Hypothetical Gene\n'
                                )
                        else:
                            o.write(
                                f'{gid}\tn/a\tn/a\tn/a\t{len(seq)}\t'
                                f'n/a\tHypothetical Gene\tn/a\tn/a\tn/a\tn/a\t'
                                f'n/a\tn/a\tn/a\tn/a\n'
                                )
                    # if using a combined database from Step 04a
                    elif len(dbused) > 1:
                        o.write(
                            f'{gid}\t{db}\tn/a\tn/a\tn/a\t{len(seq)}\t'
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
        '-a', '--annotation_file',
        help='Please specify the Annotation file from Step 03 or 04a!',
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
        '-d', '--databases_used',
        help='Please specify which databases you used: KEGG TrEMBL SwissProt',
        metavar='',
        type=str,
        nargs='+',
        required=True
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
    Add_Hypotheticals_to_Annotations(
                            args['annotation_file'],
                            args['representative_protein_fasta'],
                            args['databases_used'],
                            args['out_file']
                            )

if __name__ == "__main__":
        main()
