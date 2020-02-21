#!/usr/bin/env python

'''Combine Annotations Results.

This script was written for the Salinibacter ruber diversity project to
combine the results from 3 different annotation runs: Using the Uniprot
databases of SwissProt verse TrEMBL and comparing KEGG KO number
annotations with KofamScan. This script collects the results from each
run and writes a new TSV combining all results and keeping track of
which result is from which database.

This tool takes the following input parameters:

    * The besthit filtered tsv results file for each annotation (3)

This script returns the following files:

    * Compare Annotations summary tsv file

This script requires the following packages:

    * argparse
    * os.path
    * collections.defaultdict

This file can also be imported as a module and contains the follwing 
functions:

    * gather_data - appends annotation to list in dict of gene names
    * compare_annotations - the central function of this script
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
import os.path
from collections import defaultdict

def write_file(entry, o):
    """ Writes the output file of combined annotations. """

    for i in entry:
        if i.split('\t')[0] == 'KofamScan':
            X = i.rstrip().split('\t')
            db = X[0]
            nm = X[1]
            KO = X[2]
            thrshld = X[3]
            score = X[4]
            Evalue = X[5]
            KOdef = X[6]
            i = (
                f'{db}\t{nm}\tn/a\tn/a\t{thrshld}\t{score}\t'
                f'{Evalue}\t{KOdef}\tn/a\t{KO}\tn/a\n'
                )
            
        o.write(i)

def triple(name, entry, td, cd, header, basename):
    """ Parses annotations with results from three databases. """
    cd['Triple Annotation'] += 1

    subname = f'{basename}_triples.tsv'

    if not os.path.isfile(subname):
        with open(subname, 'w') as sub_o:
            _ = sub_o.write(header)

    with open(subname, 'a') as sub_o:
        _ = write_file(entry, sub_o)
        sub_o.write('\n')

    return td, cd


def double(name, entry, td, cd, header, basename):
    """ Parses annotations with results from two databases. """
    db1 = entry[0].split('\t')[0]
    db2 = entry[1].split('\t')[0]
    cd[f'{db1}_{db2} Double Annotation'] += 1

    subname = f'{basename}_doubles.tsv'

    if not os.path.isfile(subname):
        with open(subname, 'w') as sub_o:
            _ = sub_o.write(header)

    with open(subname, 'a') as sub_o:
        _ = write_file(entry, sub_o)
        sub_o.write('\n')

    return td, cd


def single(name, entry, td, cd, header, basename):
    """ Parses annotations with result from a single database only. """
    db = entry[0].split('\t')[0]
    cd[f'{db} Single Annotation'] += 1

    subname = f'{basename}_singles.tsv'

    if not os.path.isfile(subname):
        with open(subname, 'w') as sub_o:
            _ = sub_o.write(header)

    with open(subname, 'a') as sub_o:
        _ = write_file(entry, sub_o)
        sub_o.write('\n')

    return td, cd


def gather_data(file, dbname, d):
    """ Collects data from the annotation file and adds it to d."""

    with open(file, 'r') as f:
        _ = f.readline()
        for l in f:
            nm = l.split('\t')[0]
            d[nm].append(f'{dbname}\t{l}')

    return d


def compare_annotations(trembl, sprot, kofamscan, outfile):
    """This is the central function of this script.

    It Coordinates the collection of annotation information from the
    three separate database annotations and compares their results.

    Parameters
    ----------
    trembl : str - *_trembl_cmbnd.blast.best.mtchd
        Best hit filtered and matched annotation results for trembl db
    sprot : str - *_sprot_cmbnd.blast.best.mtchd
        Best hit filtered and matched annotation results for swissprot db
    kofamscan : str - *_kofamscan_results.tsv.best
        Best hit filtered annotation results for kofamscan KO db
    outfile : str
         Name/location for the output file

        Returns
    -------
    tsv file
        tsv file with comparison results from all databases
    """

    # Build dictionary with gene names as keys and values as a list of
    # the annotation line from each database method
    d = defaultdict(list)
    d = gather_data(trembl, 'TrEMBL', d)
    d = gather_data(sprot, 'SwissProt', d)
    d = gather_data(kofamscan, 'KofamScan', d)

    # Read through the dictionary, write out a new combined annotation
    # file, compute the comparisons, and output comparisions results
    single_d = {}
    double_d = {}
    triple_d = {}
    count_d = defaultdict(int)

    with open(outfile, 'w') as o:

        basename = outfile.split('.')[0]

        header = (
            f'DataBase\tGene_Name\tUniprot_ID\tPercent_Match\t'
            f'Alignment_Length/thrshld\tQuery_Length/score\t'
            f'Subject_Length/Evalue\tLong_Gene_Name\tShort_Gene_Name\t'
            f'KO\tInterPro\tPfam\tCOGG(eggNOG)\tGO_F\tGO_CAnnotations\t'
            f'Etc...\n'
            )

        _ = o.write(header)

        for k,v in d.items():

            if len(v) == 1:
                single_d, count_d = single(
                                    k, v, single_d, count_d, header, basename
                                    )

            elif len(v) == 2:
                double_d, count_d = double(
                                    k, v, double_d, count_d, header, basename
                                    )

            elif len(v) == 3:
                triple_d, count_d = triple(
                                    k, v, triple_d, count_d, header, basename
                                    )

            else:
                print(f'Trouble with:\n{v}')

            _ = write_file(v, o)

    for k,v in count_d.items():
        print(f'Total {k}: {v}')


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-trb', '--TrEMBL_best_file',
        help='Specify the trembl_cmbnd.blast.best.mtchd file!',
        metavar='',
        type=str,
        #required=True
        )
    parser.add_argument(
        '-spb', '--SwissProt_best_mtchd_file',
        help='Specify the sprot_cmbnd.blast.best.mtchd file!',
        metavar='',
        type=str,
        #required=True
        )
    parser.add_argument(
        '-kfs', '--KofamScan_best_mtchd_file',
        help='Specify the kofamscan_results.tsv.best file!',
        metavar='',
        type=str,
        #required=True
        )
    parser.add_argument(
        '-o', '--out_file',
        help='What do you want to name the output file?',
        metavar='',
        type=str,
        #required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('Running Script...')
    compare_annotations(
                        args['TrEMBL_best_file'],
                        args['SwissProt_best_mtchd_file'],
                        args['KofamScan_best_mtchd_file'],
                        args['out_file']
                        )


if __name__ == "__main__":
    main()
