#!/usr/bin/env python

'''

Build Plots to summarize and compare annotation results from the three
databases in used in the Annotation Pipeline.

Writes stacked barplots as .png files.
Writes gene types counts to .tsv file.
Writes annotations to tsv files for each database.

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
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd

def get_summaries(infile, gene_list, out):
    """ Read input table and generate counts """

    data = defaultdict(list)

    df = pd.read_csv(infile, sep='\t', index_col=0)

    dbs = ['KEGG', 'TrEMBL', 'SwissProt']
    gene_legend = []
    gene_match = []

    # gn is also the row names or index
    # This is a list of key words to look for in the gene annotations
    with open(gene_list) as f:
        header = f.readline()
        for l in f:
            X = l.rstrip().split(', ')
            gene_legend.append(X[0])
            gene_match.append(X[1])
    gene_legend.append('Other Genes')

    gene_counts = open(f'{out}_gene_counts.tsv', 'w')

    Master_other = df
    total = df.shape[0]

    for d in dbs: # for each database

        # the other category is built by subtracting results of each
        # search from the initial databases. Everything remaining after
        # searching all of the gn list is the other category.
        other = df

        print(d)
        gene_counts.write(d)

        # for each gene category in gn list calculation for total
        for g in gene_match:

            # Select gene category
            select = other[d].str.contains(g, case=False, regex=True)
            Mselect = Master_other[d].str.contains(g, case=False, regex=True)

            geneD = other[select]
            # Select remaining genes
            other = other[~select]
            Master_other = Master_other[~Mselect]

            glen = geneD.shape[0]
            data[g].append(glen)

            print(g, glen, f'{glen/total*100:.2f}')
            gene_counts.write(f'{g}\t{glen}\t{glen/total*100:.2f}\n')

        print('\n\n')
        gene_counts.write('\n\n')
        ### Calculate "other" gene category ###
        data['other'].append(float(len(other)))

        ### Write other category for each db to file ###
        other.to_csv(f'{out}_{d}_OtherGenes.tsv', sep='\t')

    ### Write master other category to file ###
    Master_other.to_csv(f'{out}_Master_OtherGenes.tsv', sep='\t')

    ### Convert Dictionaries to DataFrames ####
    df = pd.DataFrame.from_dict(data, orient='index', columns=dbs)

    return df, gene_legend


def Plot_Annotation_Summary(df, gene_legend, legend_columns, out):
    """ Plot Stacked Bar Charts by PanCats for each DB """

    colors = [
        #'#a6cee3','#b2df8a','#fb9a99','#8c510a','#ff7f00','#6a3d9a','#b15928',
        #'#1f78b4','#33a02c','#e31a1c','#fdbf6f','#cab2d6','#ffff99','#01665e'
        '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c', '#8c510a',
        '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', '#01665e'
            ]

    normed_df = df.div(df.sum(axis=0), axis=1)

    fig, ax = plt.subplots(figsize=(20,12))

    #ax1 = df.T.plot.bar(stacked=True, ax=ax1, color=colors)
    ax = normed_df.T.plot.bar(stacked=True, ax=ax, color=colors)

    # set the axis parameters / style
    ax.minorticks_on()
    ax.tick_params(axis='both', labelsize=18)
    ax.tick_params(axis='x', labelrotation=45)
    # set grid style
    ax.yaxis.grid(
        which="minor", color='#f0f0f0', linestyle='--', linewidth=2.5
        )
    ax.yaxis.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=3
        )
    ax.set_axisbelow(True)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles[::-1],
        gene_legend[::-1],
        ncol=legend_columns,
        loc='center left',
        bbox_to_anchor=(0.98, 0.5),
        fancybox=True,
        shadow=True,
        fontsize=18
        )

    # adjust layout, save, and close
    fig.set_tight_layout(True)
    plt.savefig(out)
    plt.close() 


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--transformed_annotation_file',
        help='The Transformed Annotation  file from previous step.',
        metavar='',
        type=str,
        #required=True
        )
    parser.add_argument(
        '-l', '--gene_types_list',
        help='List of gene types to string match.',
        metavar='',
        type=str,
        #required=True
        )
    parser.add_argument(
        '-c', '--legend_columns',
        help='Number of columns for the legend.',
        metavar='',
        type=int,
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
    print('\n\nRunning Script...\n\n')

    df, gene_legend = get_summaries(
                args['transformed_annotation_file'],
                args['gene_types_list'],
                args['out_file']
                )

    Plot_Annotation_Summary(
                        df,
                        gene_legend,
                        args['legend_columns'], 
                        args['out_file']
                        )


if __name__ == "__main__":
    main()
