#!/usr/bin/env python

'''

Build Plots to compare annotation results for different samples.

A Sample can be proteins predicted from a MAG, genome, or metagenome.

writes stacked barplots as .png files.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: March 10th, 2020
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

def get_summaries(file_list, gene_list, dbused, out):
    """ Read input files and count str.match of gene type in gene list"""

    # depending on the database used, set which column has the gene annotations
    if dbused == 'KEGG':
        anno_col = 'KO_def'
        num_col = [5]
    if dbused in ['SwissProt', 'TrEMBL']:
        anno_col = 'Long_Gene_Name'
        num_col = [6]

    # initialize list objects to store gene_list values
    gene_legend = []
    gene_match = []
    legend_color = []

    # Read in the gene list
    with open(gene_list) as f:
        header = f.readline()
        for l in f:
            X = l.rstrip().split(', ')
            gene_legend.append(X[0])
            gene_match.append(X[1])
            if len(X) == 3: legend_color.append(X[2])
    gene_legend.append('Other Genes')

    # initialize dictionary to store count values for gene list
    data = defaultdict(list)
    # initialize list of file names to use as column names
    colnames = []

    for file in file_list: # for each file
        # Print progress
        print(f'\n\nParsing file: {file}')
        # Print header to screen
        print(f'\nGeneType\tCount\tPercent')
        # get file basename to use as column name
        colname = file.split('/')[-1].split('.')[0]
        colnames.append(colname)
        # read in each file as a tab separated dataframe
        df = pd.read_csv(file, sep='\t', usecols=num_col)
        # Count total genes in file
        total = df.shape[0]
        # the other category is built by subtracting results of each
        # search from the initial databases. Everything remaining after
        # searching all of the gn list is the other category.
        other = df
        # count total genes matching each gene category in gene list
        for g in gene_match:
            # String match each g and select matching rows
            select = other[anno_col].str.contains(g, case=False, regex=True)
            # set geneD equal to matches
            geneD = other[select]
            # Select all remaining gene that do not match for next round.
            other = other[~select]
            # Count number of genes that match g
            lenD = geneD.shape[0]
            # add count to data dict for each file for each g
            data[g].append(lenD)
            # values to screen
            print(
                g, lenD, f'{lenD/total*100:.2f}',
                )
        ### Calculate "other" gene category ###
        data['other'].append(float(len(other)))
        # values to screen
        print(
            'Other Genes', lenD, f'{lenD/total*100:.2f}',
            )

    ### Convert Dictionaries to DataFrames ####
    df = pd.DataFrame.from_dict(data, orient='index', columns=colnames)

    return df, gene_legend, legend_color


def Plot_Annotation_Summary(df, gene_legend, legend_color, legend_columns, out):
    """ Plot Stacked Bar Charts by PanCats for each DB """

    print('\n\nBuilding the plot...')

    # Select colors to use. User defined or default.
    if len(legend_color) > 0:
        colors = legend_color
        # append color for other category
        colors.append('#d9d9d9')
    else:
        colors = [
            '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
            '#8c510a','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99',
            '#b15928','#01665e'
                ]

    normed_df = df.div(df.sum(axis=0), axis=1)

    fig, ax = plt.subplots(figsize=(20,12))

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
        '-i', '--input_file_list',
        help='List of annotation files to plot - each file gets a stacked bar.',
        metavar='',
        type=str,
        nargs='+',
        required=True
        )
    parser.add_argument(
        '-l', '--gene_types_list',
        help='List of gene types to string match.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-c', '--legend_columns',
        help='(Optional) Number of columns for the legend (Default=1).',
        metavar='',
        type=int,
        default=1
        )
    parser.add_argument(
        '-d', '--database_used',
        help='Please specify which databases you used: KEGG TrEMBL SwissProt',
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
    print('\n\nRunning Script...')

    # Read in the files and build a dataframe df of gene types counts.
    df, gene_legend, legend_color = get_summaries(
                                        args['input_file_list'],
                                        args['gene_types_list'],
                                        args['database_used'],
                                        args['out_file']
                                        )

    # Build the plot
    Plot_Annotation_Summary(
                        df,
                        gene_legend,
                        legend_color,
                        args['legend_columns'],
                        args['out_file']
                        )

    print('\nLooks like the script finished successfully!\n\n')

if __name__ == "__main__":
    main()
