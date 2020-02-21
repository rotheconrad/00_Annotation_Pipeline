#!~/.conda/envs/rothon/bin/python

## USAGE :: python scriptname.py -i infile.name -o outfile.name
## This script is to convert this file: https://www.genome.jp/kegg-bin/get_htext?ko00001.keg
## ( Select "Download htext" from the top of the page )
## into a DataFrame style tab separated file containing A B C D on each line.
## Author :: Roth Conrad :: rotheconrad@gatech.edu :: https://github.com/rotheconrad
## Date Created :: Tuesday, June 18, 2019
## Date Updated :: N/A first version

def parse_kegg_hierarchy(i, o):
	'''
	Takes Hierarchical file and prints TSV with A B C D on each line.
	'''

	with open(i, 'r') as f, open(o, 'w') as q:

		header = 'PATHWAY\tMODULE\tPATH\tK-number\tShort_Gene_Name\tLong_Gene_Name\tEC-number\n'
		q.write(header)

		A, B, C, D = '', '', '', ''

		for l in f:
			X = l.rstrip().split('  ')
			
			if l[0] == 'A': A = X[0][1:]
			if (X[0] == 'B') & (len(X) > 1): B = X[1]
			if X[0] == 'C': C = X[2]
			if X[0] == 'D':
				x = X[4].split(';')
				z = x[1].split('[')
				D = f"{X[3]}\t{x[0]}\t{z[0]}\t{z[1].split(']')[0]}" if (len(z) > 1) else f"{X[3]}\t{x[0]}\t{z[0]}\tn/a"
				q.write(f'{A}\t{B}\t{C}\t{D}\n')


def main():

	# Configure Argument Parser
	parser = argparse.ArgumentParser(description='This script is to convert this file: https://www.genome.jp/kegg-bin/get_htext?ko00001.keg \
												  into a DataFrame style tab separated file containing A B C D on each line.')
	parser.add_argument('-i', '--in_file', help='Please specify the ko00001.keg input file!', required=True)
	parser.add_argument('-o', '--out_file', help='What do you want to name the output TSV file?', required=True)
	args=vars(parser.parse_args())

	# Run this scripts main function
	print('Running Script...')
	parse_kegg_hierarchy(args['in_file'], args['out_file'])

if __name__ == "__main__":
	import argparse
	main()
