#!~/.conda/envs/rothon/bin/python

## USAGE :: python scriptname.py uniprot_trembl.dat
## Reads UniProt dat file and returns a tsv of ID Gene_Name DR_Predictions KW_Predictions
## Author :: Roth Conrad :: rotheconrad@gatech.edu :: https://github.com/rotheconrad
## Date Created :: Monday, June 03, 2019
## Date Modified :: Monday, June 17, 2019

def parse_DRs(X, DR):
	''' This function further parses DR annotation entries '''	

	if X[0] == 'KO;':
		DR[0] = X[1].split(';')[0] if (DR[0] == 'n/a') else DR[0] + '|' + X[1].split(';')[0]
	if X[0] == 'InterPro;':
		DR[1] = ' '.join(X[1:]) if (DR[1] == 'n/a') else DR[1] + '|' + ' '.join(X[1:])
	if X[0] == 'Pfam;':
		DR[2] = ' '.join(X[1:]) if (DR[2] == 'n/a') else DR[2] + '|' + ' '.join(X[1:])
	if X[0] == 'eggNOG;':
		DR[3] = ' '.join(X[1:]) if (DR[3] == 'n/a') else DR[3] + '|' + ' '.join(X[1:])
	if X[0] == 'GO;':
		if X[2][0] == 'C':
			DR[5] = ' '.join(X[1:]) if (DR[5] == 'n/a') else DR[5] + '|' + ' '.join(X[1:])
		if X[2][0] == 'F':
			DR[4] = ' '.join(X[1:]) if (DR[4] == 'n/a') else DR[4] + '|' + ' '.join(X[1:])
	else:
		DR[6] = '<' + ' '.join(X[1:]) + '>' if (DR[6] == 'n/a') else DR[6] + ',<' + ' '.join(X[1:]) + '>'

	return DR

def parse_UniProt_datFile(nf, of):
	'''
	This function reads through the UniProt style dat file and rewrites it as a tsv file
	with columns of:
	UniProtID, Full_Gene_Name, Short_Gene_Name, KO, InterPro, Pfam, GOf, GOc, Extended Annotations, NCBI Taxonomy
	'''

	# Simply read through the file, grab what is needed, and write it out.
	with open(nf, 'r') as f, open(of, 'w') as o:

		# Write Header
		header = 'UniProtID\tFull_Gene_Name\tShort_Gene_Name\tKO\tInterPro\tPfam\teggNOG\tGOf\tGOc\tExtended_Annotations\tNCBI Taxonomy\n'
		o.write(header)

		OP = [] # Initialize list to store output
		DE, GN = 0, 0 # Initialize entry counters because we only want the first entry of these
		OH = [] # Initialize list to store taxonomy
		DR = ['n/a','n/a','n/a','n/a','n/a','n/a','n/a'] # Initialize list to store extended annotations

		for l in f:
			X = l.rstrip().split(' ')
			# select the ID
			if X[0] == 'ID':
				OP.append(X[3])
				#print(X[3])
			# select the Full Gene Name
			if X[0] == 'DE' and DE == 0:
				Full_Gene_Name = l.split('=')[1].split('{')[0][:-1]
				OP.append(Full_Gene_Name)
				DE += 1
				#print(Full_Gene_Name)
			# select the Short Gene Name
			if X[0] == 'GN' and GN == 0:
				if DE == 0: OP.append('n/a')
				Short_Gene_Name = l.split('=')[1].split('{')[0][:-1]
				OP.append(Short_Gene_Name)
				GN += 1
				#print(Short_Gene_Name)
			# select OH NCBI Taxonomy
			if X[0] == 'OH':
				tax = ' '.join(X[3:])
				OH.append(tax)
			# select DR annotations
			if X[0] == 'DR':
				DR = parse_DRs(X[3:], DR)
				#print(entry)
			# detect end of entry, write entry out, and reset for next entry
			if X[0] == '//':
				if DE == 0: OP.append('n/a')
				if GN == 0: OP.append('n/a')
				OP.append('\t'.join(DR))
				OP.append('\t'.join(OH))
				o.write('\t'.join(OP) + '\n')
				OP, OH = [], []
				DE, GN = 0, 0
				DR = ['n/a','n/a','n/a','n/a','n/a','n/a','n/a']

	print('parse_UniProt_datFile finished successfully.')

def main():

	# Configure Argument Parser
	import argparse
	parser = argparse.ArgumentParser(description='Reads UniProt dat file (sprot or trembl) and returns a tsv of ID Gene_Name DR_Predictions KW_Predictions')
	parser.add_argument('-i', '--in_file', help='UniProt database dat file required as input!', required=True)
	parser.add_argument('-o', '--out_file', help='What do you want to name the output file?', required=True)
	args=vars(parser.parse_args())

	# Run this scripts main function
	print('Parsing UniProt dat file. If the file is large, which it usually is, this is going to take a while ...')
	parse_UniProt_datFile(args['in_file'], args['out_file'])

if __name__ == "__main__":
	main()
