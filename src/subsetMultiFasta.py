import sys
import time
import argparse
import numpy as np
import pandas as pd
import pyfaidx as px

def reverseComplement(data):

	data = data[::-1]
	baseComplement = {'A':'T','C':'G','G':'C','T':'A','N':'N','-':'-'}
	letters = list(data)
	letters = [baseComplement[base] for base in letters]
	data = ''.join(letters)

	return(data)

#################BGD GROUP######################
########UNIVERSITAT AUTONOMA DE BARCELONA#######

if __name__ == '__main__':
	'''Parse arguments and show the required inputs if only name is given to command line'''
	parser = argparse.ArgumentParser(description='Subset coordinates from reference sequence, multi-FASTA file and outgroup sequences. It return a subset multi-FASTA file and txt file with transcript ranges on new multi-FASTA.')
	# Required arguments
	parser.add_argument('--reference', type = str, required = True, help = 'Reference sequence.')
	parser.add_argument('--multiFasta', type = str, required = True, help = 'Raw data to estimate polymorphism.')
	parser.add_argument('--outgroup', type = str, required = True, help = 'Outgroup sequence.')
	parser.add_argument('--coordinates', type = str, required = True, help = 'Coordinates file containing following tabulated columns: start, end, strand, transcript')
	parser.add_argument('--output', type = str, required = True, help = 'Output file without extension.')

	args = parser.parse_args()
	start = time.time()
	

	# START-END COLUMN AND SUMMARIZE BY TRANSCRIPT ALL COLUMNS TO ONE.
	df = pd.read_csv(args.coordinates,sep='\t')
	df['coordinates'] = df['start'].astype(str) + ',' + df['end'].astype(str)
	df['length'] = df.apply(lambda row: row['end']-row['start']+1 ,axis=1)
	
	dfL = df.groupby(['strand','transcript'],sort=False)['length'].sum().reset_index().sort_values('strand').reset_index(drop=True)
	dfC = df.groupby('strand')['coordinates'].apply(lambda x: ','.join(x))
	dfC = dfC.apply(lambda x: [list(map(int,x.split(',')[i:i+2])) for i in range(0, len(x.split(',')), 2)]).reset_index()

	# Open multi-Fasta
	ref = px.Fasta(args.reference, sequence_always_upper=True)
	file = px.Fasta(args.multiFasta, sequence_always_upper=True)
	out = px.Fasta(args.outgroup, sequence_always_upper=True)

	samples = list(file.keys())

	f = open(args.output + '.fa','w')

	# MultiFasta. Iter fasta to add sequence to matrix
	for i in range(0,len(samples),1):
		# Extract each sample sequence
		if(dfC.shape[0] == 2):
			if(i == 0):
				# Ref
				tmp = ref.get_spliced_seq(list(ref.keys())[0], dfC[dfC['strand']=='+']['coordinates'][0]).seq
				tmp = tmp + reverseComplement(ref.get_spliced_seq(list(ref.keys())[0], dfC[dfC['strand']=='-']['coordinates'][1]).seq)
				f.write('>' + list(ref.keys())[0] + '\n' + tmp +'\n')
				# multiFasta
				tmp = file.get_spliced_seq(samples[i], dfC[dfC['strand']=='+']['coordinates'].iloc[0]).seq
				tmp = tmp + reverseComplement(file.get_spliced_seq(samples[i], dfC[dfC['strand']=='-']['coordinates'][1]).seq)
				f.write('>' + samples[i] + '\n' + tmp +'\n')
			elif(i==len(samples)-1):
				# multiFasta
				tmp = file.get_spliced_seq(samples[i], dfC[dfC['strand']=='+']['coordinates'].iloc[0]).seq
				tmp = tmp + reverseComplement(file.get_spliced_seq(samples[i], dfC[dfC['strand']=='-']['coordinates'][1]).seq)
				f.write('>' + samples[i] + '\n' + tmp +'\n')
				# Outgroup
				tmp = out.get_spliced_seq(list(out.keys())[0], dfC[dfC['strand']=='+']['coordinates'].iloc[0]).seq
				tmp = tmp + reverseComplement(out.get_spliced_seq(list(out.keys())[0], dfC[dfC['strand']=='-']['coordinates'].iloc[0]).seq)
				f.write('>' + list(out.keys())[0] + '\n' + tmp +'\n')
			else:
				# multiFasta
				tmp = file.get_spliced_seq(samples[i], dfC[dfC['strand']=='+']['coordinates'].iloc[0]).seq
				tmp = tmp + file.get_spliced_seq(samples[i], dfC[dfC['strand']=='-']['coordinates'].iloc[0]).reverse.complement.seq
				f.write('>' + samples[i] + '\n' + tmp +'\n')
		elif(dfC['strand'][0] == '+'):
			if(i == 0):
				# Ref
				tmp = ref.get_spliced_seq(list(ref.keys())[0], dfC[dfC['strand']=='+']['coordinates'].iloc[0]).seq
				f.write('>' + samples[i] + '\n' + tmp +'\n')
				tmp = file.get_spliced_seq(samples[i], dfC[dfC['strand']=='+']['coordinates'].iloc[0]).seq
				f.write('>' + samples[i] + '\n' + tmp +'\n')
			elif(i==len(samples) -1 ):
				# Ref
				tmp = file.get_spliced_seq(samples[i], dfC[dfC['strand']=='+']['coordinates'].iloc[0]).seq
				f.write('>' + samples[i] + '\n' + tmp +'\n')
				tmp = out.get_spliced_seq(list(out.keys())[0], dfC[dfC['strand']=='+']['coordinates'].iloc[0]).seq
				f.write('>' + list(out.keys())[0] + '\n' + tmp +'\n')
			else:
				tmp = file.get_spliced_seq(samples[i], dfC[dfC['strand']=='+']['coordinates'].iloc[0]).seq
				f.write('>' + samples[i] + '\n' + tmp +'\n')
		elif(dfC['strand'][0] == '-'):
			if(i == 0):
				# Ref
				tmp = ref.get_spliced_seq(list(ref.keys())[0], dfC[dfC['strand']=='-']['coordinates'].iloc[0]).reverse.complement.seq
				f.write('>' + samples[i] + '\n' + tmp +'\n')
				tmp = file.get_spliced_seq(samples[i], dfC[dfC['strand']=='-']['coordinates'].iloc[0]).reverse.complement.seq
				f.write('>' + samples[i] + '\n' + tmp +'\n')
			elif(i == len(samples)):
				tmp = file.get_spliced_seq(samples[i], dfC[dfC['strand']=='-']['coordinates'].iloc[0]).reverse.complement.seq
				f.write('>' + samples[i] + '\n' + tmp +'\n')
				tmp = out.get_spliced_seq(list(out.keys())[0], dfC[dfC['strand']=='-']['coordinates'].iloc[0]).reverse.complement.seq
				f.write('>' + list(out.keys())[0] + '\n' + tmp +'\n')
			else:
				tmp = file.get_spliced_seq(samples[i], dfC[dfC['strand']=='-']['coordinates'].iloc[0]).reverse.complement.seq
				f.write('>' + samples[i] + '\n' + tmp +'\n')

	# Extract positions 
	start = 0
	dfL['seqPos'] = 0
	dfL['startCds'] = 0
	for index,row in dfL.iterrows():
		# print(row['length'])
		# print(str(start)+',' +str(start+row['length']-1))
		dfL.loc[index,['seqPos']] = str(start)+'..' +str(start+row['length']-1)
		dfL.loc[index,['startCds']] = str(start)
		start = start + row['length']

	dfL = dfL[['transcript','strand','startCds','seqPos','length']]
	dfL.to_csv(args.output + 'Positions.txt',header=True,index=False,sep='\t')






	# # Create ndarray with sequences
	# multiFastaMatrixPositive = sequencesToMatrix(file,dfC.iloc[0]['coordinates'],'positive')
	# multiFastaMatrixNegative = sequencesToMatrix(file,dfC.iloc[1]['coordinates'],'negative')

	# multiFastaMatrix = np.hstack([multiFastaMatrixPositive, multiFastaMatrixNegative])


		# matrixAnalysis = multiFastaMatrix[:,start:(start+row['length']-1)]
		# sfs = uSfsFromFasta(matrixAnalysis) 
		# s,d=formatSfs(matrixAnalysis,sfs,'dafFile','divFile','/var/www/html/imkt-dev/files/temporal/')
		# print(s,d)


	# rawSfs = uSfsFromFasta(multiFastaMatrix)
	# # Creating files to R
	# s,d=formatSfs(multiFastaMatrix,rawSfs,'dafFile','divFile','/var/www/html/imkt-dev/files/temporal/')


