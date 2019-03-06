import argparse
import os
import sys
import numpy as np
import pandas as pd
from numpy import array 
from pyfaidx import Fasta
sys.path.insert(0, '/home/jmurga/mkt/201902/scripts/src')
from reverseComplement import reverseComplement
from degenerate import degenerate
from dafResampling import dafWithResampling
import time

if __name__ == "__main__":
	'''Parse arguments and show the required inputs if only name is given to command line'''
	parser = argparse.ArgumentParser(description='Extract alleles frequencies from multi-FASTA aligment')
	# Required arguments
	parser.add_argument("--genes", type = str, required = True, help = "File basic gene information.")
	parser.add_argument("--cds", type = str, required = True, help = "File all transcript coordinates by genes.")
	parser.add_argument("--outgroup", type = str, required = True, help = "Select outgroup to compute diverenge and derived allele frequency",choices=['dsim','dyak'])
	# Optional arguments
	parser.add_argument("--population", type = str, required = True, help = "Select population to extract")
	parser.add_argument("--sampling", type = int, required = True, help = "Resampling size")
	parser.add_argument("--seed", type = int, required = False, help = "Input seed")

	# Default arguments
	parser.add_argument("--path", type = str, default = '/home/jmurga/mkt/201902/rawData', help = "Path to output file")

	# Parsing common arguments
	args = parser.parse_args()

	if(args.outgroup == 'dsim'):
		outgroupFastas = args.path + '/fastas/outgroup/dsim'
		outputHeader = 'Dsimulans'
	else:
		outgroupFastas = args.path + '/fastas/outgroup/dyak'
		outputHeader = 'Dyakuba'


	dfGenes = pd.read_csv(args.path + '/annotations/'+args.genes,header = 0,sep='\t')
	cds = pd.read_csv(args.path + '/annotations/'+args.cds,header = 0,sep='\t')
	cds = pd.merge(cds, dfGenes,  how='inner', left_on=['chr','name'], right_on = ['chr','name'])
	cds = cds.loc[cds.reset_index().groupby(['chr','id'])['transcriptSize'].idxmax()].reset_index(drop=True)
	# cds = cds.sort_values(['chr','startGene'])

	for index,row in cds.iterrows():
		print(index,row['id'])
		start = time.time()

		# Convert CDS list into numeric array
		coordinates = array(row['coordinates'].split(',')).astype(int).tolist()
		coordinates =  [coordinates[i:i+2] for i in range(0, len(coordinates), 2)]

		# Open ref and outgroup
		ref = Fasta(args.path + '/fastas/ref/Chr' + row['chr'] +'.fasta')
		outgroup = Fasta(outgroupFastas + '/Chr' + row['chr'] +'_dsim.fasta')

		## Extract ref and outgroup seq
		refSeq = ref.get_spliced_seq(row['chr'], coordinates).seq.upper()
		outgroupSeq = outgroup.get_spliced_seq(outputHeader, coordinates).seq.upper()
		
		if('M' in refSeq):
			continue
		else:
			if((len(refSeq)/3).is_integer()): 
				
				# Open population multifasta
				popFasta = Fasta(args.path + '/fastas/alignments/' + args.population + '_Chr' + row['chr'] +'.seq')
				
				#Extract samples
				samples = list(popFasta.keys())
				
				# All positions as list to associate nt to pos in df
				positions=[]
				for i in range(0,len(coordinates),1):
					positions.append(list(range(coordinates[i][0],coordinates[i][1]+1)))  
				allPositions = [str(item) for sublist in positions for item in sublist]

				## Empty matrix to append ref, pop
				matrix0f=np.empty([len(samples)+2,len(refSeq)],dtype='str')
				matrix4f=np.empty([len(samples)+2,len(refSeq)],dtype='str')

				# Outgroup sequences to pd.DataFrame
				outgroup = np.empty([1,len(refSeq)],dtype='str')

				# Open population fastas and introduce to matrix depending on strand
				if(row['strand'] == '-'):
					refSeq = reverseComplement(refSeq)
					refSeq0f,refSeq4f = degenerate(refSeq)
					outgroupSeq = reverseComplement(outgroupSeq)
					allPositions = allPositions[::-1]

					matrix4f[0] = list(refSeq4f)        
					matrix0f[0] = list(refSeq0f)
					# Iter by row matrix to input sequences
					deleteIndex = []
					for i in range(0,len(samples),1):
						tmp = popFasta.get_spliced_seq(samples[i], coordinates).seq.upper()
						tmp = reverseComplement(tmp)
						if(tmp == ('N' * len(tmp))):
							deleteIndex.append(i+1)
						else:
							matrix4f[i+1] = list(tmp)
							matrix0f[i+1] = list(tmp)
					matrix4f = np.delete(matrix4f,deleteIndex,0)
					matrix0f = np.delete(matrix0f,deleteIndex,0)
				else:
					refSeq0f,refSeq4f = degenerate(refSeq)
					matrix4f[0] = list(refSeq4f)        
					matrix0f[0] = list(refSeq0f)
					# Iter by row matrix to input sequences
					deleteIndex = []
					for i in range(0,len(samples),1):
						tmp = popFasta.get_spliced_seq(samples[i], coordinates).seq.upper()
						if(tmp == ('N' * len(tmp))):
							deleteIndex.append(i+1)
						else:
							matrix4f[i+1] = list(tmp)
							matrix0f[i+1] = list(tmp)

					matrix4f = np.delete(matrix4f,deleteIndex,0)
					matrix0f = np.delete(matrix0f,deleteIndex,0)
					# Iter by matrix column to degenerate all sequences based on reference sequence

					
				# OutputSeq to last item in matrix
				matrix4f[len(matrix4f)-1] = list(outgroupSeq)
				matrix0f[len(matrix4f)-1] = list(outgroupSeq)

				# 4fold matrix positions to pd.DataFrame
				df4f = pd.DataFrame(matrix4f)
				# Columns as positions
				df4f.columns = allPositions
				df4f = df4f.loc[:,df4f.iloc[0]!='N']   
				# Transpose df to iter by row
				df4f = df4f.transpose()
			
				# 0fold matrix positions to pd.DataFrame
				df0f = pd.DataFrame(matrix0f)
				# Columns as positions
				df0f.columns = allPositions
				df0f = df0f.loc[:,df0f.iloc[0]!='N'] 
				# Transpose df to iter by row
				df0f = df0f.transpose()

				# Extracting daf and div resampling by position
				dmelFourFoldDafDiv = dafWithResampling(id=row['id'],data=df4f,resamplingValue=args.sampling,type='4fold')
				dmelZeroFoldDafDiv = dafWithResampling(id=row['id'],data=df0f,resamplingValue=args.sampling,type='0fold')
				dmelFourFoldDafDiv['chr'] = row['chr']
				dmelZeroFoldDafDiv['chr'] = row['chr']

				# Save results to df
				dmelFourFoldDafDiv.to_csv(args.path + '/alleleFrequencies/' + args.population + 'FourFold.tab',sep='\t',header=False,mode='a',index=False)
				dmelZeroFoldDafDiv.to_csv(args.path + '/alleleFrequencies/' + args.population + 'ZeroFold.tab',sep='\t',header=False,mode='a',index=False)
				
			else:
				dmelFourFoldDafDiv = pd.DataFrame({'id':row['id'],'POS':0,'rawDerivedAllele':0,'div':0,'type':'4fold'},index=[0])
				dmelZeroFoldDafDiv = pd.DataFrame({'id':row['id'],'POS':0,'rawDerivedAllele':0,'div':0,'type':'0fold'},index=[0])
				dmelFourFoldDafDiv['chr'] = row['chr']
				dmelZeroFoldDafDiv['chr'] = row['chr']

				dmelFourFoldDafDiv.to_csv(args.path + '/alleleFrequencies/' + args.population + 'FourFold.tab',sep='\t',header=False,mode='a',index=False)
				dmelZeroFoldDafDiv.to_csv(args.path + '/alleleFrequencies/' + args.population + 'ZeroFold.tab',sep='\t',header=False,mode='a',index=False)

		print("--- %s seconds ---" % (time.time() - start))
