import argparse
import os
import sys
import numpy as np
import pandas as pd
import pyfaidx as px

sys.path.insert(0, '/home/jmurga/mkt/201903/scripts/src')
from reverseComplement import reverseComplement
from degenerancy import degenerate
import time

# 11589
# args.genes='flybaseGenesCleaned.tab' 
# args.cds='cdsCoordinates.tab' 
# args.population='RAL' 
# args.sampling=160
# args.outgroup='dsim'
# args.path='/home/jmurga/mkt/201903/rawData/dmel'    

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
	parser.add_argument("--singleton", type = bool, default=False, help = "Resampling size")
	parser.add_argument("--seed", type = int, required = False, help = "Input seed")

	# Default arguments
	parser.add_argument("--path", type = str, default = '/home/jmurga/mkt/201903/rawData/dmel', help = "Path to output file")

	# Parsing common arguments
	args = parser.parse_args()

	if(args.outgroup == 'dsim'):
		outgroupFastas = '/data/shared/dgn/outgroup/dsim'
	else:
		outgroupFastas = '/data/shared/dgn/outgroup/dyak'


	dfGenes = pd.read_csv(args.path + '/annotations/'+args.genes,header = 0,sep='\t')
	cds = pd.read_csv(args.path + '/annotations/'+args.cds,header = 0,sep='\t')
	cds = pd.merge(cds, dfGenes,  how='inner', left_on=['chr','name'], right_on = ['chr','name'])
	cds = cds.loc[cds.reset_index().groupby(['chr','id'])['transcriptSize'].idxmax()].reset_index(drop=True)


	# cds = cds.sort_values(['chr','startGene'])

	for index,row in cds.iterrows():
	# for index,row in cds.head(10).iterrows():
	# for index,row in cds[(cds['id']=='FBgn0000054')].iterrows():
		print(index,row['id'])
		start = time.time()

		# Convert CDS list into numeric array
		coordinates = np.array(row['coordinates'].split(',')).astype(int).tolist()
		coordinates =  [coordinates[i:i+2] for i in range(0, len(coordinates), 2)]

		# Open ref and outgroup
		ref = px.Fasta('/data/shared/dgn/ref/Chr' + row['chr'] +'.fasta',sequence_always_upper=True)

		if(args.outgroup == 'dsim'):
			outgroup = px.Fasta(outgroupFastas + '/Chr' + row['chr'] +'_dsim.fasta',sequence_always_upper=True)
		else:
			outgroup = px.Fasta(outgroupFastas + '/Chr' + row['chr'] +'_dyak.fasta',sequence_always_upper=True)


		## Extract ref and outgroup seq
		refSeq = ref.get_spliced_seq(list(ref.keys())[0], coordinates).seq
		outgroupSeq = outgroup.get_spliced_seq(list(outgroup.keys())[0], coordinates).seq
		
		# Check length divisible by 3
		if((len(refSeq) % 3) == 0): 

			# Open multifasta
			multiFasta = px.Fasta('/data/shared/dgn/alignments/'+ args.population + '_Chr' + row['chr'] +'.seq',sequence_always_upper=True)

			# Extract samples from fastas
			samples = list(multiFasta.keys())
						
			# Create empty array with ndimesions equal to multi-Fasta lines and length
			matrix = np.empty([len(samples) + 1, len(refSeq)],dtype='str')
			
			positions=[]
			for i in range(0,len(coordinates),1):
				positions.append(list(range(coordinates[i][0],coordinates[i][1]+1)))  
			positions = [item for sublist in positions for item in sublist]
			positions = np.asarray(positions)

			if(row['strand'] == '-'):
				refSeq = reverseComplement(refSeq)
				degenCode = degenerate(refSeq)
				outgroupSeq = reverseComplement(outgroupSeq)
				positions = positions[::-1]
				# List to append indexes if whole sequence at any population is len(seq) * 'N'			
				# Iter by row matrix to input sequences
				deleteIndex = []
				for i in range(1,len(samples),1):
					tmpFasta = multiFasta.get_spliced_seq(samples[i], coordinates).seq
					tmpFasta = reverseComplement(tmpFasta)

					if(tmpFasta == ('N' * len(tmpFasta))):
						deleteIndex.append(i)
					else:
						matrix[i] = list(tmpFasta)
		
				# Delete lines
				matrix = np.delete(matrix,deleteIndex,0)
				
				# Ref and outgroup
				matrix[0] = list(degenCode)
				matrix[matrix.shape[0]-1] = list(outgroupSeq)

				# NEED TO SOLVE THIS. C-contigous change web subset, need true in order to inter properly witih nditer
				# Create and index with 0 and 4 fold positions in order to include only these positions at matrix and positions variables
				mIndex = (matrix[0]=='0') | (matrix[0]=='4')
				matrix = np.asarray(matrix[:,(matrix[0]=='0') | (matrix[0]=='4')],order='C')
				positions = positions[mIndex]

			elif(row['strand'] == '+'):

				degenCode = degenerate(refSeq)
				# List to append indexes if whole sequence at any population is len(seq) * 'N'			
				deleteIndex = []
				# Iter by row matrix to input sequences
				for i in range(1,len(samples),1):
					tmpFasta = multiFasta.get_spliced_seq(samples[i], coordinates).seq
					if(tmpFasta == ('N' * len(tmpFasta))):
						deleteIndex.append(i)
					else:
						matrix[i] = list(tmpFasta)
		
				# Delete lines
				matrix = np.delete(matrix,deleteIndex,0)
				
				# Ref and outgroup
				matrix[0] = list(degenCode)								
				matrix[matrix.shape[0]-1] = list(outgroupSeq)

				# NEED TO SOLVE THIS. C-contigous change web subset, need true in order to inter properly witih nditer
				mIndex = (matrix[0]=='0') | (matrix[0]=='4')
				matrix = np.asarray(matrix[:,mIndex],order='C')
				positions = positions[mIndex]

			output = list()
			# If not enough lines then save 0 values
			if(matrix.shape[0] < args.sampling):
				# id,pos,div,daf,class
				tmp = [row['id'],row['chr'],0,0,0,'0fold',args.population]
				output.append(tmp)
				# id,pos,div,daf,class
				tmp = [row['id'],row['chr'],0,0,0,'4fold',args.population]
				output.append(tmp)

				# Save in df
				dafDiv = pd.DataFrame(output)
				dafDiv.to_csv(args.path + '/alleleFrequencies/' + args.outgroup + '/' + args.outgroup +'DmelSites'+args.population+'.tab' ,sep='\t',header=False,mode='a',index=False)
				# dafDiv.to_csv(args.path + '/alleleFrequencies/' + args.outgroup + '/' + args.outgroup +'tmp.tab'  ,sep='\t',header=False,mode='a',index=False)

			else:

				# Manual iterator to extract position
				iter = 0
				# Iter each gene matrix to calculate daf and div by class resampling positions
				for x in np.nditer(matrix, order='F',flags=['external_loop']):
					# print(iter)
					# Check if degen is correct. Necesary to save position (manual iter over same length list as matrix)
					if((x[0]!='0') and (x[0]!='4')):
						next
					# If no enough positions then next position
					elif(x[x!='N'][1:-1].shape[0] < args.sampling):
						next
					else:
						degen = x[0]
						AA = x[-1]
						pol = np.random.choice(x[x!='N'][1:-1],args.sampling,replace=False)
						
						# Undefined Ancestra Allele. Try to clean out of the loop
						if(AA == 'N' or AA == '-'):
							next
						# Monomorphic sites. Try to clean out of the loop
						elif(np.unique(pol).shape[0] == 1 and np.unique(pol)[0] == AA):
							next
						else:

							if(degen == '4'):
								functionalClass = '4fold'
							elif(degen == '0'):
								functionalClass = '0fold'
							else:
								degen

							# Check if pol != AA and monomorphic
							if(np.unique(pol).shape[0] == 1 and np.unique(pol)[0] != AA):
								div = 1; AF = 0
								tmp = [row['id'],row['chr'],positions[iter],div,AF,functionalClass,args.population]
								output.append(tmp)
							else:
								AN = pol.shape[0]
								AC = pd.DataFrame(data=np.unique(pol, return_counts=True)[1],index=np.unique(pol, return_counts=True)[0])
								div = 0
								if(AA not in AC.index):
									next
								else:
									# GET DERIVED ALLELE
									AC = AC[AC.index!=AA]
									if(len(AC) == 0):
										next
									# IN ORDER TO NOT CHECK SINGLETONS
									elif((args.singleton==True) and (AC.iloc[0][0]==1)):
										next	
									elif(len(AC) < 2):
										AF = AC.iloc[0]/AN
										AF = AF.iloc[0]
									else:
										next
								tmp = [row['id'],row['chr'],positions[iter],div,AF,functionalClass,args.population]
								output.append(tmp)

					iter+=1
				
				print("--- %s seconds ---" % (time.time() - start))
				# Save results to df
				dafDiv = pd.DataFrame(output)
				# dafDiv.to_csv(args.path + '/alleleFrequencies/' + args.outgroup + '/' + args.outgroup +'tmp.tab'  ,sep='\t',header=False,mode='a',index=False)
				dafDiv.to_csv(args.path + '/alleleFrequencies/' + args.outgroup + '/' + args.outgroup +'DmelSites'+args.population+'.tab' ,sep='\t',header=False,mode='a',index=False)
