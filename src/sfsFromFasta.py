import sys
import time
import numpy as np
import pandas as pd
import pyfaidx as px
def degenerancy(data,codonDict):
	
	#DEGENERANCY DICTIONARIES
	standardDict = {
		'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202',
		'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004',
		'TAT': '002', 'TAC': '002', 'TAA': '022', 'TAG': '002',
		'TGT': '002', 'TGC': '002', 'TGA': '020', 'TGG': '000',
		'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204',
		'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004',
		'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002',
		'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204',
		'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000',
		'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004',
		'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002',
		'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202',
		'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004',
		'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004',
		'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002',
		'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'}

	if(codonDict == 'standard'):
		degenerateCodonTable = standardDict

	degenerancy = ''
	for i in range(0, len(data),3):
		codon = data[i:i+3]
		if('N' in codon or '-' in codon):
			degenerancy += codon
		else:
			degenerancy += degenerateCodonTable[codon]

	return(degenerancy)
def sequencesToMatrix(multiFasta,split=None):

	# Extract samples from fastas
	samples = list(multiFasta.keys())
	
	if(split is None):
		seqLen = len(multiFasta[samples[0]][:].seq)
		
		if((seqLen % 3) != 0):
			print('cdsLength')
			sys.exit('cdsLength')

		# Create empty array with ndimesions equal to multi-Fasta lines and length
		matrix = np.empty([len(samples),len(multiFasta[samples[0]][:].seq)],dtype='str')
		
		# List to append indexes if whole sequence at any population is len(seq) * 'N'
		deleteIndex = list()
		
		# Iter fasta to add sequence to matrix
		for i in range(1,len(samples),1):
			# Extract each sample sequence
			tmp = multiFasta[samples[i]][:].seq
			if(len(tmp) != seqLen):
				print('errorAlign')
				sys.exit('errorAlign')
			if('N' in tmp):
				deleteIndex.append(i)
			else:
				matrix[i] = list(tmp)
		
		# Delete lines
		matrix = np.delete(matrix,deleteIndex,0)
		
		degenCode = degenerancy(multiFasta[samples[0]][:].seq,codonTable)
		# Put degenerancy in first ndarray element
		matrix[0] = list(degenCode)
		# NEED TO SOLVE THIS. C-contigous change web subset, need true in order to inter properly witih nditer

		matrix = np.asarray(matrix[:,(matrix[0]=='0') | (matrix[0]=='4')],order='C')

		return(matrix)

	else:
		# Create empty array with ndimesions equal to multi-Fasta lines and length
		split = [split[0]+1,split[1]+1]
		matrix = np.empty([len(samples),(split[1]-split[0]+1)],dtype='str')
		
		# List to append indexes if whole sequence at any population is len(seq) * 'N'
		deleteIndex = list()
		
		# Iter fasta to add sequence to matrix
		for i in range(1,len(samples),1):
			# Extract each sample sequence
			tmp = multiFasta.get_spliced_seq(samples[i], [split]).seq
			if(tmp == ('N' * len(tmp))):
				deleteIndex.append(i)
			else:
				matrix[i] = list(tmp)
		
		# Delete lines
		matrix = np.delete(matrix,deleteIndex,0)
		
		degenCode = degenerancy( multiFasta.get_spliced_seq(samples[i], [split]).seq.upper(),codonTable)
		# Put degenerancy in first ndarray element
		matrix[0] = list(degenCode)
		# NEED TO SOLVE THIS. C-contigous change web subset, need true in order to inter properly witih nditer

		matrix = np.asarray(matrix[:,(matrix[0]=='0') | (matrix[0]=='4')],order='C')  

		return(matrix)
def uSfsFromFasta(sequenceMatrix):
	
	output = list()


	for x in np.nditer(sequenceMatrix, order='F',flags=['external_loop']): 
		degen = x[0]
		AA = x[-1]

		# Undefined Ancestra Allele. Try to clean out of the loop
		if(AA == 'N' or AA == '-'):
			next
		elif('N' in x[1:-1] or '-' in x[1:-1]):
			next
		# Monomorphic sites. Try to clean out of the loop
		elif(np.unique(x[1:][np.where(x[1:]!='N')]).shape[0] == 1):
			next
		else:
			pol = x[1:-1]
			if(degen == '4'):
				functionalClass = '4fold'
			else:
				functionalClass = '0fold'

			# Check if pol != AA and monomorphic
			if((np.unique(pol).shape[0] == 1) and (np.unique(pol)[0] != AA)):
				div = 1; AF = 0
				tmp = [AF,div,functionalClass]
				output.append(tmp)
			else:
				AN = x[1:-1].shape[0]
				AC = pd.DataFrame(data=np.unique(x[1:-1], return_counts=True)[1],index=np.unique(x[1:-1], return_counts=True)[0])
				div = 0
				if(AA not in AC.index):
					next
				else:
					AC = AC[AC.index!=AA]
					if(len(AC) == 0):
						next
					else:
						AF = AC.iloc[0]/AN
						AF = AF.iloc[0]
				tmp = [AF,div,functionalClass]
				output.append(tmp)
	return(output)
def formatSfs(sequenceMatrix,rawSfsOutput,dafFile,divFile,path,append=True):

	df = pd.DataFrame(rawSfsOutput)
	df['id'] = 'uploaded'
	df.columns = ['derivedAlleleFrequency','d','functionalClass','id']

	# Extract divergence data
	div = df[['id','functionalClass','d']]
	div = div[div['d']!=0]
	div = div.groupby(['id','functionalClass'])['d'].count().reset_index()
	div = div.pivot_table(index=['id'],columns=['functionalClass'],values='d').reset_index()
	try:
		div = div[['0fold','4fold']]
	except:
		if('4fold' in div.columns):
			div = div[['4fold']]
			div['0fold'] = 0
		elif('0fold' in div.columns):
			div = div[['0fold']]
			div['4fold'] = 0

	div['mi'] =  sequenceMatrix[0][np.where(sequenceMatrix[0]=='0')].shape[0]
	div['m0'] =  sequenceMatrix[0][np.where(sequenceMatrix[0]=='4')].shape[0]
	div.columns = ['Di','D0','mi','m0']
	# div = div.pivot_table(index=['functionalClass'],columns=['functionalClass'],values='div').reset_index()

	# Create SFS pd.DataFrame by functionClass and 20 frequency bin
	daf = df[df['d']!=1][['derivedAlleleFrequency','functionalClass','id']]
	bins = np.arange(0.025,1.05,0.05)
	labels = np.arange(0.025,1.0,0.05).tolist()
	daf['categories'] = pd.cut(daf['derivedAlleleFrequency'],bins=bins,labels=labels)
	daf = daf.groupby(['functionalClass','id','categories']).count().reset_index()
	sfs = pd.DataFrame({'daf':daf['categories'].unique(),'P0':daf[daf['functionalClass']=='4fold']['derivedAlleleFrequency'].reset_index(drop=True),'Pi':daf[daf['functionalClass']=='0fold']['derivedAlleleFrequency'].reset_index(drop=True)})

	sfs = sfs[['daf','P0','Pi']]
	sfs['P0'] = sfs['P0'].fillna(0)
	sfs['Pi'] = sfs['Pi'].fillna(0)
	sfs['daf'] = sfs['daf'].apply(lambda x: round(x,3))

	if(append is True):
		sfs.to_csv(path + dafFile,sep='\t',header=True,index=False,mode='a')
		div.to_csv(path + divFile,sep='\t',header=True,index=False,mode='a')
	else:
		sfs.to_csv(path + dafFile,sep='\t',header=True,index=False)
		div.to_csv(path + divFile,sep='\t',header=True,index=False)

start = time.time()
analysis=sys.argv[1]
multiFasta = sys.argv[2]
dafFile = sys.argv[3]
divFile = sys.argv[4]
codonTable = sys.argv[5]


# Open multi-Fasta
file = px.Fasta('/var/www/html/imkt-dev/files/temporal/' + multiFasta,duplicate_action='first',sequence_always_upper=True,read_long_names=True)

# Create ndarray with sequences
multiFastaMatrix = sequencesToMatrix(file)

if(multiFastaMatrix.shape[0] < 4):
	print('numberOfLines')
	sys.exit('numberOfLines')

# Estimating SFS
rawSfs = uSfsFromFasta(multiFastaMatrix)

formatSfs(multiFastaMatrix,rawSfs,dafFile,divFile,'/var/www/html/imkt-dev/files/temporal/')

