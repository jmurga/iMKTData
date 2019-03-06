import pandas as pd
import allel
import argparse
import os
import sys
import re

if __name__ == "__main__":
	'''Parse arguments and show the required inputs if only name is given to command line'''
	parser = argparse.ArgumentParser(description='Extract alleles frequencies transforming VCF positions to pandas dataframe. Subset 0fold and 4fold positions through manual recoding results. Depending on options selected vcf will be transformed directly to pd.DataFrame or it will be read as vcf dictionary, both using scikit-allel package. Please check scikit-allel to understand the whole code.')
	# Required arguments
	parser.add_argument("--data", type = str, required = True, help = "Genes coordinates to extract")
	parser.add_argument("--cds", type = str, required = True, help = "File with largest CDS coordinates by genes concatenated")
	parser.add_argument("--vcf", type = str, choices=['1000GP','Alns'], required = True, help = "Select raw data to extract variants positions")
	parser.add_argument("--chromosomes", type = str, required = True, help = "Select chromosomes to extract. To check specific chromosomes please introduce comma-separated values",choices=['all',str])
	# Optional arguments
	parser.add_argument("--populations", type = str, required = False, help = "Select populations to extract",choices=['Phase1','Phase3'])
	parser.add_argument("--seed", type = int, required = False, help = "Input seed")


	# Default arguments
	parser.add_argument("--path", type = str, default = '/home/jmurga/positiveSelectionHuman/201901/rawData/humans', help = "Path to output file")
	parser.add_argument("--sampling", type = int, default = None ,help = "None default value. Please introduce a number to resampling the original vcf file. Sampling value need to be minor than total individuals at raw vcf file.")
	# Parsing common arguments
	args = parser.parse_args()


	if (args.chromosomes == 'all'):
		args.chromosomes = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']
	else:
		args.chromosomes = [args.chromosomes]

	if (args.populations == 'Phase1'):
		# populations = ['CEU','CHB','YRI']
		args.populations = ['YRI']
	else:
		populations = ['ACB','ASW','BEB','CDX','CHS','CLM','ESN','FIN','GBR','GIH','GWD','IBS','ITU','JPT','KHV','LWK','MSL','MXL','PEL','PJL','PUR','STU','TSI','CEU','CHB','YRI']
	
	dfGenes = pd.read_csv(args.path + '/annotations/'+args.data,header = 0,sep='\t',usecols=['chr','startGene','endGene'])
	coordinates = pd.read_hdf(args.path+'/annotations/'+args.cds,'positions')
	coordinates['CHROM'] = coordinates['CHROM'].apply(lambda x: re.sub('chr','',x))


	if(args.vcf == '1000GP' and args.sampling is None):
		
		if(args.populations == None):
			parser.error('Please specify a populations options')

		vcfPath = '/data/shared/1000GP/'
		print('Analyzing all pops:')
		for pop in populations:
			print(pop)

			for i in args.chromosomes:
				print(i)
				subsetCoordinates = coordinates[coordinates['CHROM']==str(i)]
				subsetCoordinates['CHROM'] = subsetCoordinates['CHROM'].astype(str)
				subsetCoordinates['POS'] = subsetCoordinates['POS'].astype(int)

				pos = dfGenes[dfGenes['chr']=='chr'+str(i)].iloc[[0,-1]].reset_index()
				pos = str(i) + ':' + str(pos['startGene'][0]) + '-' + str(pos['endGene'][1])

				#Subseting VCF from first start to last end coordinate gene				
				rawVcf = allel.vcf_to_dataframe(vcfPath + 'chr' + str(i) + str(pop) + '.vcf.gz',fields=['CHROM','POS','AC','AN','REF','ALT','AA'],alt_number=1,region=pos,tabix='tabix').reset_index(drop=True)

				# Changing types to merge correctly
				if(len(rawVcf['AA']) > 1):
					rawVcf['AA'] = rawVcf['AA'].apply(lambda x: str(x)[:-3]).str.upper()
				else:
					rawVcf['AA'] = rawVcf['AA']

				rawVcf = rawVcf[(rawVcf['ALT'].str.len()==1) & (rawVcf['REF'].str.len()==1)]
				# Include population
				rawVcf['pop'] = pop
				
				# Extracting only variant positions
				alleleFreq = pd.merge(rawVcf,subsetCoordinates,on=['CHROM','POS'],how='inner')
				# Saving in hdf5 format
				# alleleFreq.to_hdf(args.path + '/alleleFrequencies/manualRawFrequencies'+p+'.h5','af',append=True)
				if(args.populations == 'Phase3'):
					alleleFreq.to_csv(args.path + '/alleleFrequencies/manualRawFrequencies.tab',mode='a',index=False,header=False,sep='\t')
				else:
					alleleFreq.to_csv(args.path + '/alleleFrequencies/manualRawFrequenciesPhase1.tab',mode='a',index=False,header=False,sep='\t')

	elif(args.vcf == '1000GP' and args.sampling is not None ):
		print('Resampling VCF file')
		os.mkdir(args.path+'/alleleFrequencies/sampling/'+str(args.sampling))

		variants = pd.read_csv(args.path + '/alleleFrequencies/'+ args.vep,sep='\t',header=None,names=['snp','CHROM','POS','id','transcript','type'])
		variants['CHROM'] = 'chr' + variants['CHROM'].astype(str)

		if(args.populations == None):
			parser.error('Please specify a populations options')

		vcfPath = '/data/shared/1000GP/'
		args.chromosomes = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]

		for pop in args.populations:
			print(pop)
			df = pd.DataFrame()

			for i in args.chromosomes:
				print(i)
				pos = dfGenes[dfGenes['chr']=='chr'+str(i)].iloc[[0,-1]].reset_index()
				pos = str(i) + ':' + str(pos['startGene'][0]) + '-' + str(pos['endGene'][1])

				#Subseting VCF from first start to last end coordinate gene		
				subsetSamples = pd.read_csv(vcfPath + '/samples/samples' + str(pop) + '.txt',header=None,names=['ind'])
				subsetSamples = subsetSamples['ind'].sample(args.sampling,random_state=args.seed).sort_values().tolist()

				callset = allel.read_vcf(vcfPath + 'chr' + str(i) + str(pop) + '.vcf.gz',fields=['CHROM','POS','AC','AN','REF','ALT','AA','GT','samples'],alt_number=1,samples=subsetSamples,region=pos,tabix='tabix')

				chrom = callset['variants/CHROM'].tolist()
				posVariant = callset['variants/POS'].tolist()
				ref = callset['variants/REF'].tolist()
				alt = callset['variants/ALT'].tolist()
				aa = callset['variants/AA'].tolist()

				gt = allel.GenotypeArray(callset['calldata/GT'])
				ac = gt.count_alleles(max_allele=1).tolist() 
				ac = pd.DataFrame(ac)
				ac = ac[1].tolist()
				rawVcf = pd.DataFrame({'CHROM':chrom,'POS':posVariant,'REF':ref,'ALT':alt,'AA':aa,'AC':ac})
				rawVcf['AN'] = args.sampling*2

				# Changing types to merge correctly
				subsetVariants = variants[variants['CHROM']=='chr'+str(i)]
				rawVcf['CHROM'] = 'chr' + rawVcf['CHROM'].astype(str)  
				#Subseting only biallelic positions
				
				#Extract only 0 and 4 fold variant positions
				tmp = pd.merge(rawVcf,subsetVariants,on=['CHROM','POS'],how='right').dropna()
				tmp['AA'] = tmp['AA'].apply(lambda x: str(x)[:-3]).str.upper()
				print(tmp.shape,subsetVariants.shape)
				df = df.append(tmp)

			df['pop'] = pop

			df.to_csv(args.path + '/alleleFrequencies/sampling/'+str(args.sampling)+'/alleleFrequency'+str(pop)+str(args.sampling)+'.tab',sep='\t',index=False,header=True,na_rep='NA')

	else:

		vcfPath = '/data/shared/1000GP/Alns/'
		
		for i in args.chromosomes:
			subsetCoordinates = coordinates[coordinates['CHROM']==str(i)]
			subsetCoordinates['POS'] = subsetCoordinates['POS'].astype(int)

			pos = dfGenes[dfGenes['chr']=='chr'+str(i)].iloc[[0,-1]].reset_index()
			pos = str(i) + ':' + str(pos['startGene'][0]) + '-' + str(pos['endGene'][1])

			callset = allel.read_vcf('/data/shared/1000GP/Alns/' + 'chr' + str(i) + '_aln.vcf.gz',tabix='tabix',region=pos,alt_number=1)
			chrom = callset['variants/CHROM'].tolist()
			posVariant = callset['variants/POS'].tolist()
			ref = callset['variants/REF'].tolist()
			alt = callset['variants/ALT'].tolist()
			# gt = allel.GenotypeArray(callset['calldata/GT'])

			# Multidimensional array to list, in order to include GT at df
			gt = callset['calldata/GT'].reshape(len(callset['calldata/GT']),2).tolist()

			rawVcf = pd.DataFrame({'CHROM':chrom,'POS':posVariant,'REF':ref,'ALT':alt,'GT':gt})
			rawVcf = rawVcf[rawVcf['ALT']!='']
			rawVcf = rawVcf[rawVcf['ALT']!='N']

			rawVcf['GT'] = rawVcf['GT'].apply(lambda x: '/'.join(map(str,x)))

			divPositions = pd.merge(rawVcf,subsetCoordinates,on=['CHROM','POS'],how='inner')
			divPositions.to_hdf(args.path+'/alleleFrequencies/manualDivergenceVariants.h5','div',append=True)