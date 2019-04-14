import pandas as pd
import allel
import argparse
import os
import sys
import re
import pysam
import time 

if __name__ == "__main__":
	'''Parse arguments and show the required inputs if only name is given to command line'''
	parser = argparse.ArgumentParser(description='Extract alleles frequencies transforming VCF positions to pandas dataframe. Subset 0fold and 4fold positions through manual recoding results. Depending on options selected vcf will be transformed directly to pd.DataFrame or it will be read as vcf dictionary, both using scikit-allel package. Please check scikit-allel to understand the whole code.')
	# Required arguments
	parser.add_argument("--data", type = str, required = True, help = "Genes coordinates to extract")
	parser.add_argument("--vcf", type = str, choices=['1000GP','Alns'], required = True, help = "Select raw data to extract variants positions")
	# Optional arguments
	parser.add_argument("--populations", type = str, required = True, help = "Select populations to extract",choices=['Phase1','Phase3'])
	parser.add_argument("--gpFiles", type = bool, required = False, default=False, help = "Analyze whole 1000GP by Phase1 populations or Phase3",choices=[True,False])
	parser.add_argument("--seed", type = int, default = 1, help = "Input seed")
	# Default arguments
	parser.add_argument("--path", type = str, default = '/home/jmurga/positiveSelectionHuman/201901/rawData/humans', help = "Path to output file")
	parser.add_argument("--sampling",type=int, default = None ,help = "None default value. Please introduce a number to resampling the original vcf file. Sampling value need to be minor than total individuals at raw vcf file.")
	# Parsing common arguments
	args = parser.parse_args()

	args.seed = 1

	if (args.populations == 'Phase1'):
		args.populations = ['CEU','CHB','YRI']
		# args.populations = ['CEU']
	elif (args.populations == 'Phase3'):
		args.populations = ['ACB','ASW','BEB','CDX','CHS','CLM','ESN','FIN','GBR','GIH','GWD','IBS','ITU','JPT','KHV','LWK','MSL','MXL','PEL','PJL','PUR','STU','TSI','CEU','CHB','YRI']
	
	dfGenes = pd.read_csv(args.path + '/annotations/'+args.data,header = 0,sep='\t')
	dfGenes['start'] = dfGenes['start'] + 1

	if(args.vcf == '1000GP' and args.sampling is None and args.gpFiles is False):
		
		if(args.populations == None):
			parser.error('Please specify a populations options')

		vcfPath = '/data/shared/1000GP/'
		print('Analyzing all pops:')
		for pop in args.populations:
			print('\t' + pop)
			start_time = time.time()
			for index,row in dfGenes.iterrows():
				# print(index)
				#Subseting VCF from first start to last end coordinate gene             
				vcfFile = pysam.VariantFile(vcfPath+row['chr'] + str(pop) + '.vcf.gz').fetch(re.sub('chr','',row['chr']),row['start'],row['end'])
				# vcfFile = pysam.VariantFile(vcfPath+row['chr'] + '_gp.vcf.gz').fetch(re.sub('chr','',row['chr']),row['startGene'],row['endGene'])
				rawVcf = list()
				for variant in vcfFile:
					# print(variant)
					siteType = ''.join(variant.info['VT'])
					try:
						if(variant.chrom == 'Y'):
							AA = variant.info['AA'].upper()
						else:
							AA = variant.info['AA'][:-3].upper()

						if(siteType == 'SNP' and AA != '.'):
							if(variant.info['AC'][0] == 0):
								next
							else:
								
								REF = variant.ref
								ALT = variant.alts[0]
								AF  = variant.info['AC'][0]/variant.info['AN']

								# Set derived sequence and freq
								if(ALT == AA):
									dervideAllele = REF
									DAF = 1 - AF
								elif(REF == AA):
									dervideAllele = ALT
									DAF = AF    
								rawVcf.append([variant.chrom,variant.pos,REF,ALT,AA,variant.info['AC'][0],variant.info['AN'],DAF,pop,row['id']])
								# print([variant.chrom,variant.pos,REF,ALT,AA,variant.info['AC'][0],variant.info['AN'],DAF,pop])
						else:
							next
					except:
						next
				alleleFreq = pd.DataFrame(rawVcf)
				# alleleFreq = pd.merge(pd.DataFrame(rawVcf),subsetCoordinates,on=['CHROM','POS'],how='inner')
				if(len(args.populations) > 3):
					alleleFreq.to_csv(args.path + '/alleleFrequencies/sfsFromVcf.tab',mode='a',index=False,header=False,sep='\t')
				else:
					alleleFreq.to_csv(args.path + '/alleleFrequencies/sfsFromVcfPhase1.tab',mode='a',index=False,header=False,sep='\t')
			print("--- %s seconds ---" % (time.time() - start_time))
	elif(args.vcf == '1000GP' and args.gpFiles is True):

		print('+++++++++++++++++\nResampling VCF file\n+++++++++++++++++\n')

		vcfPath = '/data/shared/1000GP'
	
		df = pd.DataFrame()
		allSamples = pd.read_csv(vcfPath + '/samples/allPanel.txt',header=None,names=['ind','pop','meta','sex'],sep='\t')

		###SAMPLING
		if(args.sampling is None):
			os.makedirs(args.path+'/alleleFrequencies/sampling/',exist_ok=True)
			if(len(args.populations) == 3):
				subsetSamples = allSamples[(allSamples['pop']=='CEU') | (allSamples['pop']=='YRI') | (allSamples['pop']=='CHB')]['ind'].tolist()
			else:
				subsetSamples = allSamples['ind'].tolist()
		else:
			os.makedirs(args.path+'/alleleFrequencies/sampling/'+str(args.sampling),exist_ok=True)
			print('Resampling: ' + str(args.sampling) + ' individuals')
			subsetSamples = allSamples['ind'].sample(args.sampling,random_state=args.seed).sort_values().tolist()

		chromosomes = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
		start = time.time()
		for index,row in dfGenes.head(100).iterrows():
			print(index)
			chrNumber = re.sub('chr','',row['chr'])
			pos = str(chrNumber) + ':' + str(row['start']) + '-' + str(row['end'])

			callset = allel.read_vcf(vcfPath + '/chr' + str(chrNumber) + '_gp.vcf.gz',fields=['CHROM','POS','AC','AN','REF','ALT','AA','GT'],alt_number=1,samples=subsetSamples,region=pos,tabix='tabix')
			if(callset is None):
				next
			else:
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
				rawVcf['AN'] = len(subsetSamples)*2
				# Changing types to merge correctly
				rawVcf['AA'] = rawVcf['AA'].apply(lambda x: x[:-3].upper() if len(x) > 1 else x)
				
				# Changing types to merge correctly
				rawVcf['CHROM'] = rawVcf['CHROM'].astype(str)
				#Subseting only biallelic positions
				rawVcf = rawVcf[(rawVcf['REF'].apply(lambda x: len(x)<2)) & (rawVcf['ALT'].apply(lambda x: len(x)<2))]

				if(args.sampling is None):
					rawVcf.to_csv(args.path + '/alleleFrequencies/sampling/alleleFrequencies'+'.tab',sep='\t',index=False,header=False,na_rep='NA',mode='a')
				else:
					rawVcf.to_csv(args.path + '/alleleFrequencies/sampling/'+ str(args.sampling) +'/alleleFrequencies'+str(args.sampling)+'Phase3.tab',sep='\t',index=False,header=False,na_rep='NA',mode='a')
		
				del callset,rawVcf,chrom,posVariant,ref,alt,aa,gt,ac
		print(start - time.time())
		
	else:

		vcfPath = '/data/shared/1000GP/Alns/'
		chromosomes = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']

		for i in chromosomes:

			pos = dfGenes[dfGenes['chr']=='chr'+str(i)].iloc[[0,-1]].reset_index()
			pos = str(i) + ':' + str(pos['start'][0]) + '-' + str(pos['end'][1])

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

