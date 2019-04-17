import pandas as pd
import numpy as np
def dafWithResampling(id,data,resamplingValue,type):

	columns = ['id','POS','rawDerivedAllele','div','type']
	output = pd.DataFrame(columns=columns)

	if(data.shape[1] < resamplingValue):
		
		if(type == '4fold'):
			dafDiv = pd.DataFrame({'id':id,'POS':0,'rawDerivedAllele':0,'div':0,'type':type},index=[id])
			return(dafDiv)
		else:
			dafDiv = pd.DataFrame({'id':id,'POS':0,'rawDerivedAllele':0,'div':0,'type':type},index=[id])
			return(dafDiv)
	else:

		# print(len(data.columns.tolist()))
		for index,r in data.iterrows():
			# print(j)
			div = 0
			af = 0
			
			ref = r.iloc[[0]]
			out = r.iloc[[-1]]
			r = r.iloc[1:len(r)-1]
			
			if(r[r!='N'].dropna().shape[0] < resamplingValue):
				continue
			else:
				#Sampling
				tmp = r[r!='N'].dropna().sample(resamplingValue,replace=False)
				# Merging outgroup
				tmp = pd.concat([tmp,out])
				pos = pd.concat([ref,tmp]).reset_index(drop=True)
				
				if(pos.iloc[-1]=='N' or pos.iloc[-1]=='-'):
					continue
				elif((pos.iloc[-1] != pos.iloc[0]) & (len(pos.iloc[1:(resamplingValue+1)].unique())==1)): 
					div = 1
					af = 0
				else:
					AA = pos.iloc[resamplingValue+1]
					AN = resamplingValue
					AC = pos.iloc[1:(resamplingValue+1)].value_counts()

					if(AA not in AC.index):
						continue
					else:
						AC = AC[AC!=AC[AA]]
						if(len(AC) == 0):
							af=0
						else:
							af=AC[0]/AN
			tmp = pd.DataFrame({'id':id,'POS':index,'rawDerivedAllele':af,'div':div,'type':type},index=[id])
			tmp = tmp.reset_index(drop=True)
			output = pd.concat([output,tmp])

			
	return(output)

# def dafWithResampling(id,data,resamplingValue,type):

# 	columns = ['id','POS','rawDerivedAllele','div','type']
# 	output = pd.DataFrame(columns=columns)

# 	if(data.shape[0] < resamplingValue):
		
# 		if(type == '4fold'):
# 			dafDiv = pd.DataFrame({'id':id,'daf4f':'0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0','p0':0,'d0':0,'type':type},index=['0'])
# 			return(dafDiv)
# 		else:
# 			dafDiv = pd.DataFrame({'id':id,'daf0f':'0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0','pi':0,'di':0,'type':type},index=['0'])
# 			return(dafDiv)
# 	else:
# 		# Delete reference sequence
# 		ref = data.iloc[[0]]
# 		outgroup = data.iloc[[-1]]
# 		data = data.iloc[1:len(data)-1]
# 		# print(len(data.columns.tolist()))
		
# 		for j in data.columns.tolist():
# 			# print(j)
# 			div = 0
# 			af = 0
			
# 			if(data[[j]][data[[j]] != 'N'].dropna().shape[0] < 160):
# 				continue
# 			else:
# 				#Sampling
# 				tmp = data[[j]][data[[j]]!='N'].dropna().sample(160,replace=False)
# 				# Merging outgroup
# 				tmp = pd.concat([tmp,outgroup[[j]]])
# 				pos = pd.concat([ref[[j]],tmp]).reset_index(drop=True)
				
# 				if(pos.loc[len(pos)-1,j]=='N' or pos.loc[len(pos)-1,j]=='-'):
# 					continue
# 				elif((pos.loc[len(pos)-1,j] != pos.loc[0,j]) & (len(pos.loc[1:len(pos)-2,j].unique())==1)): 
# 					div = 1
# 					af = 0

# 				else:
# 					div = 0
# 					AA = pos.loc[len(pos)-1,j]
# 					AN = 160
# 					AC = pos.loc[1:len(pos)-2,j].value_counts()

# 					if(AA not in AC.index):
# 						af=0
# 					else:
# 						AC = AC[AC!=AC[AA]]
# 						if(len(AC) == 0):
# 							af=0
# 						else:
# 							af=AC[0]/AN
# 			tmp = pd.DataFrame({'id':id,'POS':j,'rawDerivedAllele':af,'div':div,'type':type},index=[id])
# 			tmp = tmp.reset_index(drop=True)
# 			output = pd.concat([output,tmp])

# 		# Formating output
# 		# if(type == '4fold'):
# 			div = output.groupby(['id','type','pop'])['div'].sum().reset_index()
# 			div = div[['id','div','type']]
# 			div.columns = ['id','d0','type']
# 			daf = output[['id','rawDerivedAllele','type','pop']][output['rawDerivedAllele']!=0]

# 			bins = np.arange(0,1.05,0.05)
# 			labels = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1]

# 			daf['categories'] = pd.cut(daf['rawDerivedAllele'],bins=bins,labels=labels)

# 			sfs = daf.groupby(['id','type','categories','pop']).count().reset_index()
# 			sfs['rawDerivedAllele'] = sfs['rawDerivedAllele'].fillna(0).astype(int)
# 			sfs = sfs.groupby(['id','type'])['rawDerivedAllele'].apply(list).reset_index()
# 			sfs['p'] = sum(sfs['rawDerivedAllele'][0])
# 			sfs['rawDerivedAllele'] = sfs['rawDerivedAllele'].apply(lambda x:';'.join(map(str,x)))
# 			sfs.columns = ['id','type','daf4f','p0']
# 			sfs = sfs[['id','daf4f','p0']]
# 			dafDiv = pd.merge(sfs,div,on='id')

# # # 		# else:
# 			div = output.groupby(['id','type'])['div'].sum().reset_index()
# 			div = div[['id','div','type']]
# 			div.columns = ['id','di','type']
# 			daf = output[['id','rawDerivedAllele','type']][output['rawDerivedAllele']!=0]

# 			bins = np.arange(0,1.05,0.05)
# 			labels = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1]

# 			daf['categories'] = pd.cut(daf['rawDerivedAllele'],bins=bins,labels=labels)

# 			sfs = daf.groupby(['id','type','categories']).count().reset_index()
# 			sfs['rawDerivedAllele'] = sfs['rawDerivedAllele'].fillna(0).astype(int)
# 			sfs = sfs.groupby(['id','type'])['rawDerivedAllele'].apply(list).reset_index()
# 			sfs['p'] = sfs['rawDerivedAllele'].apply(lambda x: sum(x))
# 			sfs['rawDerivedAllele'] = sfs['rawDerivedAllele'].apply(lambda x:';'.join(map(str,x)))
# 			sfs.columns = ['id','type','daf0f','pi']
# 			sfs = sfs[['id','daf0f','pi']]
# 			dafDiv = pd.merge(sfs,div,on='id')
			
# 	return(output)
