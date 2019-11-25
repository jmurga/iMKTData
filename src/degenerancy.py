def degenerate(data):
	
	#DEGENERANCY DICTIONARIES
	degenerateCodonTable = {
		"TTC":'002',"TTT":'002',#PHE 002
		"TTA":'002',"TTG":'002',#LEU 002
		"CTT":'204',"CTC":'204',"CTA":'204',"CTG":'204',#LUE 204
		"ATG":'000',#MET 000
		"ATT":'003',"ATC":'003',"ATA":'003',#ILE 003
		"GTC":'004',"GTT":'004',"GTA":'004',"GTG":'004',#VAL 004
		"TCT":'004',"TCA":'004',"TCC":'004',"TCG":'004',#SER 004
		"CCT":'004',"CCA":'004',"CCC":'004',"CCG":'004',#PRO 004
		"ACT":'004',"ACA":'004',"ACC":'004',"ACG":'004',#THR 004
		"GCT":'004',"GCA":'004',"GCC":'004',"GCG":'004',#ALA 004
		"TAT":'002',"TAC":'002',#TYR 002
		"TAA":'STP',"TAG":'STP',"TGA":'STP',#STOP 000
		"CAT":'002',"CAC":'002',#HIS 002
		"CAA":'002',"CAG":'002',#GLN 002
		"AAT":'002',"AAC":'002',#ASN 002
		"AAG":'002',"AAA":'002',#LYS 002
		"GAT":'002',"GAC":'002',#ASP 002
		"GAA":'002',"GAG":'002',#GLU 002
		"TGT":'002',"TGC":'002',#CYS 002
		"TGG":'000',#TRP 000
		"CGT":'204',"CGG":'204',"CGA":'204',"CGC":'204',#ARG 204
		"AGA":'002',"AGG":'002',#ARG 002
		"AGT":'002',"AGC":'002',#SER 002
		"GGT":'004',"GGA":'004',"GGC":'004',"GGG":'004',#GLY 004					
	}		
	
	degenerancy = ''
	for i in range(0, len(data),3):
		codon = data[i:i+3]
		if('N' in codon or '-' in codon or 'M' in codon):
			degenerancy += codon
		else:
			degenerancy += degenerateCodonTable[codon]

	return(degenerancy)

def degenerateFullPositions(data):
	
	#DEGENERANCY DICTIONARIES
	degenerateCodonTable = {
		"TTC":'002',"TTT":'002',#PHE 002
		"TTA":'002',"TTG":'002',#LEU 002
		"CTT":'204',"CTC":'204',"CTA":'204',"CTG":'204',#LUE 204
		"ATG":'000',#MET 000
		"ATT":'003',"ATC":'003',"ATA":'003',#ILE 003
		"GTC":'004',"GTT":'004',"GTA":'004',"GTG":'004',#VAL 004
		"TCT":'004',"TCA":'004',"TCC":'004',"TCG":'004',#SER 004
		"CCT":'004',"CCA":'004',"CCC":'004',"CCG":'004',#PRO 004
		"ACT":'004',"ACA":'004',"ACC":'004',"ACG":'004',#THR 004
		"GCT":'004',"GCA":'004',"GCC":'004',"GCG":'004',#ALA 004
		"TAT":'002',"TAC":'002',#TYR 002
		"TAA":'002',"TAG":'002',"TGA":'000',#STOP 000
		"CAT":'002',"CAC":'002',#HIS 002
		"CAA":'002',"CAG":'002',#GLN 002
		"AAT":'002',"AAC":'002',#ASN 002
		"AAG":'002',"AAA":'002',#LYS 002
		"GAT":'002',"GAC":'002',#ASP 002
		"GAA":'002',"GAG":'002',#GLU 002
		"TGT":'002',"TGC":'002',#CYS 002
		"TGG":'000',#TRP 000
		"CGT":'204',"CGG":'204',"CGA":'204',"CGC":'204',#ARG 204
		"AGA":'002',"AGG":'002',#ARG 002
		"AGT":'002',"AGC":'002',#SER 002
		"GGT":'004',"GGA":'004',"GGC":'004',"GGG":'004',#GLY 004					
	}		
	
	degenerancy = ''
	for i in range(0, len(data),3):
		codon = data[i:i+3]
		if('N' in codon or '-' in codon or 'M' in codon):
			degenerancy += codon
		else:
			degenerancy += degenerateCodonTable[codon]

	return(degenerancy)
	# zeroFold = degenerancy.count('0')
	# fourFold = degenerancy.count('4')

	# return(zeroFold,fourFold)
				
