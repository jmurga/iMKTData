#!/usr/bin/python
def degenerate(seq):

	codonTable4f={
		"TTC":'NNN',"TTT":'NNN',#PHE 002
	    "TTA":'NNN',"TTG":'NNN',#LEU 202
	   	"CTT":'NNT',"CTC":'NNC',"CTA":'NNA',"CTG":'NNG',#LUE 204
	   	"ATG":'NNN',#MET 000
	   	"ATT":'NNN',"ATC":'NNN',"ATA":'NNN',#ILE 003
		"GTC":'NNC',"GTT":'NNT',"GTA":'NNA',"GTG":'NNG',#VAL 004
		"TCT":'NNT',"TCA":'NNA',"TCC":'NNC',"TCG":'NNG',#SER 004
		"CCT":'NNT',"CCA":'NNA',"CCC":'NNC',"CCG":'NNG',#PRO 004
		"ACT":'NNT',"ACA":'NNA',"ACC":'NNC',"ACG":'NNG',#THR 004   
		"GCT":'NNT',"GCA":'NNA',"GCC":'NNC',"GCG":'NNG',#ALA 004
		"TAT":'NNN',"TAC":'NNN',#TYR 002
		"TAA":'NNN',"TAG":'NNN',"TGA":'NNN',#STOP 000
		"CAT":'NNN',"CAC":'NNN',#HIS 002
		"CAA":'NNN',"CAG":'NNN',#GLN 002
		"AAT":'NNN',"AAC":'NNN',#ASN 002
		"AAG":'NNN',"AAA":'NNN',#LYS 002
		"GAT":'NNN',"GAC":'NNN',#ASP 002
		"GAA":'NNN',"GAG":'NNN',#GLU 002
		"TGT":'NNN',"TGC":'NNN',#CYS 002 
		"TGG":'NNN',#TRP 000
		"CGT":'NNT',"CGG":'NNG',"CGA":'NNA',"CGC":'NNC',#ARG 204
		"AGA":'NNN',"AGG":'NNN',#ARG 002
		"AGT":'NNN',"AGC":'NNN',#SER 002
		"GGT":'NNT',"GGA":'NNA',"GGC":'NNC',"GGG":'NNG',#GLY 004
		"NNN":'NNN',
	   	"---":'---',
	   }
	codonTable0f={
	  	"TTC":'TTN',"TTT":'TTN',#PHE 002
		"TTA":'NTN',"TTG":'NTN',#LEU 202
		"CTT":'NTN',"CTC":'NTN',"CTA":'NTN',"CTG":'NTN',#LUE 204
		"ATG":'ATG',#MET 000
		"ATT":'ATN',"ATC":'ATN',"ATA":'ATN',#ILE 003
		"GTC":'GTN',"GTT":'GTN',"GTA":'GTN',"GTG":'GTN',#VAL 004
		"TCT":'TCN',"TCA":'TCN',"TCC":'TCN',"TCG":'TCN',#SER 004
		"CCT":'CCN',"CCA":'CCN',"CCC":'CCN',"CCG":'CCN',#PRO 004
		"ACT":'ACN',"ACA":'ACN',"ACC":'ACN',"ACG":'ACN',#THR 004
		"GCT":'GCN',"GCA":'GCN',"GCC":'GCN',"GCG":'GCN',#ALA 004
		"TAT":'TAN',"TAC":'TAN',#TYR 002
		"TAA":'NNN',"TAG":'NNN',"TGA":'NNN',#STOP 000
		"CAT":'CAN',"CAC":'CAN',#HIS 002
		"CAA":'CAN',"CAG":'CAN',#GLN 002
		"AAT":'AAN',"AAC":'AAN',#ASN 002
		"AAG":'AAN',"AAA":'AAN',#LYS 002
		"GAT":'GAN',"GAC":'GAN',#ASP 002
		"GAA":'GAN',"GAG":'GAN',#GLU 002
		"TGT":'TGN',"TGC":'TGN',#CYS 002
		"TGG":'TGG',#TRP 000
		"CGT":'NGN',"CGG":'NGN',"CGA":'NGN',"CGC":'NGN',#ARG 204
		"AGA":'AGN',"AGG":'AGN',#ARG 002
		"AGT":'AGN',"AGC":'AGN',#SER 002
		"GGT":'GGN',"GGA":'GGN',"GGC":'GGN',"GGG":'GGN',#GLY 004
		"NNN":'NNN',
		}
	
	degenerate0fold=''
	degenerate4fold=''
	
	for i in range(0, len(seq),3):

		codon = seq[i:i+3]
		# print(codon)
		# print(codonTable0f[codon])
		degenerate4fold += codonTable4f[codon]
		degenerate0fold += codonTable0f[codon]
		
	return(degenerate0fold,degenerate4fold)