def degenerate(data):
	
	#DEGENERANCY DICTIONARIES
	degenerateCodonTable = {
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
				
