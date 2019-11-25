def reverseComplement(seq):
	# Sequence on dgn/sp contains M. No idea why
	seq = seq[::-1]
	letters = list(seq)
	baseComplement = {'A':'T','C':'G','G':'C','T':'A','N':'N','-':'-'}
#	letters = [baseComplement[base] for base in letters]
	letters = [baseComplement[base] if base in baseComplement.keys() else base for base in letters]
	seq = ''.join(letters)

	return(seq)
