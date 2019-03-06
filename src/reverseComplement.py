def reverseComplement(seq):
	baseComplement = {'A':'T','C':'G','G':'C','T':'A','N':'N','-':'-'}
	seq = seq[::-1]
	letters = list(seq)
	letters = [baseComplement[base] for base in letters]
	seq = ''.join(letters)

	return(seq)