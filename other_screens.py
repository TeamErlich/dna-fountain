from Bio.Seq import Seq
import struct
import numpy as np




class Screen:
	#this class contains method to conduct the biological screens.

	def __init__(self, 
                exDNA = False, 
                max_homopolymer = 3,
                gc = 0.05,
                orf = None):

		self.exDNA = exDNA
		self.max_homopolymer = max_homopolymer
		self.gc = gc
		self.orf = orf
		self.prepare()
		self.gc_high = 1
		self.gc_low = 0


	def prepare(self):
		#making screening faster by memoization of simple patterns and parameters
		
		if self.exDNA: #we a 6 nt alphabet
			max_alphabet = 6
		else: #we have a regular 4nt alphabet
			max_alphabet = 4

		patterns = [None] * max_alphabet
		for t in xrange(max_alphabet):
			patterns[t] = str(t) * (self.max_homopolymer + 1)
			#e.g. pattern is ["0000", "1111", "2222", "3333"]
			#for max_alphabet = 4, and max_homopolymer = 3.

		self.patterns = patterns
		self.max_alphabet = max_alphabet

		if not self.exDNA:
			#gc content if not defined for exDNA.
			self.gc_low = self.gc - 0.5
			self.gc_high = self.gc + 0.5

		return 1


	def screen(self, seq, l):
		#main method to screen a sequence

		if not self.screen_homopolymers(seq):
			return 0

		if not self.screen_gc(seq, l):
			return 0

		#if not self.screen_orf(seq):
		#	return 0

		return 1

	def screen_homopolymers(self, data):
		#detects homopoltmer patterns in data.
		#data needs to be an str. 
		#return 1 for succcess or 0 for a failed attempt.

		for pattern in self.patterns:
			if pattern in data:
				return 0
		return 1


	def screen_gc(self, data, l):
		
		gc = (data.count('1') + data.count('2') + 0.0)/ l
		if (gc < self.gc_low) or (gc > self.gc_high):
			return 0

		return 1


	def screen_orf(self, seq):
		#screen for orfs above a certain thershold
		#seq needs to be an str
		#retrun 1 for succuess or 0 foar a failed attempt

		if self.orf is None:
			return 1

		max_orf = self.orf
		seq = Seq(seq)

		for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
			for frame in range(3):
				res = (l - frame) % 3
				for pro in nuc[frame:l-res].translate(1).split("*"): 
					#translate(1) means to translate according to the regular genetic code			
					if len(pro) >= max_orf and len(pro) > 0: #avoid end of line empty string
						return 0

	
		return 1


def expandable_alphabet(array_of_bytes, l, n_symbols,  n_bytes = 8, alphabet_size = 6):
	""" takes chunks of n_bytes from array_of_bytes (of lenght l) and convert them to 
	concatanted chunks (strings) of n_symbols using alphabet_size letters. alphabet_size must be <=10"""

	res = ''
	for i in xrange(0,l, n_bytes):
		slice_array = array_of_bytes[i:i+n_bytes]
		
		number_int = long(0)
		for pos,val in enumerate(slice_array):
			number_int += long(val) * 256**(n_bytes-(pos+1)) #little endian
		
		

		res += _toDigits(n = number_int, b = alphabet_size, width = n_symbols)
		#print slice_array, number_int, toDigits(n = number_int, b = alphabet_size, width = n_symbols)
	
		#res += '\n'
	return res


def screen_repeats_x(dna, homo_length, alphabet_size):

	for t in xrange(alphabet_size):
		homo = str(t) * (homo_length + 1)

		if homo in dna:
			return 0

	return 1



def _toDigits(n, b, width):
	"""Convert a positive number n to its digit representation in base b.
	   width is the number of overall digits.
	   base MUST BE SMALLER THAN 10.
	   """

	digits = ''
	while n > 0:
		digits += str(n % b)
		n  = n // b

	digits = digits[::-1] #revsersing to little endian
	return digits.rjust(width, "0")


def _toDigits_array(n, b, width):
	"""Convert a positive number n to its digit representation in base b.
	   width is the number of overall digits.
	   returns an array of int.
	   """

	digits = list()
	while n > 0:
		digits.insert(0,int(n % b)) #little endian.
		n  = n // b

	digits = [0] * (width-len(digits)) + digits
	return digits

def dexpandable_alphabet(dna, l, n_symbols, n_bytes, alphabet_size = 6):
	
	res = list()
	for i in xrange(0,l, n_symbols):
		slice_array = dna[i:i+n_symbols]
		
		number_int = long(0)
		for pos,val in enumerate(slice_array):
			number_int += long(val) * 6**(n_symbols-(pos+1)) #little endian
		

		res += _toDigits_array(n = number_int, b = 256, width = n_bytes)
		#print slice_array, number_int, toDigits(n = number_int, b = alphabet_size, width = n_symbols)
	
		#res += '\n'
	return res


def test(size, n_symbols, n_bytes):

	
	array_of_bytes = list(np.random.randint(0, 255, size = size))
	print array_of_bytes
	ex = expandable_alphabet(array_of_bytes, size, n_symbols, n_bytes, alphabet_size = 6)
	o = dexpandable_alphabet(ex, len(ex), n_symbols, n_bytes)
	print o
	if len(set(o)-set(array_of_bytes)):
		print "Error"
	else:
		print "All good"


#test(21,65,21)