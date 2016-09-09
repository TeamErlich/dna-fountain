"""
Copyright (C) 2016 Yaniv Erlich
License: GPLv3-or-later. See COPYING file for details.
"""
import argparse
import logging
import Colorer
import os
import sys
import re
import json
from glass import Glass
from collections import defaultdict
from utils import dna_to_byte, split_header 
import md5
from preprocessing import read_file
from aggressive import Aggressive
from shutil import copyfile
from tqdm import tqdm


logging.basicConfig(level=logging.DEBUG)
sys.setrecursionlimit(10000000)

def read_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--file_in", help="file to decode", required = True)
	parser.add_argument("-d", "--header_size", help="number of bytes for the header", default = 4, type = int)
	parser.add_argument("-n", "--chunk_num", help="the total number of chunks in the file", required = True, type = int)
	parser.add_argument("--rs", help="number of bytes for rs codes", required = True, type = int)
	parser.add_argument("--delta", help="Degree distribution tuning parameter", default = 0.05, type = float)
	parser.add_argument("--c_dist", help = "Degree distribution tuning parameter", default = 0.1, type = float)
	parser.add_argument("--out", help = "Output file", required = True, type = str)
	parser.add_argument("--fasta", help = "Input file is FASTA", action = 'store_true')
	parser.add_argument("--no_correction", help = "Skip error correcting", action = 'store_true')
	parser.add_argument("--debug_barcodes", help = "Compare input barcodes to output", type = str)
	parser.add_argument("--gc", help="range of gc content", required = True, type = float)
	parser.add_argument("-m", "--max_homopolymer", help="the largest number of nt in a homopolymer", default = 4, type = int, required = True)
	parser.add_argument("--mock", help="Don't decode droplets. Just evaluate correctness of data", default = False, action = 'store_true')
	parser.add_argument("--max_hamming", help="How many differences between sequenced DNA and corrected DNA to tolerate", type = int, default = 100)
	parser.add_argument("--max_line", help="If defined, number of lines to read", type = int, default = None)
	parser.add_argument("--size", help = "The number of bytes of the data payload in each DNA string", default = 32, type = int)
	parser.add_argument("--rand_numpy", help = "Uses numpy random generator. Faster but not compatible with older versions", default = False, action = 'store_true')
	parser.add_argument("--truth", help = "Reading the `true` input file. Good for debuging", default = None, type = str)
	parser.add_argument("--aggressive", help = "Aggressive correction of errors using consensus in file building. Not tested", default = None, type = int)


	args = parser.parse_args()

	return(args)



def load_barcodes(args):


	valid_barcodes = dict()
	try:
		f = open(args.debug_barcodes, 'r')
	except:
		logging.error("%s file not found", args.text_file)
		sys.exit(0)

	for dna in f:
		if (re.search(r"^>", dna)):
				continue

		valid_barcodes[dna.rstrip("\n")] = 1
	return valid_barcodes




def main():

	args = read_args()

	

	if args.debug_barcodes:
		valid_barcodes = load_barcodes(args)

	truth = None
	if args.truth is not None:
		truth, file_size = read_file(args.truth, args.size)


	g = Glass(args.chunk_num, 
			  header_size = args.header_size, 
			  rs = args.rs, 
			  c_dist = args.c_dist, 
			  delta = args.delta, 
			  flag_correct = not (args.no_correction),
			  gc = args.gc,
			  max_homopolymer = args.max_homopolymer,
			  max_hamming = args.max_hamming,
			  decode = not(args.mock),
			  chunk_size = args.size,
			  np = args.rand_numpy,
			  truth = truth,
			  out = args.out
			  )

	line = 0
	errors = 0
	seen_seeds = defaultdict(int);


	#pbar = tqdm(total= args.chunk_num, desc = "Valid oligos")
	if args.file_in == '-':
		f = sys.stdin
	else:    
		try: 
			f = open(args.file_in, 'r')
		except:
			logging.error("%s file not found", args.text_file)
			sys.exit(0)



	aggressive = None
	if args.aggressive:
		
		aggressive = Aggressive(g = g, file_in = f, times = args.aggressive)


	######## Main loop
	while True:
   
		try:     
			dna = f.readline().rstrip('\n')
		except:
			logging.info("Finished reading input file!")
			break
		

		if len(dna) == 0:
			logging.info("Finished reading input file!")
			break

		if (args.fasta and re.search(r"^>", dna)):
				continue

		coverage = 0
		#when the file is in the format of coverage \t DNA
		if (len(dna.split()) == 2):
			coverage, dna = dna.split()
			####Aggresive mode
			if aggressive is not None and aggressive.turn_on(int(coverage), seen_seeds):
				best_file, value = aggressive.start()
				if best_file is not None:
					copyfile(best_file, args.out)
					logging.info("Done!")
				else:
					logging.error("Could not decode all file...")
				
				sys.exit(1)
			### End of aggressive mode

		if 'N' in dna:
			continue



		line += 1
		seed, data = g.add_dna(dna)


		if seed == -1: #reed-solomon error!
			errors += 1
		else:
			#pbar.update()
			if args.debug_barcodes:
				if not dna in valid_barcodes:
					logging.error("Seed or data %d in line %d are not valid:%s", seed, line, dna)
				else:
					seen_seeds[dna] += 1
			else:
				seen_seeds[seed] += 1       


		if line % 1000 == 0:
			logging.info("After reading %d lines, %d chunks are done. So far: %d rejections (%f) %d barcodes", line, g.chunksDone(), errors, errors/(line+0.0), g.len_seen_seed())
			pass

		if line == args.max_line:
			logging.info("Finished reading maximal number of lines")
			break
	
		if g.isDone():
			logging.info("After reading %d lines, %d chunks are done. So far: %d rejections (%f) %d barcodes", line, g.chunksDone(), errors, errors/(line+0.0), g.len_seen_seed())
			logging.info("Done!")
			break

	if not g.isDone():
		 logging.error("Could not decode all file...")
		 sys.exit(1)

	outstring = g.getString()
	f = open(args.out, 'wb')
	f.write(outstring)
	f.close()

	logging.info("MD5 is %s", md5.new(outstring).hexdigest())

	json.dump(seen_seeds, open("seen_barocdes.json",'w'), sort_keys = True, indent = 4)



main()
