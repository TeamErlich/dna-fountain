"""
Copyright (C) 2016 Yaniv Erlich
License: GPLv3-or-later. See COPYING file for details.
"""
from distutils.core import setup
from Cython.Build import cythonize


from fountain import DNAFountain
from utils import screen_repeat, restricted_float, prepare
import argparse
import logging
import Colorer
import os
import sys
from sets import Set
from lfsr import lfsr
import json
from tqdm import tqdm
from preprocessing import read_file, read_file_np
import timeit

logging.basicConfig(level=logging.DEBUG)


def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file_in", help="file to encode", required = True)
    parser.add_argument("-l", "--size", help="number of information bytes per message", default = 128, type = int)
    parser.add_argument("-m", "--max_homopolymer", help="the largest number of nt in a homopolymer", default = 4, type = int)
    parser.add_argument("--gc", help="the fraction of gc content above/below 0.5 (example:0.1 means 0.4-0.6)", default = 0.2, type = restricted_float)
    parser.add_argument("--rs", help="Number of bytes for rs codes", default = 0, type = int)
    parser.add_argument("--delta", help="Degree distribution tuning parameter", default = 0.05, type = float)
    parser.add_argument("--c_dist", help = "Degree distribution tuning parameter", default = 0.1, type = float)
    parser.add_argument("--out", help = "File with DNA oligos", required = True)
    parser.add_argument("--stop", help = "Maximal number of oligos", default = None, type = int)
    parser.add_argument("--alpha", help = "How many more fragments to generate on top of first k (example: 0.1 will generate 10 percent more fragments)", default = 0.07, type = float)
    parser.add_argument("--no_fasta", help = "Print oligo without a fasta header", default = False, action='store_true')
    parser.add_argument("--rand_numpy", help = "Uses numpy random generator. Faster but not compatible with older versions", default = False, action = 'store_true')


    args = parser.parse_args()
    args.orf= None

    return(args)




def main():

    args = read_args()
    logging.info("Reading the file. This may take a few mintues")
    
    f_in, file_size = read_file(args.file_in, args.size)
    
    f = DNAFountain(file_in = f_in, 
                    file_size = file_size, 
                    chunk_size = args.size , 
                    rs = args.rs, 
                    max_homopolymer = args.max_homopolymer,
                    gc = args.gc,
                    delta = args.delta, 
                    c_dist = args.c_dist,
                    np = args.rand_numpy,
                    alpha = args.alpha, 
                    stop = args.stop)

    logging.info("Upper bounds on packets for decoding is %d (x%f)  with %f probability\n", int(json.loads(f.PRNG.debug())['K_prime']), 
                                                                                           json.loads(f.PRNG.debug())['Z'],
                                                                                           json.loads(f.PRNG.debug())['delta'])
    if (args.out == '-'):
        out = sys.stdout

    else: 
        out = open (args.out, 'w')
        pbar = tqdm(total= f.final, desc = "Valid oligos")

    prepare(args.max_homopolymer)
    
    used_bc = dict()

    
    while f.good < f.final:
        d = f.droplet()


        if f.screen(d):
            if not args.no_fasta:
                out.write(">packet {}_{}\n".format(f.good, d.degree))
            out.write("{}\n".format(d.to_human_readable_DNA()))

            if d.seed in used_bc:
                logging.error("Seed %d has been seen before\nDone", d.seed)
                sys.exit(1)

            used_bc[d.seed] = 1
         
            if (args.out != '-'):
                pbar.update()

    if (args.out != '-'):
        pbar.close()
    logging.info("Finished. Generated %d packets out of %d tries (%.3f)", f.good, f.tries, (f.good+0.0)/f.tries)

    out.close()
main()
