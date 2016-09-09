"""
Copyright (C) 2016 Yaniv Erlich
License: GPLv3-or-later. See COPYING file for details.
"""
from droplet import Droplet
from math import ceil
from utils import screen_repeat
from lfsr import lfsr, lfsr32p, lfsr32s
from robust_solition import PRNG
from reedsolo import RSCodec
import operator
import sys

import random

class DNAFountain:

    def __init__(self, 
                file_in, 
                file_size, 
                chunk_size,
                alpha, 
                stop = None,
                rs = 0, 
                c_dist = 0.1, 
                delta = 0.5, 
                np = False,
                max_homopolymer = 3,
                gc = 0.05
                ):

        #alpha is the redundency level
        #stop is whether we have a limit on the number of oligos
        #chunk_size and file_size are in bytes
        #rs is the number of bytes for reed-solomon error correcting code over gf(2^8).
        #c_dist is a parameter of the degree distribution
        #delta is a parameter of the degree distribution
        #np: should we use numpy random number generator? Faster, but incompatible for previous versions
        #max_homopolymer: the largest homopolymer allowed
        #gc: the allowable range of gc +- 50%

        #things realted to data:
        self.file_in = file_in
        self.chunk_size = chunk_size
        self.num_chunks = int(ceil(file_size / float(chunk_size)))
        self.file_size = file_size
        self.alpha = alpha
        self.stop = stop
        self.final = self.calc_stop()

        #things related to random mnumber generator
        self.lfsr = lfsr(lfsr32s(), lfsr32p()) #starting an lfsr with a certain state and a polynomial for 32bits.
        self.lfsr_l = len(    '{0:b}'.format( lfsr32p() )   ) - 1 #calculate the length of lsfr in bits 
        self.seed = self.lfsr.next()


        self.PRNG = PRNG(K = self.num_chunks, delta = delta, c = c_dist, np = np) #creating the solition distribution object
        self.PRNG.set_seed(self.seed)

        #things related to error correcting code:
        self.rs = rs #the number of symbols (bytes) to add
        self.rs_obj = RSCodec(self.rs)#initalizing an reed solomon object

        #things related to biological screens:
        self.gc = gc
        self.max_homopolymer = max_homopolymer
        self.tries = 0 #number of times we tried to create a droplet
        self.good = 0 #droplets that were screened successfully.
        self.oligo_l = self.calc_oligo_length()




    def calc_oligo_length(self):
        #return the number of nucleotides in an oligo:
        bits = self.chunk_size * 8 + self.lfsr_l + self.rs * 8
        return bits/4


    def calc_stop(self):

        if self.stop is not None:
            return self.stop
        
        stop = int(self.num_chunks*(1+self.alpha))+1
        return stop

    def droplet(self):
        #creating a droplet.
        data = None

        d, num_chunks = self.rand_chunk_nums() #creating a random list of segments.

        for num in num_chunks: #iterating over each segment
            if data is None: #first round. data payload is empty.
                data = self.chunk(num) #just copy the segment to the payload.
            else: #more rounds. Starting xoring the new segments with the payload.
                data = map(operator.xor, data, self.chunk(num))

        self.tries +=  1 #upadte counter.

        #we have a droplet:
        return Droplet(data = data, 
                       seed = self.seed, 
                       rs = self.rs,
                       rs_obj = self.rs_obj,
                       num_chunks = num_chunks,
                       degree = d)

    def chunk(self, num):
        #return the num-th segment from the file
        return self.file_in[num]

    def updateSeed(self):
        #This function creates a fresh seed for the droplet and primes the solition inverse cdf sampler
        self.seed = self.lfsr.next() #deploy one round of lfsr, and read the register.
        self.PRNG.set_seed(self.seed) #update the seed with the register

    def rand_chunk_nums(self):
        #This funcation returns a subset of segments based on the solition distribution.
        #It updates the lfsr to generates a new seed.

        self.updateSeed() #get a fresh seed and prime the solition inverse cdf sampler.
        blockseed, d, ix_samples = self.PRNG.get_src_blocks_wrap()
        return d, ix_samples #return a list of segments.

    def screen(self, droplet):

        if screen_repeat(droplet, self.max_homopolymer, self.gc):
        #if self.screen_obj.screen(droplet.toDNA(), self.oligo_l):
            self.good += 1
            return 1
        return 0



