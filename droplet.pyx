"""
Copyright (C) 2016 Yaniv Erlich
License: GPLv3-or-later. See COPYING file for details.
"""
from utils import int_to_four, four_to_dna
import random
import struct
import numpy as np
from sets import Set

class Droplet:
    def __init__(self, data, seed, num_chunks = None, rs = 0, rs_obj = None, degree = None):
        #num_chunks is a list of the orignal packets numbers used to xor 
        #rs is the number of Reed Solomon symbols to add to the message

        self.data = data
        self.seed = seed
        self.num_chunks = Set(num_chunks)
        self.rs = rs
        self.rs_obj = rs_obj
        self.degree = degree

        self.DNA = None

    def chunkNums(self):
        return self.num_chunks

    def toDNA(self, flag = None):
        #this function wraps the seed, data payload, and Reed Solomon.

        if self.DNA is not None:
            return self.DNA
        
        self.DNA = int_to_four(self._package())
        return self.DNA


    def to_human_readable_DNA(self):
        #converts the DNA into a human readable [A,C,G,T
        return four_to_dna(self.toDNA())
        
    def _package(self):
        #this function converts the seed to a list of 4bytes HARD CODED!!!
        #adds the seed to the data (list of integers)
        #computes a reed solomon on the seed+data.
        #returns everything.

        seed_ord =  [ ord(c) for c in struct.pack("!I", self.seed) ]
            #converting the seed into exectly four bytes.
        message = seed_ord + self.data
        
        if self.rs > 0:
            message = self.rs_obj.encode(message) #adding RS symbols to the message

        return message
        


