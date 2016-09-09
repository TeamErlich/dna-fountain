"""
Copyright (C) 2016 Yaniv Erlich
License: GPLv3-or-later. See COPYING file for details.
"""
from utils import *
from droplet import Droplet
from reedsolo import RSCodec
from robust_solition import PRNG
import numpy as np
import operator
import sys
from sets import Set
from collections import defaultdict
import cPickle as pickle

class Glass:
    def __init__(self, num_chunks, out, header_size = 4, 
                 rs = 0, c_dist = 0.1, delta = 0.05, 
                flag_correct = True, gc = 0.2, max_homopolymer = 4, 
                max_hamming = 100, decode = True, chunk_size = 32, np = False, truth = None):
        
        self.entries = []
        self.droplets = Set()
        self.num_chunks = num_chunks
        self.chunks = [None] * num_chunks
        self.header_size = header_size
        self.decode = decode
        self.chunk_size = chunk_size
        self.np = np
        self.chunk_to_droplets = defaultdict(set)
        self.done_segments = Set()
        self.truth = truth
        self.out = out

        self.PRNG = PRNG(K = self.num_chunks, delta = delta, c = c_dist, np = np)

        self.max_homopolymer = max_homopolymer
        self.gc = gc
        prepare(self.max_homopolymer)
        self.max_hamming = max_hamming

        self.rs = rs
        self.RSCodec = None
        self.correct = flag_correct
        self.seen_seeds = Set()
        


        if self.rs > 0:
            self.RSCodec = RSCodec(rs)
       

    def add_dna(self, dna_string):
        #header_size is in bytes
 
        data = dna_to_int_array(dna_string)


        #error correcting:
        if self.rs > 0:
            #there is an error correcting code
            if self.correct: #we want to evaluate the error correcting code
                try:
                    data_corrected = list(self.RSCodec.decode(data))
                    
                except:
                    return -1, None #could not correct the code

                #we will encode the data again to evaluate the correctness of the decoding
                data_again = list(self.RSCodec.encode(data_corrected)) #list is to convert byte array to int

                if np.count_nonzero(data != list(data_again)) > self.max_hamming: #measuring hamming distance between raw input and expected raw input
                    #too many errors to correct in decoding                    
                    return -1, None

            else: #we don't want to evaluate the error correcting code (e.g. speed)
                data_corrected  = data[0:len(data) - self.rs] #just parse out the error correcting part

        else:
            data_corrected = data

        #seed, data = split_header(data, self.header_size)
        seed_array = data_corrected[:self.header_size]
        seed = sum([   long(x)*256**i        for i, x in enumerate(seed_array[::-1])   ])
        payload = data_corrected[self.header_size:]

        #more error detection (filter seen seeds)
        if seed in self.seen_seeds:
            return -1, None
        self.add_seed(seed)

        if self.decode:
            
            #create droplet from DNA
            self.PRNG.set_seed(seed)
            blockseed, d, ix_samples = self.PRNG.get_src_blocks_wrap()
            d = Droplet(payload, seed, ix_samples)
            
            #more error detection (filter DNA that does not make sense)
            if not screen_repeat(d, self.max_homopolymer, self.gc):
                return -1, None



            self.addDroplet(d)
            

        return seed, data

    def addDroplet(self, droplet):


        self.droplets.add(droplet)
        for chunk_num in droplet.num_chunks:
            self.chunk_to_droplets[chunk_num].add(droplet) #we document for each chunk all connected droplets        
        
        self.updateEntry(droplet) #one round of message passing

        
    def updateEntry(self, droplet):

        #removing solved segments from droplets
        for chunk_num in (droplet.num_chunks & self.done_segments):
            #if self.chunks[chunk_num] is not None:
                #we solved already this input segment. 
                
            droplet.data = map(operator.xor, droplet.data, self.chunks[chunk_num])
            #subtract (ie. xor) the value of the solved segment from the droplet.
            droplet.num_chunks.remove(chunk_num)
            #cut the edge between droplet and input segment.
            self.chunk_to_droplets[chunk_num].discard(droplet)
            #cut the edge between the input segment to the droplet               

        #solving segments when the droplet have exactly 1 segment
        if len(droplet.num_chunks) == 1: #the droplet has only one input segment
            lone_chunk = droplet.num_chunks.pop() 

            self.chunks[lone_chunk] = droplet.data #assign the droplet value to the input segment (=entry[0][0])
            self.done_segments.add(lone_chunk) #add the lone_chunk to a data structure of done segments.
            if self.truth:
                self.check_truth(droplet, lone_chunk)
            self.droplets.discard(droplet) #cut the edge between the droplet and input segment
            self.chunk_to_droplets[lone_chunk].discard(droplet) #cut the edge between the input segment and the droplet
            
            #update other droplets
            for other_droplet in self.chunk_to_droplets[lone_chunk].copy():
                self.updateEntry(other_droplet)


    def getString(self):
        #return ''.join(x or ' _ ' for x in self.chunks)
        res = ''
        for x in self.chunks:
            res += ''.join(map(chr, x))
        return res

    def alive(self):
        return True

    def check_truth(self, droplet, chunk_num):
        try:
            truth_data = self.truth[chunk_num]
        except:
            print "Error. chunk:", chunk_num, " does not exist."
            quit(1)

        
        if not droplet.data == truth_data:
            #error
            print "Decoding error in ", chunk_num, ".\nInput is:", truth_data,"\nOutput is:", droplet.data,"\nDNA:", droplet.to_human_readable_DNA()
            quit(1)
        else:
            #print chunk_num, " is OK. ", self.chunksDone, " are done"
            return 1

    def add_seed(self, seed):
        self.seen_seeds.add(seed)

    def len_seen_seed(self):
        return len(self.seen_seeds)

    def isDone(self):
        if self.num_chunks - len(self.done_segments) > 0:
            return None 
        return True

    def chunksDone(self):
        return len(self.done_segments)

    def save(self):

        name =  self.out + '.glass.tmp'
        with open(name, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
        return name


