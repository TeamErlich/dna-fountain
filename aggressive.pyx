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
from random import shuffle
import md5
import logging
import Colorer

class Aggressive:
    def __init__(self, g, file_in, times, min_coverage = 1):
        
        self.entries = []
        self.glass = g
        self.file_in = file_in
        self.times = times
        self.on = False
        self.min_coverage = min_coverage
        self.glass_file = None
        self.md5_dict = defaultdict(list)

    def turn_on(self, coverage, seen_seed):
        self.seen_seed = seen_seed
        #coverage is too high
        if coverage > self.min_coverage:
            return 0

        return 1

    def start(self, ):

        #turnning on!
        logging.debug("We are in aggressive mode. Hold tight")
        self.glass_file = self.glass.save()
        logging.debug("Frozen glass at %s", self.glass_file)
        self.get_all_remaining_lines()
        logging.debug("Read all remaining lines (%d)", len(self.lines))
        self.saveme('del.agg.1.tmp')
        return self.loop()


    def loop(self):
        random.seed(1)
        for i in xrange(self.times):

            logging.debug("Try %d out of %d", i+1, self.times)
            g = self.load_glass(self.glass_file)
            logging.debug("Loaded glass successfully from %s", self.glass_file)
            logging.debug("Glass status is alive? %r", g.alive())
            lines = self.shuffle_lines()
            logging.debug("Shuffled lines [we have %d lines]", len(lines))
            outstring, errors, n_line = self.reciever(lines, g)
            logging.debug("Finished. Read %d additional lines. %d were rejected", n_line, errors)

            if outstring is not None:
                self.save(i, outstring)
            else:
                logging.debug("Can't decode file in this try")

        logging.debug("Finished aggressive decoding. Let's evaluate...")
        return self.find_best()

    def find_best(self):

        best = 0
        best_file = None
        for md5_str in self.md5_dict.keys():
            value = len(self.md5_dict[md5_str])
            logging.debug("For MD5 %s, we have %d successes. For example: %s", md5_str, value, self.md5_dict[md5_str][0])
            if best < value:
                best = value
                best_file = self.md5_dict[md5_str][0]

        logging.debug("Best is %s with %d successes", best_file, best)
        return best_file, best

    def save(self, index, outstring):

        outname = ''.join([self.glass_file, str(i)])
        with open(outname, "w") as o:
            o.write(outstring)
        logging.debug("Results of decoding are at %s", outname)
        o.close()
        
        md5_str = md5.new(outstring).hexdigest()
        logging.debug("Results of decoding are at %s with MD5: %s", outname, md5_str)
        self.md5_dict[md5_str].append(outname)
        return 1

    def get_all_remaining_lines(self):
        self.lines = self.file_in.readlines()

    def shuffle_lines(self):
        lines = self.lines
        shuffle(lines)
        return lines

    def load_glass(self, name):
        with open(name, 'rb') as input:
            return pickle.load(input)
        
    def reciever(self, lines, g):

        errors = 0
        n_line = 0

        logging.debug("Starting reciever. Already %d chunks are done", g.chunksDone())
        for dna in lines:

            if 'N' in dna:
                continue
            coverage, dna = dna.rstrip('\n').split()
            seed, data = g.add_dna(dna)
            n_line += 1


            if seed == -1: #reed-solomon error!
                errors += 1

            if n_line % 100 == 0:
                logging.info("After reading %d additional lines, %d chunks are done. So far: %d rejections %d barcodes", n_line, g.chunksDone(), errors, g.len_seen_seed())
            
            if g.isDone():
                logging.debug("Done![don't get too excited...]")
                break

        if not g.isDone():
             return None, errors, n_line

        return g.getString(), errors, n_line
       
    def saveme(self, name):
        with open(name, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
        return name

