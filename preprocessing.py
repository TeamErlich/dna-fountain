"""
Copyright (C) 2016 Yaniv Erlich
License: GPLv3-or-later. See COPYING file for details.
"""

import logging
import Colorer
import os
import sys
import numpy as np
import md5

logging.basicConfig(level=logging.DEBUG)


def read_file(file_in, size):
    
    chunk_size = size

    try:
        f = open(file_in, 'rb')
    except: 
        logging.error("%s file not found", file_in)
        sys.exit(0)

    data = f.read()    
    original_size = os.path.getsize(file_in)
    logging.debug("Input file has %d bytes", original_size)
    
    pad = -len(data) % chunk_size
    if pad > 0:
        logging.debug("Padded the file with %d zero to have a round number of blocks of data", pad)    
    data += "\0" * pad #zero padding.
    size = len(data)
    logging.info("File MD5 is %s", md5.new(data).hexdigest())


    data_array = [None] * (size/chunk_size) 

    logging.info("There are %d input segments", size/chunk_size)    

    for num in xrange(size/chunk_size):
        start = chunk_size * num
        end = chunk_size * (num+1)
        chunk_binary = data[start:end]

        chunk_ords = [None] * chunk_size
        for pos in xrange(chunk_size):
            chunk_ords[pos] = ord(chunk_binary[pos])

        data_array[num] = chunk_ords
    return (data_array, len(data)) 


def read_file_np(file_in, size):
    
    chunk_size = size

    try:
        f = open(file_in, 'rb')
    except: 
        logging.error("%s file not found", file_in)
        sys.exit(0)

    data = f.read()    
    original_size = os.path.getsize(file_in)
    logging.debug("Input file has %d bytes", original_size)
    
    pad = -len(data) % chunk_size
    if pad > 0:
        logging.debug("Padded the file with %d zero to have a round number of blocks of data", pad)    
    data += "\0" * pad #zero padding.
    size = len(data)
    logging.info("There are %d input segments", size/chunk_size)


    data_array = np.fromstring(data, dtype =np.uint8)
    data_array.shape = (size/chunk_size, chunk_size)
    
    return (data_array, len(data)) 
