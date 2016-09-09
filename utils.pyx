"""
Copyright (C) 2016 Yaniv Erlich
License: GPLv3-or-later. See COPYING file for details.
"""
from string import maketrans   # Required to call maketrans function.
import struct
import random
import os
import numpy as np
import argparse

intab = "0123"
outtab = "ACGT"

trantab = maketrans(intab, outtab)
revtab = maketrans(outtab, intab)


def charN(str, N):
    if N < len(str):
        return str[N]
    return 'X'
    
def xor(str1, str2, length):
    return ''.join(chr(ord(str1[i]) ^ ord(str2[i])) for i in xrange(length))


def xor_np(ord_array1, ord_array2):
    return np.bitwise_xor(ord_array1, ord_array2)

def xor_ord(ord_array1, ord_array2, length):
    #ord_array_res = [None] * length
    for i in xrange(length):
        #ord_array_res[i] = ord_array1[i] ^ ord_array2[i]
        ord_array1[i] = ord_array1[i] ^ ord_array2[i]
    return ord_array1#_res


def xor_bin(str1, str2):
    length = max(len(str1),len(str2))
    bytes = ''
    for i in xrange(length):
         byte_xor_dec = ord(charN(str1,i)) ^ ord(charN(str2,i))
         byte_xor_bin = "{0:08b}".format(byte_xor_dec)
         bytes = ''.join([bytes, byte_xor_bin])
    return bytes



def bin_to_dna(bin_str):
    s = ''.join(str(int(bin_str[t:t+2],2)) for t in xrange(0, len(bin_str),2)) #convert binary 2-tuple to 0,1,2,3
    return s.translate(trantab)
 

def byte_to_dna(s):
    #convert byte data (\x01 \x02) to DNA data: ACTC
    bin_data = ''.join('{0:08b}'.format(ord(s[t])) for t in xrange(0,len(s)))
    return bin_to_dna(bin_data)


def byte_from_bin(s):
    #convert string like 01010101 to string of bytes like \x01 \x02 \x03
    return ''.join(chr(int(s[t:t+8],2)) for t in xrange(0, len(s), 8))

def byte_to_int_array(s):
    #convert a strong like \x01\x02 to [1,2,]
    a = list()
    for t in xrange(0, len(s)):
        a.append(ord(s[t]))
    return a


def int_to_dna(a):
    #a is an array of integers between 0-255.
    #returns ACGGTC
    bin_data = ''.join('{0:08b}'.format(element) for element in a) #convert to a long sring of binary values
    s = ''.join(str(int(bin_data[t:t+2],2)) for t in xrange(0, len(bin_data),2)) #convert binary array to a string of 0,1,2,3
    return s.translate(trantab)


def int_to_four(a):
    #a is an array of integers between 0-255.
    #returns 0112322102
    bin_data = ''.join('{0:08b}'.format(element) for element in a) #convert to a long sring of binary values
    return ''.join(str(int(bin_data[t:t+2],2)) for t in xrange(0, len(bin_data),2)) #convert binary array to a string of 0,1,2,3

def four_to_dna(s):
    return s.translate(trantab)


def dna_to_byte(dna_str):
    #convert a string like ACTCA to a string of bytes like \x01 \x02
    num = dna_str.translate(revtab)
    s = ''.join('{0:02b}'.format(int(num[t])) for t in xrange(0, len(num),1))
    data = ''.join(chr(int(s[t:t+8],2)) for t in xrange(0, len(s), 8))

    return data

def dna_to_int_array(dna_str):
    #convert a string like ACTCA to an array of ints like [10, 2, 4]
    num = dna_str.translate(revtab)
    s = ''.join('{0:02b}'.format(int(num[t])) for t in xrange(0, len(num),1))
    data = [int(s[t:t+8],2) for t in xrange(0,len(s), 8)]

    return data


def split_header(data_str, header_bytes):
    data = data_str[header_bytes:]
    header_raw = data_str[:header_bytes]
    header_binary = ''.join('{0:08b}'.format(ord(header_raw[t])) for t in xrange(0,header_bytes))
    
    

    header = 0
    for t in xrange(0, len(header_binary)):

        header = header << 1
        header +=  int(header_binary[t])
    
    return (header, data)

def prepare(max_repeat):
    global As, Cs, Gs, Ts
    As = '0' * (max_repeat+1)
    Cs = '1' * (max_repeat+1)
    Gs = '2' * (max_repeat+1)
    Ts = '3' * (max_repeat+1)

def screen_repeat(drop, max_repeat, gc_dev):

    dna = drop.toDNA()
    return screen_repeat_dna(dna, max_repeat, gc_dev)
    

def screen_repeat_dna(dna, max_repeat, gc_dev):
    
    if As in dna or Cs in dna or Gs in dna or Ts in dna: 
        return 0

    gc = dna.count("1") + dna.count("2")  
    gc = gc/(len(dna)+0.0)

    if (gc < 0.5 - gc_dev) or (gc > 0.5 + gc_dev):
        return 0
    return 1


def hamming_distance(x,y):

    d = 0
    for t in xrange(0, max(len(x), len(y))):
        if(x[t] != y[t]):
            d += 1
    return d

def restricted_float(x):
    #helper function from http://stackoverflow.com/questions/12116685/how-can-i-require-my-python-scripts-argument-to-be-a-float-between-0-0-1-0-usin
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x






