"""
Copyright (C) 2016 Yaniv Erlich
License: GPLv3-or-later. See COPYING file for details.
"""


def lfsr(state, mask):
    #Galois lfsr:
    result = state
    nbits = mask.bit_length()-1
    while True:
        result = (result << 1)
        xor = result >> nbits
        if xor != 0:
            result ^= mask

        yield result

def lfsr32p():
    #this function returns a hard coded polynomial (0b100000000000000000000000011000101).
    #The polynomial corresponds to 1 + x^25 + x^26 + x^30 + x^32, which is known 
    #to repeat only after 32^2-1 tries. Don't change unless you know what you are doing.
    return 0b100000000000000000000000011000101

def lfsr32s():
    #this function returns a hard coded state for the lfsr (0b001010101)
    #this state is the inital position in the register. You can change it without a major implication.
    return 0b001010101

def test():
    #run the test to see a stream of seed by the polynomial
    for pattern in lfsr(0b001, 0b100000000000000000000000011000101):
        print pattern
