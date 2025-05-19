#!/bin/python

'''
Usage: normalise_top_blast_hits_ssw.py [G1G2.blastp] [G1G2.blastp] [G1G2.SWscore] [G1G2.SWscore]

This script is used to isolate normalised SW scores for all top blast hits
in the species comparison. 


'''

from __future__ import division
import sys
import re , numpy
import os
from collections import Counter

def usage():
    sys.stderr.write(__doc__)
    sys.exit(1)
    
#

def usage():
    sys.stderr.write(__doc__)
    sys.exit(1)
## MATH FUNCTIONS

def calc_median(lst):
    return numpy.median(numpy.array(lst))
#
def calc_mean(lst):
    return numpy.mean(numpy.array(lst))
    
#
def calc_stdev(lst):
    return numpy.std(numpy.array(lst))
#
def calc_mad(lst):
    med = calc_median(lst)
    arr = numpy.array(lst)
    return numpy.median(numpy.absolute(arr - med))
#
blosum50 = {
'*': {'*': 1, 'A': -5, 'C': -5, 'B': -5, 'E': -5, 'D': -5, 'G': -5, 'F': -5, 'I': -5, 'H': -5, 'K': -5, 'M': -5, 'L': -5, 'N': -5, 'Q': -5, 'P': -5, 'S': -5, 'R': -5, 'T': -5, 'W': -5, 'V': -5, 'Y': -5, 'X': -5, 'Z': -5}, 
'A': {'*': -5, 'A': 5, 'C': -1, 'B': -2, 'E': -1, 'D': -2, 'G': 0, 'F': -3, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -2, 'N': -1, 'Q': -1, 'P': -1, 'S': 1, 'R': -2, 'T': 0, 'W': -3, 'V': 0, 'Y': -2, 'X': -1, 'Z': -1}, 
'C': {'*': -5, 'A': -1, 'C': 13, 'B': -3, 'E': -3, 'D': -4, 'G': -3, 'F': -2, 'I': -2, 'H': -3, 'K': -3, 'M': -2, 'L': -2, 'N': -2, 'Q': -3, 'P': -4, 'S': -1, 'R': -4, 'T': -1, 'W': -5, 'V': -1, 'Y': -3, 'X': -1, 'Z': -3}, 
'B': {'*': -5, 'A': -2, 'C': -3, 'B': 6, 'E': 1, 'D': 6, 'G': -1, 'F': -4, 'I': -4, 'H': 0, 'K': 0, 'M': -3, 'L': -4, 'N': 5, 'Q': 0, 'P': -2, 'S': 0, 'R': -1, 'T': 0, 'W': -5, 'V': -3, 'Y': -3, 'X': -1, 'Z': 1}, 
'E': {'*': -5, 'A': -1, 'C': -3, 'B': 1, 'E': 6, 'D': 2, 'G': -3, 'F': -3, 'I': -4, 'H': 0, 'K': 1, 'M': -2, 'L': -3, 'N': 0, 'Q': 2, 'P': -1, 'S': -1, 'R': 0, 'T': -1, 'W': -3, 'V': -3, 'Y': -2, 'X': -1, 'Z': 5}, 
'D': {'*': -5, 'A': -2, 'C': -4, 'B': 6, 'E': 2, 'D': 8, 'G': -1, 'F': -5, 'I': -4, 'H': -1, 'K': -1, 'M': -4, 'L': -4, 'N': 2, 'Q': 0, 'P': -1, 'S': 0, 'R': -2, 'T': -1, 'W': -5, 'V': -4, 'Y': -3, 'X': -1, 'Z': 1}, 
'G': {'*': -5, 'A': 0, 'C': -3, 'B': -1, 'E': -3, 'D': -1, 'G': 8, 'F': -4, 'I': -4, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -2, 'P': -2, 'S': 0, 'R': -3, 'T': -2, 'W': -3, 'V': -4, 'Y': -3, 'X': -1, 'Z': -2}, 
'F': {'*': -5, 'A': -3, 'C': -2, 'B': -4, 'E': -3, 'D': -5, 'G': -4, 'F': 8, 'I': 0, 'H': -1, 'K': -4, 'M': 0, 'L': 1, 'N': -4, 'Q': -4, 'P': -4, 'S': -3, 'R': -3, 'T': -2, 'W': 1, 'V': -1, 'Y': 4, 'X': -1, 'Z': -4}, 
'I': {'*': -5, 'A': -1, 'C': -2, 'B': -4, 'E': -4, 'D': -4, 'G': -4, 'F': 0, 'I': 5, 'H': -4, 'K': -3, 'M': 2, 'L': 2, 'N': -3, 'Q': -3, 'P': -3, 'S': -3, 'R': -4, 'T': -1, 'W': -3, 'V': 4, 'Y': -1, 'X': -1, 'Z': -3}, 
'H': {'*': -5, 'A': -2, 'C': -3, 'B': 0, 'E': 0, 'D': -1, 'G': -2, 'F': -1, 'I': -4, 'H': 10, 'K': 0, 'M': -1, 'L': -3, 'N': 1, 'Q': 1, 'P': -2, 'S': -1, 'R': 0, 'T': -2, 'W': -3, 'V': -4, 'Y': 2, 'X': -1, 'Z': 0}, 
'K': {'*': -5, 'A': -1, 'C': -3, 'B': 0, 'E': 1, 'D': -1, 'G': -2, 'F': -4, 'I': -3, 'H': 0, 'K': 6, 'M': -2, 'L': -3, 'N': 0, 'Q': 2, 'P': -1, 'S': 0, 'R': 3, 'T': -1, 'W': -3, 'V': -3, 'Y': -2, 'X': -1, 'Z': 1}, 
'M': {'*': -5, 'A': -1, 'C': -2, 'B': -3, 'E': -2, 'D': -4, 'G': -3, 'F': 0, 'I': 2, 'H': -1, 'K': -2, 'M': 7, 'L': 3, 'N': -2, 'Q': 0, 'P': -3, 'S': -2, 'R': -2, 'T': -1, 'W': -1, 'V': 1, 'Y': 0, 'X': -1, 'Z': -1}, 
'L': {'*': -5, 'A': -2, 'C': -2, 'B': -4, 'E': -3, 'D': -4, 'G': -4, 'F': 1, 'I': 2, 'H': -3, 'K': -3, 'M': 3, 'L': 5, 'N': -4, 'Q': -2, 'P': -4, 'S': -3, 'R': -3, 'T': -1, 'W': -2, 'V': 1, 'Y': -1, 'X': -1, 'Z': -3}, 
'N': {'*': -5, 'A': -1, 'C': -2, 'B': 5, 'E': 0, 'D': 2, 'G': 0, 'F': -4, 'I': -3, 'H': 1, 'K': 0, 'M': -2, 'L': -4, 'N': 7, 'Q': 0, 'P': -2, 'S': 1, 'R': -1, 'T': 0, 'W': -4, 'V': -3, 'Y': -2, 'X': -1, 'Z': 0}, 
'Q': {'*': -5, 'A': -1, 'C': -3, 'B': 0, 'E': 2, 'D': 0, 'G': -2, 'F': -4, 'I': -3, 'H': 1, 'K': 2, 'M': 0, 'L': -2, 'N': 0, 'Q': 7, 'P': -1, 'S': 0, 'R': 1, 'T': -1, 'W': -1, 'V': -3, 'Y': -1, 'X': -1, 'Z': 4}, 
'P': {'*': -5, 'A': -1, 'C': -4, 'B': -2, 'E': -1, 'D': -1, 'G': -2, 'F': -4, 'I': -3, 'H': -2, 'K': -1, 'M': -3, 'L': -4, 'N': -2, 'Q': -1, 'P': 10, 'S': -1, 'R': -3, 'T': -1, 'W': -4, 'V': -3, 'Y': -3, 'X': -1, 'Z': -1}, 
'S': {'*': -5, 'A': 1, 'C': -1, 'B': 0, 'E': -1, 'D': 0, 'G': 0, 'F': -3, 'I': -3, 'H': -1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': 0, 'P': -1, 'S': 5, 'R': -1, 'T': 2, 'W': -4, 'V': -2, 'Y': -2, 'X': -1, 'Z': 0}, 
'R': {'*': -5, 'A': -2, 'C': -4, 'B': -1, 'E': 0, 'D': -2, 'G': -3, 'F': -3, 'I': -4, 'H': 0, 'K': 3, 'M': -2, 'L': -3, 'N': -1, 'Q': 1, 'P': -3, 'S': -1, 'R': 7, 'T': -1, 'W': -3, 'V': -3, 'Y': -1, 'X': -1, 'Z': 0}, 
'T': {'*': -5, 'A': 0, 'C': -1, 'B': 0, 'E': -1, 'D': -1, 'G': -2, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': 0, 'Q': -1, 'P': -1, 'S': 2, 'R': -1, 'T': 5, 'W': -3, 'V': 0, 'Y': -2, 'X': -1, 'Z': -1}, 
'W': {'*': -5, 'A': -3, 'C': -5, 'B': -5, 'E': -3, 'D': -5, 'G': -3, 'F': 1, 'I': -3, 'H': -3, 'K': -3, 'M': -1, 'L': -2, 'N': -4, 'Q': -1, 'P': -4, 'S': -4, 'R': -3, 'T': -3, 'W': 15, 'V': -3, 'Y': 2, 'X': -1, 'Z': -2}, 
'V': {'*': -5, 'A': 0, 'C': -1, 'B': -3, 'E': -3, 'D': -4, 'G': -4, 'F': -1, 'I': 4, 'H': -4, 'K': -3, 'M': 1, 'L': 1, 'N': -3, 'Q': -3, 'P': -3, 'S': -2, 'R': -3, 'T': 0, 'W': -3, 'V': 5, 'Y': -1, 'X': -1, 'Z': -3}, 
'Y': {'*': -5, 'A': -2, 'C': -3, 'B': -3, 'E': -2, 'D': -3, 'G': -3, 'F': 4, 'I': -1, 'H': 2, 'K': -2, 'M': 0, 'L': -1, 'N': -2, 'Q': -1, 'P': -3, 'S': -2, 'R': -1, 'T': -2, 'W': 2, 'V': -1, 'Y': 8, 'X': -1, 'Z': -2}, 
'X': {'*': -5, 'A': -1, 'C': -1, 'B': -1, 'E': -1, 'D': -1, 'G': -1, 'F': -1, 'I': -1, 'H': -1, 'K': -1, 'M': -1, 'L': -1, 'N': -1, 'Q': -1, 'P': -1, 'S': -1, 'R': -1, 'T': -1, 'W': -1, 'V': -1, 'Y': -1, 'X': -1, 'Z': -1}, 
'Z': {'*': -5, 'A': -1, 'C': -3, 'B': 1, 'E': 5, 'D': 1, 'G': -2, 'F': -4, 'I': -3, 'H': 0, 'K': 1, 'M': -1, 'L': -3, 'N': 0, 'Q': 4, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -2, 'V': -3, 'Y': -2, 'X': -1, 'Z': 5}}

#
blosum62 = {
'*': {'*': 1, 'A': -4, 'C': -4, 'B': -4, 'E': -4, 'D': -4, 'G': -4, 'F': -4, 'I': -4, 'H': -4, 'K': -4, 'M': -4, 'L': -4, 'N': -4, 'Q': -4, 'P': -4, 'S': -4, 'R': -4, 'T': -4, 'W': -4, 'V': -4, 'Y': -4, 'X': -4, 'Z': -4}, 
'A': {'*': -4, 'A': 4, 'C': 0, 'B': -2, 'E': -1, 'D': -2, 'G': 0, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 0, 'W': -3, 'V': 0, 'Y': -2, 'X': 0, 'Z': -1}, 
'C': {'*': -4, 'A': 0, 'C': 9, 'B': -3, 'E': -4, 'D': -3, 'G': -3, 'F': -2, 'I': -1, 'H': -3, 'K': -3, 'M': -1, 'L': -1, 'N': -3, 'Q': -3, 'P': -3, 'S': -1, 'R': -3, 'T': -1, 'W': -2, 'V': -1, 'Y': -2, 'X': -2, 'Z': -3}, 
'B': {'*': -4, 'A': -2, 'C': -3, 'B': 4, 'E': 1, 'D': 4, 'G': -1, 'F': -3, 'I': -3, 'H': 0, 'K': 0, 'M': -3, 'L': -4, 'N': 3, 'Q': 0, 'P': -2, 'S': 0, 'R': -1, 'T': -1, 'W': -4, 'V': -3, 'Y': -3, 'X': -1, 'Z': 1}, 
'E': {'*': -4, 'A': -1, 'C': -4, 'B': 1, 'E': 5, 'D': 2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': -2, 'L': -3, 'N': 0, 'Q': 2, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -3, 'V': -2, 'Y': -2, 'X': -1, 'Z': 4}, 
'D': {'*': -4, 'A': -2, 'C': -3, 'B': 4, 'E': 2, 'D': 6, 'G': -1, 'F': -3, 'I': -3, 'H': -1, 'K': -1, 'M': -3, 'L': -4, 'N': 1, 'Q': 0, 'P': -1, 'S': 0, 'R': -2, 'T': -1, 'W': -4, 'V': -3, 'Y': -3, 'X': -1, 'Z': 1}, 
'G': {'*': -4, 'A': 0, 'C': -3, 'B': -1, 'E': -2, 'D': -1, 'G': 6, 'F': -3, 'I': -4, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -2, 'P': -2, 'S': 0, 'R': -2, 'T': -2, 'W': -2, 'V': -3, 'Y': -3, 'X': -1, 'Z': -2}, 
'F': {'*': -4, 'A': -2, 'C': -2, 'B': -3, 'E': -3, 'D': -3, 'G': -3, 'F': 6, 'I': 0, 'H': -1, 'K': -3, 'M': 0, 'L': 0, 'N': -3, 'Q': -3, 'P': -4, 'S': -2, 'R': -3, 'T': -2, 'W': 1, 'V': -1, 'Y': 3, 'X': -1, 'Z': -3}, 
'I': {'*': -4, 'A': -1, 'C': -1, 'B': -3, 'E': -3, 'D': -3, 'G': -4, 'F': 0, 'I': 4, 'H': -3, 'K': -3, 'M': 1, 'L': 2, 'N': -3, 'Q': -3, 'P': -3, 'S': -2, 'R': -3, 'T': -1, 'W': -3, 'V': 3, 'Y': -1, 'X': -1, 'Z': -3}, 
'H': {'*': -4, 'A': -2, 'C': -3, 'B': 0, 'E': 0, 'D': -1, 'G': -2, 'F': -1, 'I': -3, 'H': 8, 'K': -1, 'M': -2, 'L': -3, 'N': 1, 'Q': 0, 'P': -2, 'S': -1, 'R': 0, 'T': -2, 'W': -2, 'V': -3, 'Y': 2, 'X': -1, 'Z': 0}, 
'K': {'*': -4, 'A': -1, 'C': -3, 'B': 0, 'E': 1, 'D': -1, 'G': -2, 'F': -3, 'I': -3, 'H': -1, 'K': 5, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -1, 'S': 0, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': -2, 'X': -1, 'Z': 1}, 
'M': {'*': -4, 'A': -1, 'C': -1, 'B': -3, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 1, 'H': -2, 'K': -1, 'M': 5, 'L': 2, 'N': -2, 'Q': 0, 'P': -2, 'S': -1, 'R': -1, 'T': -1, 'W': -1, 'V': 1, 'Y': -1, 'X': -1, 'Z': -1}, 
'L': {'*': -4, 'A': -1, 'C': -1, 'B': -4, 'E': -3, 'D': -4, 'G': -4, 'F': 0, 'I': 2, 'H': -3, 'K': -2, 'M': 2, 'L': 4, 'N': -3, 'Q': -2, 'P': -3, 'S': -2, 'R': -2, 'T': -1, 'W': -2, 'V': 1, 'Y': -1, 'X': -1, 'Z': -3}, 
'N': {'*': -4, 'A': -2, 'C': -3, 'B': 3, 'E': 0, 'D': 1, 'G': 0, 'F': -3, 'I': -3, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 6, 'Q': 0, 'P': -2, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -3, 'Y': -2, 'X': -1, 'Z': 0}, 
'Q': {'*': -4, 'A': -1, 'C': -3, 'B': 0, 'E': 2, 'D': 0, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': 0, 'L': -2, 'N': 0, 'Q': 5, 'P': -1, 'S': 0, 'R': 1, 'T': -1, 'W': -2, 'V': -2, 'Y': -1, 'X': -1, 'Z': 3}, 
'P': {'*': -4, 'A': -1, 'C': -3, 'B': -2, 'E': -1, 'D': -1, 'G': -2, 'F': -4, 'I': -3, 'H': -2, 'K': -1, 'M': -2, 'L': -3, 'N': -2, 'Q': -1, 'P': 7, 'S': -1, 'R': -2, 'T': -1, 'W': -4, 'V': -2, 'Y': -3, 'X': -2, 'Z': -1}, 
'S': {'*': -4, 'A': 1, 'C': -1, 'B': 0, 'E': 0, 'D': 0, 'G': 0, 'F': -2, 'I': -2, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 1, 'Q': 0, 'P': -1, 'S': 4, 'R': -1, 'T': 1, 'W': -3, 'V': -2, 'Y': -2, 'X': 0, 'Z': 0}, 
'R': {'*': -4, 'A': -1, 'C': -3, 'B': -1, 'E': 0, 'D': -2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 2, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -2, 'S': -1, 'R': 5, 'T': -1, 'W': -3, 'V': -3, 'Y': -2, 'X': -1, 'Z': 0}, 
'T': {'*': -4, 'A': 0, 'C': -1, 'B': -1, 'E': -1, 'D': -1, 'G': -2, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': 0, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 5, 'W': -2, 'V': 0, 'Y': -2, 'X': 0, 'Z': -1}, 
'W': {'*': -4, 'A': -3, 'C': -2, 'B': -4, 'E': -3, 'D': -4, 'G': -2, 'F': 1, 'I': -3, 'H': -2, 'K': -3, 'M': -1, 'L': -2, 'N': -4, 'Q': -2, 'P': -4, 'S': -3, 'R': -3, 'T': -2, 'W': 11, 'V': -3, 'Y': 2, 'X': -2, 'Z': -3}, 
'V': {'*': -4, 'A': 0, 'C': -1, 'B': -3, 'E': -2, 'D': -3, 'G': -3, 'F': -1, 'I': 3, 'H': -3, 'K': -2, 'M': 1, 'L': 1, 'N': -3, 'Q': -2, 'P': -2, 'S': -2, 'R': -3, 'T': 0, 'W': -3, 'V': 4, 'Y': -1, 'X': -1, 'Z': -2}, 
'Y': {'*': -4, 'A': -2, 'C': -2, 'B': -3, 'E': -2, 'D': -3, 'G': -3, 'F': 3, 'I': -1, 'H': 2, 'K': -2, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -3, 'S': -2, 'R': -2, 'T': -2, 'W': 2, 'V': -1, 'Y': 7, 'X': -1, 'Z': -2}, 
'X': {'*': -4, 'A': 0, 'C': -2, 'B': -1, 'E': -1, 'D': -1, 'G': -1, 'F': -1, 'I': -1, 'H': -1, 'K': -1, 'M': -1, 'L': -1, 'N': -1, 'Q': -1, 'P': -2, 'S': 0, 'R': -1, 'T': 0, 'W': -2, 'V': -1, 'Y': -1, 'X': -1, 'Z': -1}, 
'Z': {'*': -4, 'A': -1, 'C': -3, 'B': 1, 'E': 4, 'D': 1, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': -1, 'L': -3, 'N': 0, 'Q': 3, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -3, 'V': -2, 'Y': -2, 'X': -1, 'Z': 4}
}
#
#
#
def gene_length_filter(gene_list, PEP_LEN_DICT):
    lengths = []
    for f in gene_list:
        lengths.append(int(PEP_LEN_DICT[f]))
    mn,md,lw,tp,madn = mean(lengths), median(lengths), min(lengths), max(lengths), mad(lengths)
    lower = md + (3*madn)
    outliers = 0
    for i in lengths:
        if i > lower:
            outliers += 1 
    
    if mn >= md*3:
        return False
    if outliers >= len(gene_list)*0.3:
        return False
    else:
        return True
#
def FASTA(filename):
    try:
        f = file(filename)
    except IOError:  
        print "The file, %s, does not exist" % filename
        return

    order = []
    sequences = {}
      
    for line in f:
        if line.startswith('>'):
            name = line[1:].rstrip('\n').split()[0]
            #name = name.replace('_', ' ')
            order.append(name)
            sequences[name] = ''
        else:
            sequences[name] += line.rstrip('\n')
              
    return sequences
    
#
def parse_cigar(cigs):
    """
    input : parse_cigar('51M12I1M2I2M4I3M1D3M3I3M2D2M3D3M3D16M')
    output: (84, 30, 0.7368421052631579) 
    """
    CIG = {}
    CIGstr = {}
    for f in "I","M","D":
        CIG[f] = 0
        CIGstr[f] = []
    for cig in re.findall('(\d+)([MID])', cigs):
        CIG[cig[1]] += int(cig[0])
        CIGstr[cig[1]].append(int(cig[0]))
    matches = CIG['M']
    mismatches = CIG['I'] + CIG['D']
    alen = matches + mismatches
    pid = float(matches / alen) * 100
    return matches,mismatches, pid,sorted(CIGstr['M'],reverse=True)[0]
#
#

def get_similarity(local_aln):
    F = float(local_aln[0].match_frequency(local_aln[1],relative=True))
    similarity = F * 100
    return similarity

def sw_score(seq1,seq2,NUCDIR):
    qr = SSW(str(NUCDIR[seq1]),protein=True,substitution_matrix=blosum50)
    al = qr(str(NUCDIR[seq2]))
    #local_al = LSSW(str(NUCDIR[seq1]),str(NUCDIR[seq2]),protein=True,substitution_matrix=blosum50)
    #sim = get_similarity(local_al)
    lQR = len(NUCDIR[seq1])
    lTR = len(NUCDIR[seq2])
    #match,mismatch,pid,maxM = parse_cigar(al.cigar())
    #sw_norm = float(((score /alen)* pid) * 100)
    match,mismatch,pid,maxM = parse_cigar(al.cigar)
    return [seq1, seq2, match, al.optimal_alignment_score, pid]
# 
ODIR = os.path.dirname(os.path.realpath(__file__))
#

def parsem9(line):
    linex = line.split("\t")
    g1 = linex[0]
    g2 = linex[1]
    alen = int(linex[3])
    qstart = int(linex[6])
    qend = int(linex[7])
    evalue = float(linex[10])
    bit = float(linex[11])
    return g1, g2, alen, qstart, qend, evalue, bit
    
def parseSWm9(linex):
    #linex = line.split("\t")
    g2 = linex[0]
    alen = int(linex[5])
    SW = int(linex[1])
    pid = float(linex[2])
    return g2, alen, SW, pid

def process_hits(FILE,SWDIR,OUTF):
    try:
        f = open(FILE).readlines()
        sys.stdout.write("opened FILE %s \n " % FILE)
    except IOError:
        print "The file, %s, does not exist" % FILE
        return
    dirMatches = {} #Directory of all targets for any given query. Used to check for duplicate blast scores
    dirResults = {} #The actual results
    queries = []
    list_results = []
    QID = ''
    K = 0
    N = 0
    SW_set = set(SWDIR.keys())
    SW_DIR = SWDIR
    Topout = open(OUTF,'w')
    for line in f:
        if line.startswith("#"):
            pass
        if len(line.split("\t")) <= 10:
            pass
        else:
            g1,g2,alen,qs,qe,evalue,bit = parsem9(line)
            sys.stdout.write("\rQuery Number %d being processed - " % (N))
            sys.stdout.flush()
            if g1 != QID:
                list_results.sort(key=lambda x: int(x[2]), reverse=True)
                dirResults[QID] = list_results
                list_results = []
                QID = g1
                dirMatches[QID] = []
                dirResults[QID] = []
                queries.append(g1)
                K = 0
                N = N + 1
                #if N % 100 == 0:
                    #sys.stdout.write("\r%d Queries processed" % N)
                    #sys.stdout.flush()
            listM = []
            listM.extend((g1,g2,str(bit),str(evalue))) 
            AL = (qe - qs) * ALCUT
            if K == 5:
                pass
            else:
                if (alen > AL and bit > BCUT and evalue < ECUT):
                    if g2 in dirMatches[QID]:
                        pass
                    else:
                        g1g2 = str(g1) + ";" + str(g2)
                        if g1g2 in SW_set: ##Add another line here to account for high confidence BLASTP scores if it is 
                            sw_mod = SW_DIR[g1g2]
                            listM[2] = str(sw_mod)
                            Topout.write('\t'.join(listM) + "\n")
                            #listM.append(str(bit))
                            #dirResults[QID].append(listM)
                            list_results.append(listM)
                            dirMatches[QID].append(g2)
                            K = K + 1
                        else: pass
    Topout.close()
    return queries, dirMatches, dirResults, list_results

def process_sw_hits(FILE,OUTF):
    try:
        f = open(FILE).readlines()
        sys.stdout.write("opened FILE %s \n " % FILE)
    except IOError:
        print "The file, %s, does not exist" % FILE
        return
    dirMatches = {}
    QID = ''
    SWout = open(OUTF,'w')
    for line in f:
        if line.startswith(">"):
            if line == ">///":
                break
            QID = line.strip(">").strip()
            #print(QID)
            pass
        if line.startswith("#"):
            pass
        elif not line.startswith(">") or line.startswith("#"):
            linex = line.strip().split("\t")
            if len(linex) == 17:
                g1 = QID
                g2, alen, sw ,pid = parseSWm9(linex)
                sw_mod = float(((sw /alen)* pid) * 100)
                g1g2 = str(g1) + ";" + str(g2)
                SWout.write(g1g2 + "\t" + str(sw_mod) + "\n")
                dirMatches[g1g2] = int(sw_mod)
    return dirMatches

if not len(sys.argv) == 6:
    usage()
    
ECUT=float(sys.argv[5])
BCUT=30
ALCUT=0.50 #AL = (qend - qstart +1) * ALCUT <- length of alignment has to be greater than this
BFIN1 = sys.argv[1]
BFIN2 = sys.argv[2]
SFIN1 = sys.argv[3]
SFIN2 = sys.argv[4]
# BLIN1 = open(BFIN1).readlines()
# BLIN2 = open(BFIN2).readlines()
# SWIN1 = open(SFIN1).readlines()
# SWIN2 = open(SFIN2).readlines()

sw1M = process_sw_hits(SFIN1,"G1ToG2.SW")
#print "Successfully parsed " + SFIN1
sw2M = process_sw_hits(SFIN2,"G2ToG1.SW")
#print "Successfully parsed " + SFIN2

bl1q,bl1M,bl1R,bl1out = process_hits(BFIN1,sw1M,"G1ToG2.top")
bl2q,bl2M,bl2R,bl2out = process_hits(BFIN2,sw2M,"G2ToG1.top")

# G1T2 = open("G1ToG2.top","w")
# G2T1 = open("G2ToG1.top","w")

# for f in bl1q:
    # for x in bl1R[f]:
        # G1T2.write('\t'.join(x) + "\n")

# for f in bl2q:
    # for x in bl2R[f]:
        # G2T1.write('\t'.join(x) + "\n")
    
# G1T2.close()
# G2T1.close()








# dirMatches = {} #Directory of all targets for any given query. Used to check for duplicate blast scores
# dirResults = {} #The actual results
# queries = []
# QID = ''
# K = 0
# for line in IN1:
    # if line.startswith("#"):
        # pass
    # else:
        # g1,g2,len,qs,qe,eval,bit = parsem9(line)
        # if g1 != QID:
            # QID = g1
            # dirMatches[QID] = []
            # dirResults[QID] = []
            # queries.append(g1)
            # K = 0
        # listM = []
        # listM.extend((g1,g2,bit,eval)) 
        # dirMatches[QID].append(listM)
        # AL = (qe - qs) * ALCUT
        # if K == 5:
            # pass
        # else:
            # if (len > AL and bit > BCUT and eval < ECUT):
                # if g2 not in dirMatches[QID]:
                    # dirResults[QID}.append(listM)
                    # K = K + 1
