#!/bin/python

"""
This program will summarize fasta files for sequence length. 
Outputs total, median, mean, quartiles
"""
import sys
import numpy

def median(lst):
    return numpy.median(numpy.array(lst))
def mean(lst):
    return numpy.mean(numpy.array(lst))


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
      name = line[1:].rstrip('\n')
      name = name.replace('_', ' ')
      order.append(name)
      sequences[name] = ''
    else:
      sequences[name] += line.rstrip('\n')
            
  return sequences

lengths = {}
seqs = FASTA(sys.argv[1])
#print(len(seqs))
for f in seqs.keys():
    lengths[f] = len(seqs[f])
    
vals = lengths.values()
#print vals
med, avg, total,max,min = median(vals), mean(vals), len(vals),max(vals),min(vals)
print("median = " + str(med))
print("mean = " + str(avg))
print("sequences = " + str(total))
print("Max Length = " + str(max))
print("Min Length = " + str(min))