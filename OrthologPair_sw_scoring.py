#!/usr/bin/python

"""
usage: OrthologPair_sw_scoring.py [G1G2.SW] [G2G1.SW] [MSOAR2_result] [MS_scores]
"""
#  OrthologPair2Scores <BLASTPOUTPUT- Tab delim format> <Pairs of genes> <output>
#
#  Copyright 2014 Rahul Vivek Rane <rahulvrane@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  


import subprocess #used to call unix shell command (subprocess.call())
import os #OS routines
import os.path #useful for manipulating pathname components
import sys 
import shutil # module for copying and archiving (gz etc)

def usage():
    sys.stderr.write(__doc__)
    sys.exit(1)
    
print("opening SW input file:" + str(sys.argv[1]))
FILE=open(sys.argv[1]).readlines() # BLAST Result File
FILEs=open(sys.argv[2]).readlines() # BLAST Result File
FILEX={}
for item in FILE:
    if not item.startswith("#"):
        itemx = item.split()
        FILEX[itemx[0]] = str(itemx[1])

print("opening SW input file:" + str(sys.argv[2]))
for item in FILEs:
    if not item.startswith("#"):
        itemx = item.split()
        FILEX[itemx[0]] = str(itemx[1])
        
print("Made Score index")
#make dict with all genes as seeds
seeds = set(FILEX.keys())

rec_dict = FILEX #Dict for blastp records

print("Finished construction of Blast Result index for " + str(len(rec_dict)) + " records and " + str(len(seeds)) + " pairs")    
print("opening Ortholog Pairs input file:" + str(sys.argv[3]))
FILE2=open(sys.argv[3]).readlines() #Ortholog Pairs file
MSR_PAIRS = []

for x in range(0, len(FILE2)):
    MSR_PAIRS.append(FILE2[x].split())

print(str(len(MSR_PAIRS)) + " pairs found")
out_file = open(sys.argv[4],'w')

for pair in MSR_PAIRS:
    g1=pair[0]
    g2=pair[1]
    g1g2=str(g1) + ";" + str(g2)
    g2g1=str(g2) + ";" + str(g1)
    if g1g2 in seeds:
        sw_mod = rec_dict[g1g2]
        out_file.write(g1 + "\t" + g2 + "\t" + sw_mod + "\n")
        #out_file.write(g1 + "\t" + g2 + "\t" + sw_mod)
    elif g2g1 in seeds:
        sw_mod =  rec_dict[g2g1]
        out_file.write(g1 + "\t" + g2 + "\t" + sw_mod + "\n")
    else:
        sw_mod = '100'
        out_file.write(g1 + "\t" + g2 + "\t" + sw_mod + "\n")
    

print("Finished processing all samples")
