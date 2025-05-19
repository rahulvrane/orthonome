#!/usr/bin/python

""""  modify_info.py <info_file>
#  Method 1: Take info files and change the CHR field to comply with MSOAR2
#  How?: INPUT : INFO FILE - OUTPUT: Sorted + corrected INFO file

"""
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
import logging
import logging.handlers
import re
import natsort


if len(sys.argv) < 2:
    sys.exit('Usage: %s info_file > output_info_file' % sys.argv[0])

if not os.path.exists(sys.argv[1]):
    sys.exit('ERROR: File %s was not found!' % sys.argv[1])
    
FILE=open(sys.argv[1])
FILEX=FILE.readlines()


list_CHR = list(set(FILEX[x].split()[2] for x in range(0,len(FILEX))))
list_CHR_sorted = natsort.natsorted(list_CHR, reverse=False)
dict_CHR = {}
#print list_CHR_sorted #LIST IS OK

for x in range(1,len(list_CHR_sorted)+1):
    CHR_item = list_CHR_sorted[x-1]
    dict_CHR[CHR_item] = "chr" + str(x)
#print dict_CHR.items()  #DICT IS OK

if '2L' or '2R' or '3R' or '3L' in list_CHR_sorted:
    dmel_chrs = 'True'
else:
    dmel_chrs = 'False'





#Make a list of gene records in INFO file
MASTER_LIST=[]
for x in range(0,len(FILEX)):
    MASTER_LIST.append(FILEX[x].split())
INFOLIST = []
for fr in range(0, len(FILEX)):
    zx = FILEX[fr].split()
    chr_zx=zx[2]
    #gene_zx=zx[0]
    #gene_pre=re.sub("\d", "", gene_zx)
    #gene_suf=int(re.sub("\D", "", gene_zx))
    #zx[0]=gene_pre + str(gene_suf)
    #zx[1]=gene_pre + str(gene_suf)
    if dmel_chrs == 'True':
        REC = dict_CHR[chr_zx]
        zx[2] = REC
    else:
        REC=chr_zx
        REC_PREF=re.sub("\d", "", REC) #gives prefix for chr
        REC_SUB=re.sub("\D", "", REC) # gives only number
        if len(REC_SUB) == 1:
            if REC_SUB == "0":
                zx[2]="chr99999"
            else:
                zx[2]="chr"+REC_SUB
        elif len(REC_SUB) > 1:
            zx[2]="chr"+REC_SUB
    if len(zx) == 5:INFOLIST.append(zx)

    
INFO_SORTED = natsort.natsorted(INFOLIST, key=lambda k: (k[2], k[4]), reverse=False)
# for record in INFO_SORTED:
    # print('\t'.join(record))

bin_size= (len(INFO_SORTED) / 15) + 1

prev_chr = 1
prev_loc = 0
prev_true_loc = 0
prev_bin = 1
FINAL_INFO = []
N=1
for ZX in INFO_SORTED: 
    rec = ZX #Isolate record
    chr_no = int(re.sub("\D", "", rec[2])) #chr suffix
    loc = int(rec[4]) #scaffold location
    bin = int(N / bin_size) + 1  #bin classification
    
    if chr_no == prev_chr:
        new_loc = (loc - prev_true_loc) + int(prev_loc)
        rec[4] = str(new_loc)
        rec[2] = "chr" + str(prev_bin)
        prev_chr = chr_no
        prev_loc = new_loc
        prev_true_loc = loc
        
    elif chr_no != prev_chr:
        if bin == prev_bin:
            new_loc = int(prev_loc) + loc + 5000
            rec[2] = "chr" + str(prev_bin)
            rec[4] = str(new_loc)
            prev_chr = chr_no
            prev_loc = new_loc
            prev_true_loc = loc
        elif bin != prev_bin:
            prev_bin = bin
            new_loc = rec[4]
            rec[2] = "chr" + str(prev_bin)
            rec[4] = str(new_loc)
            prev_chr = chr_no
            prev_loc = new_loc
            prev_true_loc = loc
            
    N = N + 1
            
    FINAL_INFO.append(rec)
    
FINAL_INFO_SORTED = natsort.natsorted(FINAL_INFO, key=lambda k: (k[2], k[4]), reverse=False)
for record in FINAL_INFO_SORTED:
    print('\t'.join(record))