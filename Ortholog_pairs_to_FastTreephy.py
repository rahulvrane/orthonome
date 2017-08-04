#!/usr/bin/python

#  Copyright 2014 - 2017 Rahul Vivek Rane <rahulvrane@gmail.com>
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

"""
This program will take all pairwise analysis outputs and find most confident 1-1 orthlog groups
It will then test each group for sanity and output files for a systematic phylogeny reconstruction.
INPUT: 1. Gene list
2. path to pep's
3. path to nucs
4. path to S*_S*

"""
from __future__ import division
import sys
import numpy
import subprocess #used to call unix shell command (subprocess.call())
import os #OS routines
import sys 
import shutil # module for copying and archiving (gz etc)
import itertools 
import string
import argparse
import glob
import datetime
import time
import multiprocessing as mp
from collections import Counter
import ete2
ODIR = os.path.dirname(os.path.realpath(__file__))


## Housekeepers
def TIMESTAMP(ts):
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    return st
#
def runCMD(cmd):
    try:
        retcode = subprocess.call(cmd, shell=True)
        if retcode < 0:
            sys.stderr.write("[ERROR]: Child %s terminated by signal" % cmd, -retcode)
        # else:
            # sys.stderr.write("Child %s returned %i" % (cmd, retcode))
    except OSError as e:
        sys.stderr.write("[ERROR]: Execution of %s failed: %s" % (cmd, e))
        
    return retcode
#
def remove_files(list_to_remove):
    for f in list_to_remove:
        os.remove(f)
    
## FILE HANDLING

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
def concat_files(FILELIST,OUTFILE,type = "fasta"):
    if type == "fasta":
        with open(OUTFILE,'a') as outfile:
            for fname in FILELIST:
                with open(fname) as infile:
                    outfile.write(infile.read())
        DICTFA = FASTA(OUTFILE)
        return DICTFA
    elif type == "orthologs":
        OGF = 0
        with open(OUTFILE,'a') as outfile:
            for fname in FILELIST:
                OGF = OGF + 1
                with open(fname) as infile:
                    outfile.write(infile.read())
        return OGF
#
## MATH FUNCTIONS

def median(lst):
    return numpy.median(numpy.array(lst))
#
def mean(lst):
    return numpy.mean(numpy.array(lst))
    
#
def stdev(lst):
    return numpy.std(numpy.array(lst))
#
def mad(lst):
    med = median(lst)
    arr = numpy.array(lst)
    return numpy.median(numpy.absolute(arr - med))
#
#
#
#
def gene_length_filter(gene_list, PEP_LEN_DICT):
    lengths = []
    for f in gene_list:
        lengths.append(int(PEP_LEN_DICT[f]))
    mn,md,lw,tp,madn = mean(lengths), median(lengths), min(lengths), max(lengths), mad(lengths)
    lower = md + (2*madn)
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

#def gene_length_filter(gene_list, PEP_LEN_DICT):
#    lengths = []
#    for f in gene_list:
#        lengths.append(int(PEP_LEN_DICT[f]))
#    mn,md,lw,tp,madn = mean(lengths), median(lengths), min(lengths), max(lengths), mad(lengths)
#    lower = md - (1.5*madn)
#    outliers = 0
#    for i in lengths:
#        if i < lower:
#            outliers += 1 
#    
#    if mn >= md*3:
#        return False
#    if outliers >= len(gene_list)*0.3:
#        return False
#    else:
#        return True

## SEQUENCE ANALYSIS FUNCTIONS


def runMAFFT(infile,outfile,threads): #contains runCMD
    mafft = ODIR + "/Programs/MAFFT/bin/mafft"
    assert os.path.isfile(mafft)
    ALN = {}
    cmd = mafft + " --quiet --maxiterate 1000 --localpair --thread %d  %s > %s" % (threads,infile,outfile)
    retcode = runCMD(cmd)
    ALN = FASTA(outfile)
    return retcode, ALN
#
def runPAL2NAL(infile,nucfile,codontable = 1): 
    pal2nal = ODIR + "/Programs/pal2nal"
    revtrans = ODIR + "/Programs/RevTrans-1.4/revtrans_jarmo.py"
    assert os.path.isfile(pal2nal)
    ALN = {}
    outfile = infile + ".codon"
    cmd = pal2nal + " %s %s -output fasta -nostderr -codontable %s 2>/dev/null > %s" % (infile,nucfile,codontable,outfile)
    retcode = runCMD(cmd)
    if os.stat(outfile).st_size == 0:
        #try:
        #    cmd = revtrans + " %s %s -readthroughstop -match name -mtx %s 2>/dev/null > %s" % (nucfile,infile,codontable,outfile)
        #    retcode = runCMD(cmd)
        #    if os.stat(outfile).st_size == 0:
        #        retcode = 2
        #except:
        #    retcode = 2
        retcode = 2
        return retcode, ALN
    if retcode == 0:
        PAL = FASTA(outfile)
        lengths = []
        seqs = []
        for key in PAL.keys():
            lengths.append(len(PAL[key]))
            seqs.append(len(PAL[key].replace('-','')))
        if any(e <= 10 for e in seqs):
            retcode = 2
            return retcode, ALN
        NORM = min(lengths)
        ALN = {}
        for key in PAL.keys():
            ALN[key] = PAL[key][:NORM]
        return retcode, ALN
    else:
        return retcode, ALN
    #try:
    #    p = subprocess.Popen(cmd,shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    #    retcode = p.wait()
    #    if retcode != 0:
    #        sys.stderr.write("Child %s terminated by signal\n" % cmd, -retcode)
    #        return retcode, ALN
    #    else:
    #        sys.stdout.write("%s returned %i \n" % (cmd, retcode))
    #        p_out, p_err = p.communicate()
    #        ALN = read_pal2nal(p_out)
    #        return retcode, ALN
    #except OSError as e:
    #    sys.stderr.write("Execution of %s failed: %s" % (cmd, e))
    #    return 2, ALN
#
def read_pal2nal(item): #sits within runPAL2NAL

    order = []
    sequences = {}
      
    for line in item.split("\n"):
        if line.startswith('>'):
            name = line[1:].rstrip('\n')
            name = name.replace('_', ' ')
            order.append(name)
            sequences[name] = ''
        else:
            sequences[name] += line.rstrip('\n')
              
    return sequences
#
def prettify_clust(listin,g_dict,PEP_LEN_DICT,filter):
    DICTIN = {}
    if filter == True:
        stat_status = gene_length_filter(listin,PEP_LEN_DICT)
    elif filter == False:
        stat_status = True
    for f in listin:
        sppf = g_dict[f]
        DICTIN[sppf] = f
            
    return DICTIN, stat_status
#
def coding_aln(pepfile,nucfile,codontable,CODON_ALN_NUC,CODON_ALN_PEP,threads): # runs runMAFFT and runPAL2NAL and outputs a dictionary of alignments
    prf = pepfile.split(".")[0]
    maf_out = prf + ".mafft"
    ret_mafft , PEP_ALN = runMAFFT(pepfile, maf_out,threads)
    if ret_mafft == 0:
        ret_p2n, ALN_DIR = runPAL2NAL(maf_out, nucfile, codontable)
        if ret_p2n != 0:
            status = "failed"
            #remove_files([pepfile,nucfile])
            return status
        else:
            CODON_ALN_NUC[prf] = ALN_DIR
            CODON_ALN_PEP[prf] = PEP_ALN
            status = "success"
            #remove_files([pepfile,nucfile])
            return status
    else:
        status = "failed"
        #remove_files([pepfile,nucfile])
        return status
    
#
def process_geneset(ID,DICT_ALIGN,CODON_ALN_PEP,threads): #Administers codon alignments for successful genesets
    pepfile = "GeneSet"+str(ID) + ".pepf"
    nucfile = "GeneSet"+str(ID) + ".nucf"
    codontable = 1
    sys.stdout.write("\r [CODON ALIGNMENT]: Processing ID: %d " % ID)
    sys.stdout.flush()
    status = coding_aln(pepfile,nucfile,codontable ,DICT_ALIGN,CODON_ALN_PEP,threads)
    return status
#
def runFastTreeMP(infile,outfile):
    FT = ODIR + "/Programs/FastTreeMP"
    assert os.path.isfile(FT)
    cmd = FT + " -nt -gtr -gamma -spr 4 -mlacc 2 -slownni -bionj %s > %s " % (infile,outfile)
    retcode = runCMD(cmd)
    return retcode 
    
#
#
#
def main():
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('[ERROR]: error: %s\n' % message)
            self.print_help()
            sys.exit(2)
            
    # parse command line
    argparser = MyParser()
    argparser.usage='------------\n%(prog)s -l [FILE] -g [PATH] -p [PATH] -n [PATH] -t [THREADS] --noFilter/--Filter'
    argparser.description='Accepts MultiMSOAR2 output and all_vs_all blast to classify genes into summarized orthogroups and supergroups and gene duplications and birth.'
    argparser.epilog='--------------'
    argparser.add_argument("-l", required=True,metavar='FILE',dest = 'genelist',
                                               help="File with list of all genes in each species - format is 1 line per species")
    argparser.add_argument("-g", required=True,metavar='PATH',dest = 'ORGPATH',
                                               help="Path to all pairwise ortholog files {WKDIR/MultiMSOAR_inputs")
    argparser.add_argument("-p", required=True,metavar='PATH',dest = 'PEPPATH',
                                               help="Path to all peptide files {WKDIR")
    argparser.add_argument("-n", required=True,metavar='PATH',dest = 'NUCPATH',
                                               help="Path to all nucleotide files {WKDIR")
    argparser.add_argument("-t", default=1,metavar='NUM',dest = 'threads',
                                               help="Number of threads to run on")
    argparser.add_argument("--noFilter", action='store_false', default=False,dest = 'filter',
                                               help="Set to switch off dynamic filtering of ortholog groups")
    argparser.add_argument("--Filter", action='store_true', default=False,dest = 'filter',
                                               help="Set to switch off dynamic filtering of ortholog groups")
    argparser.add_argument("--PF-prep", action='store_true', default=False,dest = 'PF',
                                               help="Set to switch off dynamic filtering of ortholog groups")
    argparser.add_argument("--outgroup", required=True,metavar='PATH',dest = 'outgroup',
                                               help="Comma separated list of outgroups . eg. SP1 if 1 species\nOR\nSP1,SP2 for more than 1 outgroup")
    
    if len(sys.argv)==1:
        argparser.print_help()
        sys.exit(1)
    args = argparser.parse_args()
    
    PEPFILES=glob.glob(args.PEPPATH + "/*.pep")
    NUCFILES=glob.glob(args.NUCPATH + "/*.nuc")
    ORGFILES=glob.glob(args.ORGPATH + "/S*_S?*")
    
    
    
    ########################################################
    """
    Parse the list of genes in each species
    """
    FILE1=open(args.genelist).readlines()
    
    print ODIR

    spp = []
    gene_dict = {} #DICT{GENE:SPP}

    for f in FILE1:
        line = f.split("\t")
        species = line[0]
        spp.append(species)
        for x in range(1, len(line)):
            gene_rec = line[x]
            gene_dict[gene_rec] = species
            

    sys.stdout.write("#\ntotal number of species catalogued:" + str(len(spp)))
    sys.stdout.write("#\nTotal genes in gene_dict =" + str(len(gene_dict)) + "\n")
    NUMSPP = len(spp)
    
    # Protein list
    PEPS = concat_files(PEPFILES,"all_prots.fa",type = "fasta")
    PROT_LENGTH = {}
    for f in PEPS.keys():
        L = len(PEPS[f])
        PROT_LENGTH[f] = L
        
    
    # Nucleotide list
    NUCS = concat_files(NUCFILES,"all_nucs.fa",type = "fasta")
    #remove_files(["all_nucs.fa","all_prots.fa"])
    ########################################################
    """
    Create the MCL clusters
    """
    OGF = concat_files(ORGFILES,'orthologs_cat.txt',type = "orthologs")
#    with open('orthologs_cat.txt', 'a') as outfile:
#        for fname in ORGFILES:
#            with open(fname) as infile:
#               outfile.write(infile.read()) 
    
    # MCL command:
    
    sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ": Completed concatenating %d pairwise orthologs. Proceeding to run MCL\n" % OGF)
    command = "mcl orthologs_cat.txt --abc -te %d -I 2.0 -o ortholog_markov_clusters.txt" % int(args.threads)
    if os.path.isfile('ortholog_markov_clusters.txt'):
        sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ": MCL was run earlier\n")
        
    else:
        retcode_MCL = runCMD(command)
        if retcode_MCL:
            sys.stderr.write('[ERROR]: MCL did not return 0')
            sys.exit('Something went wrong while running MCL')
        else:
            sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ": Completed MCL run - ortholog_markov_clusters.txt\n")
    
    ########################################################
    ## WORKING - 5AUG-2015
    """
    Parse the MCL clusters 
    """
    sys.stdout.write("Opening MCL file: ortholog_markov_clusters.txt\n")
    MCL_FILE=open("ortholog_markov_clusters.txt")
    MK_CLUSTERS = {}
    ORX = 0
    PMX = 0
    if args.filter == False:
        sys.stdout.write("[WARNING:] " + TIMESTAMP(time.time()) + ": dynamic ortholog length distribution filter switched OFF. Only filtering by species participation\n")
    elif args.filter == True:
        sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ": dynamic ortholog length distribution filter switched ON.\n")
    for line in MCL_FILE:
        # if PMX == 0:
            # print line
        
        listor = line.strip().split("\t")
        if len(listor) == NUMSPP:
            
            #Dlistor, stat_status = prettify_clust(listor,gene_dict,PROT_LENGTH,"True")
            Dlistor = {}
            stat_status = True
            for f in listor:
                sppf = gene_dict[f]
                Dlistor[sppf] = f
            PMX += 1
            # if stat_status == False:
                # continue
            if len(Dlistor.keys()) == NUMSPP and stat_status == True:
                ORX += 1
                sys.stdout.write("\r[LOG:]     " + TIMESTAMP(time.time()) + ":%d clusters passed length filters - and %d failed" % (ORX,PMX))
                sys.stdout.flush()
                p_out = []
                n_out = []
                
                #OUTP = open("GeneSet"+str(ORX) + ".pepf",'w')
                OUTN = open("GeneSet"+str(ORX) + ".nucf",'w')
                
                listout = []
                raise_error = 0
                for s in range(0, NUMSPP):
                    sx = spp[s]
                    gx = Dlistor[sx]
                    
                    #OUTP.write('>' + sx + " " + gx + "\n" + PEPS[gx]+ "\n")
                    OUTN.write('>' + sx + " " + gx + "\n" + NUCS[gx]+ "\n")
                #    try:
                #        FP = ">" + sx + "\n" + PEPS[gx]+ "\n"
                #        p_out.append(FP)
                #    except:
                #        raise_error = 1
                #    try:
                #        FN = '>' + sx + "\n" + NUCS[gx]+ "\n"
                #        n_out.append(FN)
                #    except:
                #        raise_error = 1
                
                #OUTP.close()
                OUTN.close()
                #if len(p_out) == NUMSPP and len(n_out) == NUMSPP:
                #    OUTP = open("GeneSet"+str(ORX) + ".pepf",'w')
                #    OUTN = open("GeneSet"+str(ORX) + ".nucf",'w')
                #    for G in p_out:
                #        OUTP.write(G)
                #    for G in n_out:
                #        OUTN.write(G)
                #    OUTP.close()
                #    OUTN.close()
                    
    sys.stdout.write("\n")
    sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ": Number of markov clusters passing species number filter - %sd | %d FAILED \n" % (ORX,PMX))
    

    
    ######################################################## 

    """
    Align every one of the [ORX] number of alignments and pick the successes
    """
    sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ":Starting processing of successful clusters\n")
    CODON_ALN_NUC = {}
    CODON_ALN_PEP = {}
    SUCCESSES = 0
    FAILS = 0
    sys.stdout.write("\n")
    #serial implementation
    for IDX in range(1,ORX+1):
        status = process_geneset(IDX,CODON_ALN_NUC,CODON_ALN_PEP,int(args.threads))
        if status == "failed":
            FAILS += 1
        elif status == "success":
            SUCCESSES += 1
        sys.stdout.write("\r[CODON ALIGNMENT]: %s clusters successful || %s failed " % (SUCCESSES,FAILS))
        sys.stdout.flush()
        
    sys.stdout.write("\n")
    sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ":Finished processing all clusters. %s clusters produced successful alignments and %s failed \n" % (SUCCESSES,FAILS))
        
    sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ":Concatenating nucleotide and protein alignments for successful clusters\n")
    
    
    ######################################################## 
    """
    Combine all successful alignments
    """
    successful_sets = set(CODON_ALN_NUC.keys())
    
    CAT_P = {}
    CAT_N = {}
    
    order_genes = []
    partitions = []
    partitionsN = []
    NOS = 0
    PARTBREAK = 1
    PARTBREAKN = 1
    for sp in spp:
        CAT_P[sp] = ''
        CAT_N[sp] = ''
    for sets in successful_sets:
        sys.stdout.write("\r[[LOG:]     " + TIMESTAMP(time.time()) + ": processing ID %s" % sets)
        sys.stdout.flush()
        DICP = CODON_ALN_PEP[sets]
        DICN = CODON_ALN_NUC[sets]
        for sp in spp:
            CAT_P[sp] += DICP[sp] #+= line.
            CAT_N[sp] += DICN[sp]
            lenp = len(CAT_P[sp])
            lenN = len(CAT_N[sp])
        NOS = NOS + 1
        partition = "Gene%d = %d-%d;" % (NOS,PARTBREAK,lenp)
        for cod in [1,2,3]:
            addn = cod - 1 
            PBN = PARTBREAKN + addn
            partitionN = "Gene%d_pos%d = %d-%d\\3;" % (NOS,cod,PBN,lenN)
            partitionsN.append(partitionN)
        PARTBREAK = lenp +1
        PARTBREAKN = lenN + 1
        partitions.append(partition)
            
    sys.stdout.write("\n")
    #Gene1_pos1 = 1-789\3;
    #Gene1_pos2 = 2-789\3;
    #Gene1_pos3 = 3-789\3;
    #Gene2_pos1 = 790-1449\3;
    #Gene2_pos2 = 791-1449\3;
    #Gene2_pos3 = 792-1449\3;
    #Gene3_pos1 = 1450-2208\3;
    #Gene3_pos2 = 1451-2208\3;
    #Gene3_pos3 = 1452-2208\3;

    OUTP = open("CONCAT_align_pep.fasta",'w')
    OUTN = open("CONCAT_align_nuc.fasta",'w')
    OUTPART = open("CONCAT_align_pep.partitions","w")
    OUTPARTN = open("CONCAT_align_nuc.partitions","w")
    
    for s in range(0, NUMSPP):
        sx = spp[s]
        OUTP.write('>' + sx + "\n" + CAT_P[sx]+ "\n")
        OUTN.write('>' + sx + "\n" + CAT_N[sx]+ "\n")
    
    sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ":Writing partitions file for PartitionFinder2 by Rob Lanfear's group\n")
    
    OUTPART.write('\n'.join(partitions))
    OUTPARTN.write('\n'.join(partitionsN))
    OUTP.close()
    OUTN.close()
    OUTPART.close()
    OUTPARTN.close()
    
    
    sys.stdout.write("\n[LOG:]     " + TIMESTAMP(time.time()) + ":Concatenation completed for all combined fasta\n")
    
    sys.stdout.write("[OUTPUTS:] " + TIMESTAMP(time.time()) + ":Concatenated peptide sequence for all successful alignments: -> CONCAT_align_pep.fasta\n")
    sys.stdout.write("[OUTPUTS:] " + TIMESTAMP(time.time()) + ":Concatenated nucleotide sequence for all successful alignments: -> CONCAT_align_nuc.fasta\n")
    sys.stdout.write("[OUTPUTS:] " + TIMESTAMP(time.time()) + ":Partitions for CONCAT_align_pep.fasta: -> CONCAT_align_pep.partitions\n")
    sys.stdout.write("[OUTPUTS:] " + TIMESTAMP(time.time()) + ":Partitions for CONCAT_align_nuc.fasta: -> CONCAT_align_nuc.partitions\n")
    ######################################################## 

    """
    run FastTreeMP max cores to create a ML phylogeny topology..
    """
    sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ":Running FastTreeMP on multiple cores using CONCAT_align_nuc.fasta\n")
    retcode = runFastTreeMP("CONCAT_align_nuc.fasta","CONCAT_align_nuc.tree")
    if retcode != 0:
        sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ":FastTreeMP Failed... try command: FastTreeMP -wag -nosupport -bionj CONCAT_align_pep.fasta > CONCAT_align_pep.tree\n")
    elif retcode == 0:
        sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ":Completed running FastTreeMP on multiple cores using CONCAT_align_nuc.fasta\n")
        
    sys.stdout.write("[OUTPUTS:] " + TIMESTAMP(time.time()) + ":Raw ML newick tree based on successful alignments in CONCAT_align_nuc.fasta : -> CONCAT_align_nuc.tree\n")
    ######################################################## 
    """
    Create the newick topology of the tree and set outgroups
    """
    sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ":Rooting input tree based on outgroups\n")
    infile = open("CONCAT_align_nuc.tree").readlines()
    t = infile[0]
    tree = ete2.Tree(t)
    
    outgroup = args.outgroup.split(',')
    if len(outgroup) == 1:
        tree.set_outgroup(outgroup[0])
    if len(outgroup) >= 2:
        anc = tree.get_common_ancestor(outgroup)
        tree.set_outgroup(anc)
    
    ts = ete2.TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_branch_support = True
    
    tree.write(format=9, outfile="CONCAT_align_nuc.nwk")
    tree.render("Phylogeny_rooted.pdf",tree_style=ts,dpi=300)
    sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ": Finished rooting phylogeny\n")
    sys.stdout.write("[OUTPUTS:] " + TIMESTAMP(time.time()) + ":Rooted ML tree for phylogeny : -> CONCAT_align_nuc.nwk\n")
    sys.stdout.write("[OUTPUTS:] " + TIMESTAMP(time.time()) + ":Rooted ML tree for phylogeny : -> Phylogeny_rooted.pdf\n")
    
    ######################################################## 
    """
    Create the phylip file and also PartitionFinder config files for protien and nucleotide
    """
    sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ":Preparing PartitionFinder input\n")
    
    
    OUTPART = open("CONCAT_align_pep.partitionfinder.cfg","w")
    OUTPARTN = open("CONCAT_align_nuc.partitionfinder.cfg","w")
    OUTPART.write("alignment = CONCAT_align_pep.phy\nuser_tree_topology = CONCAT_align_nuc.nwk;\nbranchlengths = linked;\nmodels = LG+G, LG+G+F;\nmodel_selection = AICc;\n[data_blocks]\n")
    OUTPARTN.write("alignment = CONCAT_align_nuc.phy\nuser_tree_topology = CONCAT_align_nuc.nwk;\nbranchlengths = linked;\nmodels = LG+G, LG+G+F;\nmodel_selection = AICc;\n[data_blocks]\n") #search = rcluster\n
    OUTPART.write('\n'.join(partitions))
    OUTPARTN.write('\n'.join(partitionsN))
    OUTPART.write("[schemes]\nsearch=greedy;")
    OUTPARTN.write("[schemes]\nsearch=rcluster;")
    OUTPART.close()
    OUTPARTN.close()
    
    cmd = 'bash convertFasta2Phylip.sh CONCAT_align_pep.fasta > CONCAT_align_pep.phy'
    retcode = runCMD(cmd)
    if retcode != 0:
        sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ":Conversion of FASTA to phylip failed\n")
    elif retcode == 0:
        sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ":Conversion of FASTA to phylip succeeded for PEPTIDE\n")
    
    cmd = 'bash convertFasta2Phylip.sh CONCAT_align_nuc.fasta > CONCAT_align_nuc.phy'
    retcode = runCMD(cmd)
    if retcode != 0:
        sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ":Conversion of FASTA to phylip failed\n")
    elif retcode == 0:
        sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ":Conversion of FASTA to phylip succeeded for NUCLEOTIDE\n")
    
    sys.stdout.write("[OUTPUTS:] " + TIMESTAMP(time.time()) + ":Input config file for PartitionFinder : -> CONCAT_align_pep.partitionfinder.cfg\n")
    sys.stdout.write("[OUTPUTS:] " + TIMESTAMP(time.time()) + ":Input config file for PartitionFinder : -> CONCAT_align_nuc.partitionfinder.cfg\n")
    sys.stdout.write("[OUTPUTS:] " + TIMESTAMP(time.time()) + ":Input sequence file for PartitionFinder : -> CONCAT_align_pep.phy\n")
    sys.stdout.write("[OUTPUTS:] " + TIMESTAMP(time.time()) + ":Input sequence file for PartitionFinder : -> CONCAT_align_nuc.phy\n")
    
    
    
    
    #sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ":\n")
    
    ######################################################## 
    """

    """
    
    
    #sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ":\n")
    
    ######################################################## 
    """

    """
    
    
    #sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ":\n")
    
    ######################################################## 
    """

    """
    
    
    #sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ":\n")
    #sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ":\n")
    sys.stdout.write("[LOG:]     " + TIMESTAMP(time.time()) + ": COMPLETEEEEE! \n")
    
if __name__ == "__main__":
    main()
