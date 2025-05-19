#!/usr/bin/python
#
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
TO-DO
1. For each gene duplication - tag it with the superorthogroup number - basically a dictionary anew
2. Do the same.. this time characterize each orthgroup as a SuperGroup.. This will be printed alongside duplications. 
3. Figure out how we place the duplications with the orthogroups. 

##  SPP1    SPP2    SPP3    SPP4    SPP5
--  Dup1    Dup2    Dup3    NONE    NONE
--  NONE    Dup2    NONE    Dup4    Dup5?
4. Parse the entire all_vs_al blast file we input and catalogue all edges after screening for e-value > 1e-6. DIR[GENE1]={GENE2:score}

"""

from __future__ import division
import subprocess #used to call unix shell command (subprocess.call())
import os #OS routines
import os.path #useful for manipulating pathname components
import sys 
import shutil # module for copying and archiving (gz etc)
from itertools import groupby
import subprocess #used to call unix shell command (subprocess.call())
import itertools 
import string
import argparse

def proc_blast(all_blast,listBirths,output_prefix):
    print("#\nopening blast input file:" + str(all_blast))

    output_blastp = open(output_prefix + "Genebirths_all_proteins.blastp", 'w')
    output_blastpx = open(output_prefix + "Genebirths_strictly_births.blastp", 'w')

    FILEX=[]
    NB = 0
    N = 0
    NX = 0
    FILEXR = []
    ONLY_BIRTHS = []
    set_births = set(listBirths)
        
        
    with open(all_blast) as blast_file:
        for item in blast_file:
            if not item.startswith("#"):
                itemx = item.split()
                if itemx[0] in set_births or itemx[1] in set_births:
                    sppx1, sppx2 = gene_dict[itemx[0]], gene_dict[itemx[1]]
                    itemx.append(sppx1)
                    itemx.append(sppx2)
                    FILEXR.append(itemx)
                    if itemx[0] in set_births and itemx[1] in set_births:
                        ONLY_BIRTHS.append(itemx)
                        NX = NX + 1
                    N = N +1
            sys.stdout.write("\rRecords Complete - %.2f | Blast input parsed : %d hits for Orphan Gene Births, with %d strictly between Orphans " % (NB,N,NX))
            sys.stdout.flush()
            NB=NB+1
            
    sys.stdout.write("\r All Blast input parsed : %d hits for Orphan Gene Births, with %d strictly between Orphans " % (N,NX))
    sys.stdout.flush()

    print("#\nPrinting Genebirths + blastp scores out to file: Genebirths_all_proteins.blastp - AND - Genebirths_strictly_births.blastp")


    for line_item in FILEXR:
        output_blastp.write("\t".join(line_item) + "\n")
        
    for line_item in ONLY_BIRTHS:
        output_blastpx.write("\t".join(line_item) + "\n")

    print("#\nBlastp output created for records classified as GeneBirths for further use")

def main():
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)
            
    # parse command line
    argparser = MyParser()
    argparser.usage='------------\n%(prog)s -l [FILE] -a [FILE] -m [FILE] -i [FILE] -g [FILE] -o [STR]'
    argparser.description='Accepts MultiMSOAR2 output and all_vs_all blast to classify genes into summarized orthogroups and supergroups and gene duplications and birth.'
    argparser.epilog='--------------'
    argparser.add_argument("-l", required=True,metavar='FILE',dest = 'genelist',
                                               help="File with list of all genes in each species - format is 1 line per species")
    argparser.add_argument("-a", metavar='FILE',dest = 'all_blast',
                                               help="All vs All blastp result")
    argparser.add_argument("-m", required=True,metavar='FILE',dest = 'mcl_clusters',
                                               help="MCL clusters based on All vs All blastp result")
    argparser.add_argument("-i", required=True,metavar='FILE',dest = 'geneinfo',
                                               help="Geneinfo output from MultiMSOAR2")
    argparser.add_argument("-g", required=True,metavar='FILE',dest = 'orthogroups',
                                               help="Orthogroup output from MultiMSOAR2")
    argparser.add_argument("-s", required=True,metavar='FILE',dest = 'Sppindex',
                                               help="Index of species names. eg. \nSp1\t1 \n Sp2\t2")
    argparser.add_argument("-o", default='Orthogroup',metavar='Orthogroup',dest = 'output_prefix',
                                               help="output prefix")
    
    if len(sys.argv)==1:
        argparser.print_help()
        sys.exit(1)
    args = argparser.parse_args()
    
    ########################################################
    
    """
    Create supergroup ID's
    """
    alp_list = list(string.ascii_uppercase)
    SOGdf = {} # SOGdf{1:"ABCD"}
    SX = 0
    for combo in itertools.permutations(alp_list,4):
        SX = SX + 1
        lisx = list(combo)
        STRx = ''.join(lisx)
        SOGdf[SX]=STRx
    
    ########################################################
    
    """
    Parse the list of genes in each species
    """
    FILE1=open(args.genelist).readlines()

    spp = []
    gene_dict = {}

    for f in FILE1:
        line = f.split()
        species = line[0]
        spp.append(species)
        for x in range(1, len(line)):
            gene_rec = line[x]
            gene_dict[gene_rec] = species
            
    spidx = {}
    Fspdx = open(args.Sppindex).readlines()
    for line in Fspdx:
        if len(line.split())>1:
            spidx[line.strip("_").split()[0]] = line.strip("_").split()[1]
    print("#\ntotal number of species catalogued:" + str(len(spp)))
    print("#\nTotal genes in gene_dict =" + str(len(gene_dict)))
    ########################################################
    
    """
    Parse the MCL clusters 
    """
    print("Opening MCL file: %s" % args.mcl_clusters)
    MCL_FILE=open(args.mcl_clusters)
    MK_CLUSTERS = []
    
    for line in MCL_FILE:
        MK_CLUSTERS.append(line)
    
    SUPERGROUP_COUNT = len(MK_CLUSTERS)
    print("#\nTotal Ortholog supergroups =" + str(SUPERGROUP_COUNT))
    
    SUPER_GROUP = {} # DICT{GENE(FBgnXXYYAHH):SOG_ID(ABCD)}
    SUPERGROUP_LIST = []
    
    NX = 0

    for f in range(0, SUPERGROUP_COUNT):
        line = MK_CLUSTERS[f].strip("\n").split("\t")
        SX = f + 1
        SOG_ID = SOGdf[SX]
        SUPERGROUP_LIST.append(SOG_ID)
        for x in line:
            SUPER_GROUP[x] = SOG_ID
            NX = NX + 1
            
    print("#\nParsed MCL file with %d genes in clusters/singletons" % NX)

            
    
    ########################################################

    """
    Parse the GeneInfo file from MultiMSOAR2 for duplications and gene births
    """
    FILE2=open(args.geneinfo).readlines()
    
    # Initiate dictionaries for births and deaths
    gene_births = {}
    gene_duplications = {}
    
    # Initiate dictionaries - these will house per superfamily - the duplications. Therefore dict = {dict_spp{spp:['genes']}}
    SUP_births = {}
    SUP_dups = {}
    SUP_orth = {}
    
    # initiate a dict for each Supergroup
    for f in SUPERGROUP_LIST:
        SUP_births[f] = {}
        SUP_dups[f] = {}
        #SUP_orth[f] = [['####'],["SUPER_GROUP",str('SOG_'+f)]]
    
    # Initiate a list for every species...
    for sppx in spp:
        gene_births[sppx] = [sppx]
        gene_duplications[sppx] = [sppx]
        for f in SUPERGROUP_LIST:
            SUP_births[f][sppx] = ['#', sppx]
            SUP_dups[f][sppx] = ['#', sppx]
        
    #[gene_births[sppx] = [str(sppx)] for sppx in spp]
    list_births = FILE2[0].split()[2:]
    for x in range(0,len(list_births)):
        gene = list_births[x]
        species_rec = gene_dict[gene]
        sog_id = SUPER_GROUP[gene]
        gene_births[species_rec].append(gene)
        SUP_births[sog_id][species_rec].append(gene)
        
    #[gene_duplications[sppx] = [sppx] for sppx in spp]
    list_duplications = FILE2[1].split()

    for x in range(2,len(list_duplications)):
        gene = list_duplications[x]
        species_rec = gene_dict[gene]
        sog_id = SUPER_GROUP[gene]
        gene_duplications[species_rec].append(gene)
        SUP_dups[sog_id][species_rec].append(gene)
        
    output_Gbir = open(args.output_prefix + '_Genebirths.txt','w')
    output_Gdup = open(args.output_prefix + '_GeneDuplications.txt','w')

    output_Gbir.write("SPECIES\tCount\tGenes" + "\n")
    output_Gdup.write("SPECIES\tCount\tGenes" + "\n")

    for species in spp:
        dups = gene_duplications[species]
        Dups = str(len(dups) - 1)
        dups.insert(1,Dups)
        
        births = gene_births[species]
        Births = str(len(births) -1)
        births.insert(1,Births)
        output_Gbir.write("\t".join(births) + "\n")
        output_Gdup.write("\t".join(dups) + "\n")
        
    def add_lengths_sog_class(base_DIR,TYPE):
        DIR_OUT = {}
        for sog in SUPERGROUP_LIST:
            DIR_OUT[sog] = []
            dir_sog = base_DIR[sog]
            spps = dir_sog.keys()
            for spp in spps:
                list_sp = dir_sog[spp]
                count_g = len(list_sp) - 2
                list_sp[0] = str(count_g)
                list_sp.insert(0,TYPE)
                list_sp.insert(2,sog)
                DIR_OUT[sog].append(list_sp)
        return DIR_OUT
        
    SUP_NON_ORTH_OUTS = {}
    SUP_NON_ORTH_OUTS['dups'] = add_lengths_sog_class(SUP_dups,'DUP')
    SUP_NON_ORTH_OUTS['births'] = add_lengths_sog_class(SUP_births,'DNV')
        
    print("#\nGeneinfo parsed and output files created: Genebirths.txt & GeneDuplications.txt")
        
    ######################################################## 
    """
    Parse Orthgroups and write a hit_list
    """
    def sog_counter(DIR,sog_id):
        curr = DIR[sog_id]
        next = curr + 1
        DIR[sog_id] =  next
        return next
    
    Ortho_list_out = []
    Ortho_spp_freq = []
    listx = ['#','SOG_ID','OG_ID']
    for x in spp:
        listx.append(str(x))
    Ortho_list_out.append(listx)
    
    # Initiate the counter cycle
    SOG_COUNTER = {}
    for f in SUPERGROUP_LIST:
        SOG_COUNTER[f] = 0
        SUP_orth[f] = [['######'],["SUPER_GROUP",str('SOG_'+f)]]
        SUP_orth[f].append(listx)
    
    OrthGrpFILE = open(args.orthogroups).readlines()
    print("#\nOpening Orthogroups file: " + str(args.orthogroups))

    
    #sog_id = SUPER_GROUP[gene] to fill up DICT ->SUP_orth
    NXN = 1
    N = 0
    web_list = []
    for line in OrthGrpFILE:
        N = N + 1
        
        dict_spp = {}
        list_spp = []
        #Initiate empty string for each species in dict_spp
        for sppx in spp:
            dict_spp[sppx] = ''
        list_line = line.split()
        # Find out what the SOG_ID of the orthogroup is! <---!! 
        sog_id = SUPER_GROUP[list_line[0]]
        IDX=sog_id + '-' + '%04d' % sog_counter(SOG_COUNTER,sog_id)
        #Fill the dict_spp with gene for each species in the orthogroup
        for x in list_line:
            species = gene_dict[x]
            list_spp.append(species)
            dict_spp[species] = x
        #now we order the dict_spp into a final list
        ordered_list = []
        ordered_list.append(str(len(list_line)))
        ordered_list.append(str(sog_id))
        ordered_list.append(str(IDX))
        
        for x in range(0,len(spp)):
            specx = spp[x]
            if dict_spp[specx]:
                ordered_list.append(dict_spp[specx])
                webent = ['ORG',str(sog_id),str(IDX),specx,dict_spp[specx],str(NXN)]
                web_list.append(webent)
                NXN = NXN+1
            else:
                ordered_list.append('###NONE###')
                webent = ['ORG',str(sog_id),str(IDX),specx,'NONE',str(NXN)]
                web_list.append(webent)
                NXN = NXN+1
        list_spp.sort()
        Ortho_spp_freq.append(list_spp)
        #Ortho_spp_freq.append(str('\t'.join(list_spp)))
        Ortho_list_out.append(ordered_list)
        SUP_orth[sog_id].append(ordered_list)
        
    Ortho_spp_freq.sort()
    OrthoGrp_out = open(args.output_prefix + "_tabulated.tab", "w")
    OrthoGrp_FREQ_out = open(args.output_prefix + "_spp_set_frequency.tab", "w")
    OrthoGrp_sort_out = open(args.output_prefix + "_tabulated_sorted.tab", "w")
    OrthoGrp_web_out = open(args.output_prefix + "_web_sorted.tab", "w")
    for rec in web_list:
        OrthoGrp_web_out.write('\t'.join(rec) + "\n")

    for x in Ortho_spp_freq:
        OrthoGrp_FREQ_out.write('\t'.join(x) + "\n")

    # for k,g in groupby(Ortho_spp_freq):
        # OrthoGrp_FREQ_out.write(str(k) + "\t" + str(len(list(g))) + "\n")

    for orthogroup in Ortho_list_out:
        OrthoGrp_out.write("\t".join(orthogroup) + "\n")
        
    for f in SUPERGROUP_LIST:
        #start outputting for each SOG first SUP_orth, then SUP_NON_ORTH_OUTS['dups'] and SUP_NON_ORTH_OUTS['births']
        dup_out = SUP_NON_ORTH_OUTS['dups'][f]
        bir_out = SUP_NON_ORTH_OUTS['births'][f]
        ort_out = SUP_orth[f]
        for ort in ort_out:
            OrthoGrp_sort_out.write("ORG" + "\t" + "\t".join(ort) + "\n")
        for ort in dup_out:
            OrthoGrp_sort_out.write("\t".join(ort) + "\n")
            gene_lst = ort[4:]
            for g in gene_lst:
                OrthoGrp_web_out.write('\t'.join([str(ort[0]),str(ort[2]),"#",str(ort[3]),g,str(NXN)]) + "\n")
                NXN = NXN+1
        for ort in bir_out:
            OrthoGrp_sort_out.write("\t".join(ort) + "\n")
            gene_lst = ort[4:]
            for g in gene_lst:
                OrthoGrp_web_out.write('\t'.join([str(ort[0]),str(ort[2]),"#",str(ort[3]),g,str(NXN)]) + "\n")
                NXN = NXN+1

    print("#\nTabulated Orthogroup file created for further use with gene order :" + str(args.output_prefix) + "_tabulated.tab")

    ######################################################## 
    """
    Get Blastp records for these Gene Births
    """
    if args.all_blast:
        proc_blast(args.all_blast,list_births,args.output_prefix)
    else:
        print("No blast file provided - skipping blast output scanning\n")
    
    ######################################################## 

    """

    """


    ######################################################## 
    """

    """

if __name__ == "__main__":
    main()
    
