from collections import defaultdict
from Bio.Alphabet import generic_nucleotide
from Bio.Seq import Seq
import re
import sys
from tqdm import tqdm

#!/usr/bin/env python3

def FASTA(filename):
    try:
        f = open(filename)
    except IOError:
        print("The file, %s, does not exist" % filename)
        return

    order = []
    sequences = {}

    for line in tqdm(f):
        if line.startswith('>'):
            name = line[1:].rstrip('\n')
            order.append(name)
            sequences[name] = ''
        else:
            sequences[name] += line.rstrip('\n')

    return sequences

def parse_gff(gffile):
    rna_to_gene = {}
    gff_feats = {}
    for line in open(gffile).readlines():
        linel = line.strip().split("\t")
        if len(linel) == 9:
            # Run for only genes to get coordinate data
            if linel[2] == "gene":
                try:
                    # characterise all descriptions to grab the 'Name'
                    dumdict = {}
                    for rec in linel[8].split(";"):
                        k, v = rec.split("=")
                        dumdict[k] = v
                    chr, start, phase = linel[0],linel[3],linel[6]
                    name = dumdict["gene"].split("%")[0]
                    gff_feats[name] = "%s\t%s\t%s" % (chr, phase, start)
                except:
                    print("failed on line: %s" % line)
            elif linel[2] == "mRNA":
                # characterise all descriptions to grab the 'Name'
                dumdict = {}
                for rec in linel[8].split(";"):
                    k, v = rec.split("=")
                    dumdict[k] = v
                rna_to_gene[dumdict["ID"]] = dumdict["gene"]
    return gff_feats, rna_to_gene
    
def ret_break(s):
    a,b = re.split("\s",s)[:2]
#    b1 = b.split("=")[1]
    return a, b

def main():
    prf = sys.argv[1]
    pgff = sys.argv[2]
    print("Processing: "+ prf)
    
    #Get the GFF into and mapping dict
    print("Processing GFF for: "+ prf)
    gff_info, rnadict = parse_gff(pgff)
    
    #Get the sequences
    print("Processing CDS for: "+ prf)
    spp_cds = FASTA(prf+".cds")
    
    print("Total seqs = %d" % len(spp_cds))
    
    #Get the list of RNA we want to keep
    idxs = set([line.strip().split(" ")[0] for line in open(prf+".idx").readlines()])
    
    #Start the selection process 
    print("Total selected seqs = %d" % len(idxs))
    print(list(idxs)[0:10])
    print(spp_cds.keys()[:10])




    spp_cds_sel = {}
    for k, v in tqdm(spp_cds.items()):
        if re.split("\t",k)[0] in idxs:
            spp_cds_sel[k] = v
    print("Total selected seqs extracted = %d" % len(spp_cds_sel))
    
    # Prepare output files: nuc, pep and preinfo
    print("outputting into files : %s.{nuc,pep,preinfo}" % prf)
    nuc = open(prf+".nuc","w")
    pep = open(prf+".pep","w")
    preinfo = open(prf+".preinfo","w")
    
    
    for k, v in tqdm(spp_cds_sel.items()):
        rna, gene = ret_break(k)
        pc = "%s_%s_%s" % (prf, gene, rna)

        if gene in gff_info.keys():
            nuc.write(">%s\n%s\n" % (pc,v))
            pep.write(">%s\n%s\n" % (pc,str(Seq(v,generic_nucleotide).translate())))
            preinfo.write("%s\n" % "\t".join([pc,pc,gff_info[gene]])) 
        else:
            print("%s not found in the valid gene list. could be a pseudogene classification. excluded from analysis" % gene)
    
    nuc.close()
    pep.close()
    preinfo.close()
#
if __name__ == "__main__":
    main()
