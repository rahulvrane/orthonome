# README #

Please cite: Rane, RV, Oakeshott, JG, Nguyen, T., Hoffmann, AA, & Lee, SF (2017). Ortho-a new pipeline for predicting high quality orthologous gene sets applicable to complete and draft genomes. BMC genomics , 18 (1), 673. 

https://doi.org/10.1186/s12864-017-4079-6

For searchable Drosophila orthologue predictions for 20 Flybase + modENCODE genomes please visit http://www.orthonome.com 


### Who do I talk to? ###

* orthonome@gmail.com

##########################################################################

# Installation #



Orthonome depends upon the following binaries to be installed and accessible in your PATH. 
1. MCL (micans.org/mcl/)

2. MAFFT (http://mafft.cbrc.jp/alignment/software/linux.html)
3. Ssearch36 (FASTA toolkit: http://faculty.virginia.edu/wrpearson/fasta/CURRENT/)
4. NCBI BLAST++ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
5. Python 2.7
6. FastTree2 (compiled using -DUSE_DOUBLE, http://www.microbesonline.org/fasttree/)
7. OpenMPI 
8. Gffread (https://github.com/gpertea/gffread)

Orthonome also requires the following Python libraries

1. Numpy
1. Pandas

To get started, run the following commands in the directory where you want to install Orthonome. 

## Clone from github ##
```
git clone https://github.com/rahulvrane/orthonome.git
cd orthonome
```
## OR Get the latest release ##
```
wget www.orthonome.com/help/orthonome_current_release.zip
unzip orthonome_current_release.zip
cd orthonome
```
## Make all scripts executable ##
```
chmod a+x *.sh *.py ./*/*sh ./*/*.py
cd Programs/
make
cd ./phylip-3.67/
make install

```
You are now ready to run the program
##########################################################################


##########################################################################


# Inputs files and preparing run inputs #

## Required inputs: ##
For every species to be analysed, the pipeline requires three input files starting with the same prefix and with the indicated suffixes

1. Gene nucleotide sequences (PREFIX.cds.fasta)

1. Gene protein sequences (PREFIX.pep.fasta)

1. Gene annotations in GFF3 format – (PREFIX.gff3)

Additionally, the pipeline also needs a file listing all prefixes named ```run_DB```


##If starting from a genome assembly and annotation in GFF3 format (NCBI format)

The genome and gff3 formats can be converted into the above inputs. This requires the GFFREAD utility listed above.
```
#Extract all cds sequences
gffread PRF_ncbi.gff -g PRF_genome.fasta -x PRF.cds

#Extract one transcript per gene if multiple isoforms have been annotated
grep ">" PRF.cds |tr ' ' '\t'|sort -k2,2 -u|cut -f1|tr -d ">" > PRF.idx

#Convert the genomic coordinates and cds sequences to the required formats (see below for detail)
python orthonome_inputs_after_gffread.py PRF PRF_ncbi.gff
```

In some cases - if the GFF has a non-NCBI-like format - the script `orthonome_inputs_after_gffread.py` might fail. For such situations - we recommend the following steps:

```
CODE=PRF
gffread $gff -g $genome -x "$CODE".cds -y "$CODE".prot 2>/dev/null && echo "GFF converted"
grep ">" "$CODE".prot|tr ' ' '\t'|sort -k2,2 -u|tr -d ">" > $CODE.selected_genes.IDX
seqtk subseq "$CODE".prot $CODE.selected_genes.IDX|sed 's/>.*=/>/g' > "$CODE".pep && echo "PEP created"
seqtk subseq "$CODE".cds $CODE.selected_genes.IDX|sed 's/>.*=/>/g' > "$CODE".nuc && echo "NUC created"
awk '$3=="gene"' $gff|sed 's/ID=.*Name=//g'|cut -f1 -d ';'|awk '{FS=OFS="\t"}{print $9,$9,$1,$7,$4}' > $CODE.preinfo
python modify_info.py $CODE.preinfo > $CODE.info

#Check that the number of genes is the same in all files.. Sometimes depending on the format - $CODE.info may have more genes if ncRNA's are included - which is OK to ignore since we focus on protein coding genes.
echo `grep -c ">" $CODE.pep` `grep -c ">" $CODE.nuc` `cat $CODE.info|wc -l`

```

The new `$CODE.{pep,nuc,info}` are now ready for use.

##########################################################################

# The entire pipeline is designed to run with resume function in the form of a wrapper script: #

```
$ORTHONOME_HOME/Orthonome_run.sh -p ORTHONOME_HOME -t 12 -f run_DB


Usage: Orthonome_run.sh [-p|--on_dir path_to_orthonome] [-t|--threads #threads] [-f|--prefix_file file_with_species_prefixes]
```
