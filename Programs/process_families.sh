#!/bin/bash

export DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

PDIR=`pwd`
i=$1

####################################
# Create MultipleSequenceAlignment of all family members
#####################################

mafft --quiet Families/gf$i.pep |tr "." "*" > Families/gf$i.fasta

#####################################
# Generate Codon alignment : Use trancation sensitive version of RevTrans
#####################################

#$DIR/RevTrans-1.4/revtrans_jarmo.py Families/gf$i.nuc Families/gf$i.fasta -readthroughstop -match name -O phylip > Families/gf$i.codon
$DIR/pal2nal Families/gf$i.fasta Families/gf$i.nuc -output paml > Families/gf$i.codon

mkdir $i
cd $i
#echo moved to new DIR
#####################################
# 1. Calculate DNA distance matrix (dnadist)
# 2. Generate probable inparalog sets
#####################################

$DIR/mapToPhylip $PDIR/Families/gf$i.codon gf$i.codon gf$i.map
#echo phylip mapped

$DIR/dnadist.sh $i $DIR 1> /dev/null

#echo dnadist run
mv outfile $PDIR/gf$i.dist
mv gf$i.map $PDIR
rm *
cd ..
rm -fr $i
if [[ ! -s gf$i.dist ]];then
    echo "family processing failed for family number "$1"; check alignment input"
    rm -f gf$i.codon
    touch gf$i.tags
    mv gf$i.dist Families/
    mv gf$i.map Families/
    mv gf$i.tags Families/
    exit
fi
#$DIR/addOutgroup gf$i.distx > gf$i.dist
$DIR/addOutgroup gf$i.dist | fnj -I phylip -O newick -m BIONJ|$DIR/treeTran 0000000000 /dev/stdin |tee gf$i.tree| $DIR/tagGenerator /dev/stdin Families/gf$i.pos gf$i.map gf$i.tags
rm -f gf$i.codon

#####################################
# Move files
#####################################
mv gf$i.dist Families/
mv gf$i.map Families/
mv gf$i.tree Families/
mv gf$i.tags Families/
