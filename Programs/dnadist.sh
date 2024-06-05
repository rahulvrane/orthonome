#!/bin/bash

i=$1
DIR=$2

$DIR/phylip-3.67/exe/dnadist << END1
gf$i.codon
I
2
Y
END1
echo " "