#!/bin/bash

if [ $# != 4 ]; then
	echo "<Usage>: blast_runner.sh Spp1 Spp2 THREADS OUTDIR"
	exit
fi

FA="$1"
FB="$2"
THR="$3"
OUTDIR="$4"

if [ ! -s $OUTDIR/"$FA"_"$FB".blastp ];
then
    cat $FA.pep | parallel --gnu -j"$THR" --recstart '>' --block 100k --pipe blastp -db $FB -query - -outfmt 7 -evalue 0.5 > $OUTDIR/"$FA"_"$FB".blastp 2> /dev/null && echo Blast Complete $FA vs $FB && touch "$FA"_"$FB".success
else
    echo Blast Complete $FA vs $FB && touch "$FA"_"$FB".success
fi
