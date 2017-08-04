#!/bin/bash

if [ $# != 4 ]; then
	echo "<Usage>: sw_runner.sh Spp1 Spp2 THREADS OUTDIR"
	exit
fi

FA="$1"
FB="$2"
THR="$3"
OUTDIR="$4"

if [ ! -s $OUTDIR/"$FA"_"$FB".SWscore ];
then
    ssearch36 -T 32 -d 0 -m 9 -z -1 $FA.pep $FB.pep |tail -n +9|egrep -v '(Algorithm|Parameters|Library|residues|Query|<)'|sed '/^$/d'|sed '/>/ s/^.*>>>/>/g;/>/ s/ .*//g;s/(.*)//g;s/The best scores are:/#TopScores/g;s/ \+/\t/g' > $OUTDIR/"$FA"_"$FB".SWscore 2> /dev/null && echo SW search Complete $FA vs $FB && touch "$FA"_"$FB".success
else
    echo SW search Complete $FA vs $FB && touch "$FA"_"$FB".success
fi
