#!/bin/bash

if [ $# == 0 ]; then
        echo "Usage: Orthonome_run.sh [-p|--on_dir path_to_orthonome] [-t|--threads #threads] [-f|--prefix_file file_with_species_prefixes]"
        exit 1
fi
# Update USAGE (USAGE1, USAGE2, USAGE3 may remain unchanged):
USAGE='Usage: Orthonome_run.sh [-p|--on_dir path_to_orthonome] [-t|--threads #threads] [-f|--prefix_file file_with_species_prefixes] '
USAGE1='
Ambiguously abbreviated long option:'
USAGE2='
No such option:'
USAGE3='
Missing argument for'

# List all long options here (including leading --):
LONGOPTS=(--prefix_file --on_dir --phylogeny)

# List all short options that take an option argument here
# (without separator, without leading -):
SHORTARGOPTS=ftp

while [[ $# -ne 0 ]] ; do
  # This part remains unchanged
  case $1 in
  --) shift ; break ;;  ### no more options
  -)  break ;;          ### no more options
  -*) ARG=$1 ; shift ;;
  *)  break ;;          ### no more options
  esac

  # This part remains unchanged
  case $ARG in
  --*)
    FOUND=0
    for I in "${LONGOPTS[@]}" ; do
      case $I in
      "$ARG")  FOUND=1 ; OPT=$I ; break ;;
      "$ARG"*) (( FOUND++ )) ; OPT=$I ;;
      esac
    done
    case $FOUND in
    0) echo -e "$USAGE$USAGE2 $ARG" 1>&2 ; exit 1 ;;
    1) ;;
    *) echo -e "$USAGE$USAGE1 $ARG" 1>&2 ; exit 1 ;;
    esac ;;
  -["$SHORTARGOPTS"]?*)
    OPT=${ARG:0:2}
    set dummy "${ARG:2}" "$@"
    shift ;;
  -?-*)
    echo -e "$USAGE" 1>&2 ; exit 1 ;;
  -??*)
    OPT=${ARG:0:2}
    set dummy -"${ARG:2}" "$@"
    shift ;;
  -?)
    OPT=$ARG ;;
  *)
    echo "OOPS, this can't happen" 1>&2 ; exit 1 ;;
  esac

  # Give both short and long form here.
  # Note: If the option takes an option argument, it it found in $1.
  # Copy the argument somewhere and shift afterwards!
  case $OPT in
  -f|--prefix_file)  [[ $# -eq 0 ]] && { echo "$USAGE$USAGE3 $OPT" 1>&2 ; exit 1 ; }
              FILE=$1 ; shift ;;
  -p|--on_dir)  [[ $# -eq 0 ]] && { echo "$USAGE$USAGE3 $OPT" 1>&2 ; exit 1 ; }
              _DIR=$1 ; shift ;;
  -t|--threads) [[ $# -eq 0 ]] && { echo "$USAGE$USAGE3 $OPT" 1>&2 ; exit 1 ; }
              THR=$1 ; shift ;;
  *)          echo "$USAGE$USAGE2 $OPT" 1>&2 ; exit 1 ;;
  esac
done
#############################################################################
RED='\033[0;31m'
NC='\033[0m'
GREEN='\033[0;32m'


export WKDIR=`pwd`
export _DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "${GREEN}[`date`]${NC}: Using Orthonome installation in: $_DIR"

export NUMJOBS=$((THR/3))

echo "${GREEN}[`date`]${NC}: Number of threads: $THR"
echo "${GREEN}[`date`]${NC}: Using 3 threads per job, number of concurrent jobs run in parallel: $NUMJOBS"

# Remaining arguments are now in "$@":
echo "${GREEN}[`date`]${NC}: Testing dependencies"

foo="mafft";command -v $foo >/dev/null 2>&1 || { echo >&2 "I require '$foo' but it's not installed.  Aborting."; exit 1; }
foo="gffread";command -v $foo >/dev/null 2>&1 || { echo >&2 "I require '$foo' but it's not installed.  Aborting."; exit 1; }
foo="mcl";command -v $foo >/dev/null 2>&1 || { echo >&2 "I require '$foo' but it's not installed.  Aborting."; exit 1; }
foo="makeblastdb";command -v $foo >/dev/null 2>&1 || { echo >&2 "I require '$foo' but it's not installed.  Aborting."; exit 1; }
foo="blastp";command -v $foo >/dev/null 2>&1 || { echo >&2 "I require '$foo' but it's not installed.  Aborting."; exit 1; }
foo="ssearch36";command -v $foo >/dev/null 2>&1 || { echo >&2 "I require '$foo' but it's not installed.  Aborting."; exit 1; }
foo="parallel";command -v $foo >/dev/null 2>&1 || { echo >&2 "I require '$foo' but it's not installed.  Aborting."; exit 1; }
#$foo="mafft";command -v $foo >/dev/null 2>&1 || { echo >&2 "I require $foo but it's not installed.  Aborting."; exit 1; }


################################################################################
echo "${GREEN}[`date`]${NC}: Inspecting inputs"

for f in *.preinfo; do 
    $_DIR/modify_info.py $f > ${f%%.preinfo}.info
done

#

for pr in `cat $FILE`; do 
    for suf in pep nuc info; do
        if [[-s $pr.$suf]];
            echo "${GREEN}[`date`]${NC}: $pr.$suf OK"
        else:
            echo "${RED}[`date`]${NC}: ERROR: $pr.$suf is empty. Please check input"
            exit 1
        fi
    done
done



################################################################################
echo "${GREEN}[`date`]${NC}: Generating Blast databases" 
mkdir -p $WKDIR/blast_pairs
cd $WKDIR/blast_pairs
parallel --gnu -j$THR 'makeblastdb -in ../{}.pep -dbtype prot -out {}' :::: $WKDIR/$FILE

echo "${GREEN}[`date`]${NC}: Running blastp runs" 

ln -s $_DIR/blast_runner.sh .

parallel --gnu -j1 --bar 'bash blast_runner.sh {1} {2} {3} {4}' :::: $WKDIR/$FILE :::: $WKDIR/$FILE ::: $THR ::: `pwd`

for x in `cat $WKDIR/$FILE`; do 
    for y in `cat $WKDIR/$FILE`; do
        #$_DIR/blast_runner.sh $x $y $THR `pwd`
        if [ ! -f $x"_"$y".success" ];
        then
            echo "${RED}[`date`]${NC}: ERROR: Blast run between $x and $y failed"
            exit 1
        fi
    done
done

################################################################################

echo "${GREEN}[`date`]${NC}: Clustering Blast comparisons between species" 

cat *.blastp|grep -v "#"|parallel --gnu --pipe -q awk '{OFS="\t"}{if ($11<=0.5) print $1, $2, $12;else print $1,$2,0}' | mcl - --abc -q x -V all -te $THR -o Allruns.clusters

################################################################################
echo "${GREEN}[`date`]${NC}: Finished running Blast comparisons between species" 

mkdir -p $WKDIR/sw_scores
cd $WKDIR/sw_scores

################################################################################
echo "${GREEN}[`date`]${NC}: Running Smith-Waterman alignments between species" 


ln -s $_DIR/sw_runner.sh .

parallel --gnu -j1 --bar 'bash sw_runner.sh {1} {2} {3} {4}' :::: $WKDIR/$FILE :::: $WKDIR/$FILE ::: $THR ::: `pwd`

for x in `cat $WKDIR/$FILE`; do 
    for y in `cat $WKDIR/$FILE`; do
        #$_DIR/sw_runner.sh $x $y $THR `pwd`
        if [ ! -f $x"_"$y".success" ];
        then
            echo "${RED}[`date`]${NC}: ERROR: Smith-Waterman run between $x and $y failed"
            exit 1
        fi
    done
done

echo "${GREEN}[`date`]${NC}: Finished running Smith-Waterman alignments between species" 

cd $WKDIR

################################################################################
echo "${GREEN}[`date`]${NC}: Preparing batch scripts for each pairwise comparison" 

for f in `cat $WKDIR/$FILE`; do 
    echo -e $f"\t"`cut -f1 ../$f.info|tr '\n' '\t'`|tr ' ' '\t'
done > genelists.txt


awk '{print $0"_"NR}' $FILE > Spp_list.idx
python "$_DIR"/combinations.py Spp_list.idx | tr -d '() ' > combinations.txt

N=1
for f in `cat combinations.txt`; do 
    P=`echo $f | cut -f1 -d,`
    S=`echo $f | cut -f2 -d,`
    SPP1=`awk -v var="$P" 'NR == var' $FILE`
    SPP2=`awk -v var="$S" 'NR == var' $FILE`
    #rm -fr $SPP1"_"$SPP2
    echo $SPP1"_"$SPP2 >> combination_prf.txt
    echo "mkdir -p "$SPP1"_"$SPP2"; cd "$SPP1"_"$SPP2"; bash "$_DIR"/Process_Species_pairs.sh "$WKDIR"/"$SPP1".pep "$WKDIR"/"$SPP2".pep "$WKDIR"/"$SPP1".nuc "$WKDIR"/"$SPP2".nuc "$WKDIR"/"$SPP1".info "$WKDIR"/"$SPP2".info "$SPP1","$SPP2" "$P" "$S" 3 all 1> PairComparison.log 2>&1 && touch "$WKDIR"/PairComparison"$N".success" > PairComparison"$N".sh
    #echo Batch file for $SPP1 vs $SPP2 comparison created
    N=$((N+1))
done

################################################################################
mkdir -p $WKDIR/MultiMSOAR_inputs

echo "${GREEN}[`date`]${NC}: Using 3 threads per job, number of concurrent jobs run in parallel: $NUMJOBS"

################################################################################
parallel --gnu --bar -j$NUMJOBS 'bash {}' ::: PairComparison*.sh

ERR=0
NOS=0
for prcomb in `cat combination_prf.txt`; do 
    if [ ! -s $prcomb/MSOAR2_result ];
    then
        echo "${RED}[`date`]${NC}: ERROR: Orthonome pairwise run $prcomb failed. Check $prcomb/PairComparison.log for details"
        ERR=$(($ERR+1))
    else
        echo "${GREEN}[`date`]${NC}: Orthonome pairwise run $prcomb successfully completed."
        NOS=$(($NOS+1))
    fi
done

if [ "$ERR" ! = 0 ];then
    echo "${RED}[`date`]${NC}: Errors found in earlier runs. See log above"
else
    echo "${GREEN}[`date`]${NC}: Orthonome pairwise runs completed and moving on to multispecies stages."
fi

#

################################################################################
cd $WKDIR/MultiMSOAR_inputs

echo "${GREEN}[`date`]${NC}: Generating consensus phylogeny using pairwise orthologues"

$_DIR/Ortholog_pairs_to_FastTreephy.py -l $WKDIR/genelists.txt -g $WKDIR/MultiMSOAR_inputs/ -p $WKDIR -n $WKDIR -t $THR 

################################################################################
if [ ! -s CONCAT_align_nuc.nwk ];then
    echo "${RED}[`date`]${NC}: $_DIR/Ortholog_pairs_to_FastTreephy.py not run successfully. See log above"
    exit 1
fi

cp CONCAT_align_nuc.nwk Tree
################################################################################

echo "${GREEN}[`date`]${NC}: Modifying phylogeny for MultiMSOAR"

cat $WKDIR/$FILE|awk '{print $1 "S"$NR-1}'|while read line; do S=$(echo $line|cut -f1 -d ' '); PS=$(echo $line|cut -f2 -d ' '); sed -i "s:$S:$PS:g" Tree;done

echo "${GREEN}[`date`]${NC}: Running MultiMSOAR"

################################################################################
$_DIR/MultiMSOAR `wc -l $WKDIR/$FILE||awk '{print $1}'` Tree $WKDIR/blast_pairs/Allruns.clusters Geneinfo Orthogroups 

echo "${GREEN}[`date`]${NC}: Creating final Orthonome output"
$_DIR/summarise_orthogroups_internet_OUT.py -l genelists.txt -m $WKDIR/blast_pairs/Allruns.clusters -i GeneInfo -g OrthoGroups -s $WKDIR/Spp_list.idx -o Orthonome_out

echo "${GREEN}[`date`]${NC}: All done!"
