#  Process_Species_pairs.sh
#  
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
 
#!/bin/bash

if [ $# != 11 ]; then
	echo "<Usage>: Process_Species_pairs.sh G1.pep G2.pep G1.nuc G2.nuc G1.info G2.info Spp1,Spp2 Spp_1_index Spp_2_index THREADS [msprep?]"
	exit
fi
export _DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

check_file () {
    FILE=$1
    if [[ ! -s $1 ]];then
    echo "$1 file empty"
    exit
    else 
    echo "$1 file OK"
    fi
}


echo "$7"
THRDS=${10}
PRF1=`echo $7|cut -f1 -d ','`
PRF2=`echo $7|cut -f2 -d ','`
G1G2="$PRF1"_"$PRF2"
G2G1="$PRF2"_"$PRF1"
G1G1="$PRF1"_"$PRF1"
G2G2="$PRF2"_"$PRF2"

echo "Input files $1 $2 $3 $4 $5 $6"

echo Current Directory is `pwd`
echo Working directory is "$WKDIR"
#################################################################
## Prepare Inputs
#################################################################
check_file $1
check_file $2
check_file $3
check_file $4
check_file $5
check_file $6
check_file ../blast_pairs/"$G1G2".blastp
check_file ../sw_scores/"$G1G2".SWscore
check_file ../blast_pairs/"$G2G1".blastp
check_file ../sw_scores/"$G2G1".SWscore
check_file ../blast_pairs/"$G1G1".blastp
check_file ../blast_pairs/"$G2G2".blastp


ln -fs ../blast_pairs/"$G1G2".blastp ../blast_pairs/"$G2G1".blastp .
ln -fs ../sw_scores/"$G1G2".SWscore ../sw_scores/"$G2G1".SWscore .
ln -fs ../blast_pairs/"$G1G1".blastp ../blast_pairs/"$G2G2".blastp .
ln -fs $5 G1.info
ln -fs $6 G2.info


#####################################
# MCL Clustering
#####################################
if [[ -e G1G2.cluster ]]; then
    echo Clustering successfully completed before! YAY! Making soft link
else
    cat *.blastp|grep -v "#"|sort -k1,1 |parallel --gnu --pipe -q awk '{OFS="\t"}{if ($11<=0.5) print $1, $2, $12;else print $1,$2,0}' | mcl - --abc -te "$THRDS" -I 2.0 -o G1G2.cluster
fi
echo [`date`]:"Clustering complete"
#################################################################
##  Process Families identified using MCL clusters
#################################################################

##################################################
# Separate Different Gene Families
##################################################
rm -Rf Families/
mkdir Families/
if [[ ! -e "$5"x ]];then 
    cut -f1 $5 |sed 's/^/>/'|sed -r '/>/s/$/\nDUMY/g' > "$5"x && \
    echo [`date`]:"Faux_index created "$5"x"
else
    echo [`date`]:"Faux_index found "$5"x"
fi

if [[ ! -e "$6"x ]];then 
    cut -f1 $6 |sed 's/^/>/'|sed -r '/>/s/$/\nDUMY/g' > "$6"x && \
    echo [`date`]:"Faux_index created "$6"x"
else
    echo [`date`]:"Faux_index found "$6"x"
fi

echo [`date`]:"Faux_indices catalogued"


if [[ ! -e family_process.success ]];then
NumOfFamilies=`$_DIR/Programs/familyAssembler G1G2.cluster $1 $2 $3 $4 $5"x" $6"x"`

echo [`date`]:"Family data separated for $NumOfFamilies families"

######################################################################
# Construct Gene Tree for each Family and Find positional inparalogs
######################################################################

parallel --gnu -j"$THRDS" $_DIR/Programs/process_families.sh {} ::: `seq $NumOfFamilies` 1>famile_proc.log 2>famile_proc.err

echo [`date`]:"$NumOfFamilies families processed"
rm -f TAGs 

if [ $NumOfFamilies -ne 0 ]; then
	cat Families/gf*.tags > TAGs
else
	touch TAGs
fi
touch family_process.success && rm -fr Families
else
echo [`date`]:"Family data processed"
cp Results/TAGs .
fi
#################################################################
##  Prepare and run MSOAR
#################################################################

################################################
# Remove Inparalogs in TAGs and Invoke Pair_comparisons
# Input: TAGs G1G2.blastp
# Output: msoar_result
################################################
if [[ ! -e msoar_prep.success ]];then

$_DIR/Programs/normalise_top_blast_hits.py "$G1G2".blastp "$G2G1".blastp "$G1G2".SWscore "$G2G1".SWscore 0.05 1> normalise.log 2>normalise.err

$_DIR/Programs/removeInparalogs2 TAGs G1ToG2.top NI_G1ToG2.top
$_DIR/Programs/removeInparalogs2 TAGs G2ToG1.top NI_G2ToG1.top

rm -f G1G2 G2G1 G1G2*_* G2G1*_* pre_res msoar_output

for f in G*.top;do echo $f;comm -3 <(cat $f|cut -f1|sort|uniq) <(cat NI_$f|cut -f1|sort|uniq) > ${f%%T}.inparalogs;done

$_DIR/tools/cutbyname.pl -f 1,2,4,3 NI_G1ToG2.top|sort -k1,2 -u > G1G2_bh_top
$_DIR/tools/cutbyname.pl -f 1,2,4,3 NI_G2ToG1.top|sort -k1,2 -u > G2G1_bh_top
$_DIR/tools/joinby.pl $5 1 G1G2_bh_top 1 > G1G2
$_DIR/tools/joinby.pl $6 1 G2G1_bh_top 1 > G2G1
$_DIR/tools/sortposition_whole.pl G1G2 >G1G2_pos
$_DIR/tools/sortposition_whole.pl G2G1 >G2G1_pos
$_DIR/tools/joinby.pl -o G1G2 6 G2G1_pos 2 > G1G2_s
$_DIR/tools/joinby.pl G1G2_s 2 G1G2_pos 2 >G1G2_ss
$_DIR/tools/cutbyname.pl -f 10,1,3-6,9,7-8 G1G2_ss > G1G2_whole
sort -k1n G1G2_whole > G1G2_msoar
$_DIR/tools/joinby.pl -o G2G1 6 G1G2_pos 2 > G2G1_s
$_DIR/tools/joinby.pl G2G1_s 2 G2G1_pos 2 > G2G1_ss
$_DIR/tools/cutbyname.pl -f 10,1,3-6,9,7-8 G2G1_ss > G2G1_whole
sort -k1n G2G1_whole >G2G1_msoar
sort -u G1G2_msoar>G1G2_msoar_u
sort -u G2G1_msoar>G2G1_msoar_u
sort -k1n G2G1_msoar_u >G2G1_msoar
sort -k1n G1G2_msoar_u >G1G2_msoar 
check_file G1G2_msoar 
check_file G2G1_msoar 

touch msoar_prep.success

else
echo [`date`]:"MSOAR prepped earlier."
check_file G1G2_msoar
fi


######################################################
if [[ ${11} = "msprep" ]];then
echo [`date`]: "MSOAR prep ready"
exit
fi


if [[ ! -e msoar.success ]];then
$_DIR/Programs/MSOAR/MSOAR G1G2_msoar G2G1_msoar pre_res msoar_output 1>msoar.log 2>msoar.err && touch msoar.success
$_DIR/tools/cutbyname.pl -f 1,3 msoar_output > msoar_result
#rm -f G1G2*_* G2G1*_* pre_res msoar_output
check_file msoar_result
else
echo [`date`]:"MSOAR run earlier."
cp Results/msoar_result .
check_file msoar_result
fi


if [[ ! -e Orthology_pair_classifications.success ]];then
$_DIR/Programs/normalizeScores NI_G1ToG2.top > G1ToG2.score
$_DIR/Programs/normalizeScores NI_G2ToG1.top > G2ToG1.score
if [[ ! -s MSOAR2_result ]];then
echo [`date`] "MSOAR2 run successfully."
touch Orthology_pair_classifications.success
fi
else
echo [`date`]:"MSOAR2 was run earlier... what do you want me to do? Here  - I will paste stuff now!."
fi

mkdir -p Results/
cp G1*G2* Results/
cp G2*G1* Results/
cp TAGs Results/
cp msoar_result Results/
cp MSOAR2_result Results/Ortholog_Pairs.txt

P=$8
S=$9
PX=$((P-1))
SX=$((S-1))

$_DIR/Programs/postprocessing $5 $6 msoar_result G1ToG2.score G2ToG1.score 1.0 > MSOAR2_result
rm -f MSOAR2_out.scores
$_DIR/OrthologPair_sw_scoring.py G1ToG2.SW G2ToG1.SW MSOAR2_result MSOAR2_out.scores
wc -l MSOAR2_out.scores MSOAR2_result
cp  MSOAR2_out.scores ../MultiMSOAR_inputs/S"$PX"_S"$SX"

echo -e "`wc -l MSOAR2_out.scores`\tpairs in $SPP1 vs $SPP2"

