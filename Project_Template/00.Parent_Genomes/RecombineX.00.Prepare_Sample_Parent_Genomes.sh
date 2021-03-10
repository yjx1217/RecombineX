#!/bin/bash
set -e -o pipefail

#######################################
# load environment variables for RecombineX
source ./../../env.sh

#######################################
# set project-specific variables

# none

#######################################
# process the pipeline

echo "retrieve sample parental genome data ..."
for i in S288C SK1
do
    cp $RECOMBINEX_HOME/data/$i.genome.fa .
    cp $RECOMBINEX_HOME/data/$i.all_feature.gff .
    perl $RECOMBINEX_HOME/scripts/filter_gff_by_feature.pl -i $i.all_feature.gff -o $i.centromere.gff -f centromere -m keep
done

# echo "retrieve sample subtelomere GFF files ..."
# for i in S288C SK1
# do
#     cp $RECOMBINEX_HOME/data/Saccharomyces_cerevisiae_subtelomere_gff3/$i.subtelomere.gff .
# done

echo ""
echo "removing intermediate files and directories ..."
for i in S288C SK1
do
    rm $i.all_feature.gff
done


############################
# checking bash exit status
if [[ $? -eq 0 ]]
then
    echo ""
    echo "RecombineX message: This bash script has been successfully processed! :)"
    echo ""
    echo ""
    exit 0
fi
############################
