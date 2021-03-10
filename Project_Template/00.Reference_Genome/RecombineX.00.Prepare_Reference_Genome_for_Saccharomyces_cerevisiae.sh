#!/bin/bash
set -e -o pipefail

#######################################
# load environment variables for RecombineX
source ./../../env.sh

#######################################
# set project-specific variables

debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".

#######################################



######################################
# process the pipeline
echo "retrieve sample reference genome data ..."
wget -c https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-2-1_20150113.tgz
tar -xvzf S288C_reference_genome_R64-2-1_20150113.tgz
cp ./S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113.fsa  SGDref.genome.raw.fa
cp ./S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.gff SGDref.all_feature.gff
perl $RECOMBINEX_HOME/scripts/tidy_SGDref_genome.pl -i SGDref.genome.raw.fa -o SGDref.genome.tidy.fa
perl $RECOMBINEX_HOME/scripts/select_fasta_by_list.pl -i SGDref.genome.tidy.fa -l $RECOMBINEX_HOME/data/Saccharomyces_cerevisiae.chr_list.txt -o SGDref.genome.fa -m normal
gzip SGDref.genome.fa
perl $RECOMBINEX_HOME/scripts/filter_gff_by_feature.pl -i SGDref.all_feature.gff -o SGDref.centromere.gff -f centromere -m keep

# echo "retrieve sample subtelomere GFF files ..."
# cp $RECOMBINEX_HOME/data/Saccharomyces_cerevisiae_subtelomere_gff3/SGDref.subtelomere.gff .

if [[ $debug = "no" ]]
then
    echo ""
    echo "removing intermediate files and directories ..."

    rm -rf S288C_reference_genome_R64-2-1_20150113*
    rm SGDref.genome.raw.fa
    rm SGDref.genome.tidy.fa
    rm SGDref.all_feature.gff
fi

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
