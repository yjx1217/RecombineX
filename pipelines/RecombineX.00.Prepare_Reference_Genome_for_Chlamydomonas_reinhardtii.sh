#!/bin/bash
set -e -o pipefail

#######################################
# load environment variables for Varathon
source ./../../env.sh

#######################################
# set project-specific variables
ref_genome_prefix="Chlamydomonas_reinhardtii" # The file name prefix of the reference genome. Default = "Chlamydomonas_reinhardtii". 
ref_genome_download_URL="ftp://ftp.ensemblgenomes.org/pub/plants/release-49/fasta/chlamydomonas_reinhardtii/dna/Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.dna_sm.toplevel.fa.gz" # The URL for downloading the reference genome. Default = "ftp://ftp.ensemblgenomes.org/pub/plants/release-49/fasta/chlamydomonas_reinhardtii/dna/Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.dna_sm.toplevel.fa.gz".
chr_list="./../../data/Chlamydomonas_reinhardtii.chr_list.txt" # The single-column list defining chromosomes/scaffolds/contigs to be included. Default = ./../../data/Chlamydomonas_reinhardtii.chr_list.txt".
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
#######################################




#######################################
# process the pipeline

download_and_extract() {
    url=$1
    echo "Downloading $url"
    if [[ $url =~ \.gz$ ]]; 
    #if [[ $url =~ \.fa.gz$ || $url =~ \.fasta.gz$ ]]; 
    then
	download_location="$ref_genome_prefix.raw.fa.gz"
        extract_command="gunzip"
	wget -c --no-check-certificate $url -O $download_location
	gunzip $download_location
    else
	download_location="$ref_genome_prefix.raw.fa"
	wget -c --no-check-certificate $url -O $download_location
    fi
}

echo ""
echo "Retrieve the sample reference genome assembly ..."
download_and_extract $ref_genome_download_URL 
echo ""
echo "Tidy the sample reference genome assembly ..."
$RECOMBINEX_HOME/scripts/tidy_fasta.pl -i $ref_genome_prefix.raw.fa -o $ref_genome_prefix.tidy.fa
sed -i "s/>/>chr/gi" $ref_genome_prefix.tidy.fa
$RECOMBINEX_HOME/scripts/select_fasta_by_list.pl -i $ref_genome_prefix.tidy.fa -l $chr_list -m normal -o $ref_genome_prefix.tidy.lite.fa.gz

if [[ $debug = "no" ]]
then
    echo ""
    echo "Removing intermediate files ..."
    rm $ref_genome_prefix.raw.fa
    rm $ref_genome_prefix.tidy.fa
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
