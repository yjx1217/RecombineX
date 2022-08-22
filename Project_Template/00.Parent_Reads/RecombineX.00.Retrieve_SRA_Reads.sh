#!/bin/bash
set -e -o pipefail

#######################################
# load environment variables for RecombineX
source ./../../env.sh

#######################################
# set project-specific variables
sra_run_list="sample2reads_map.txt" # A simple tab separated file with two columns, in which the first column contains the sample name and the sencond column contains the corresponding SRR id. Lines started with "#" will be ignored. Default = "sample2reads_map.txt".  
#######################################





#######################################
# process the pipeline

while read -r line
do
    [[ $line == \#* ]] && continue
    [[ $line == "" ]] && continue
    IFS=$'\t' read -r sample_id srr_id <<<"$line"
    echo "retrieve reads by the SRR_id: $srr_id for the sample $sample_id ..."
    $sra_dir/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --defline-qual '+$sn[_$rn]/$ri' --gzip --split-files --skip-technical --dumpbase --read-filter pass --clip $srr_id
    mv ${srr_id}_pass_1.fastq.gz $sample_id.R1.fq.gz
    mv ${srr_id}_pass_2.fastq.gz $sample_id.R2.fq.gz
done < $sra_run_list

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
