#!/bin/bash
set -e -o pipefail

##########################################
# load environment variables for RecombineX
source ./../../env.sh

###########################################
# set project-specific variables
batch_id="Batch_S288C-SK1" # The batch id for this analysis. This batch name will also be used as the output directory name. Default = "Batch_S288C-SK1".
master_sample_table="Master_Sample_Table.${batch_id}.txt" # The master sample table for this batch. Default = "Master_Sample_Table.${batch_id}.txt".
window_size=250 # The window size for the non-overlaping sliding-window-based CNV profiling. Default = 250 (i.e. 250 bp).
threads=4 # The number of threads to use. Default = "4".
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
###########################################



# process the pipeline
###########################################
# Normally no need to change the following parameters
parent_genome_assembly_dir="./../11.Parent_Genome_Preprocessing"
gamete_reads_dir="./../00.Gamete_Reads"
output_dir="$batch_id" # The output directory 
min_mappability=0.85 # The minimal mappability for sliding-window-based CNV profiling. Default = "0.85".
mapping_quality_cutoff_for_mpileup=30 # The mapping quality cutoff for filtering the resulting bam file. Default = "30".
##########################################

perl $RECOMBINEX_HOME/scripts/batch_read_mapping_to_parent_genomes.pl \
    -i $master_sample_table \
    -t $threads \
    -q $mapping_quality_cutoff_for_mpileup \
    -parent_genome_assembly_dir $parent_genome_assembly_dir \
    -gamete_reads_dir $gamete_reads_dir \
    -output_dir $output_dir \
    -min_mappability $min_mappability \
    -window $window_size \
    -step $window_size \
    -ploidy 1

# clean up intermediate files
# if [[ $debug = "no" ]]
# then

# fi

############################
# checking bash exit status
if [[ $? -eq 0 ]]
then
    echo ""
    echo "##########################################################################"
    echo ""
    echo "RecombineX message: This bash script has been successfully processed! :)"
    echo ""
    echo "##########################################################################"
    echo ""
    exit 0
fi
############################

