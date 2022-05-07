#!/bin/bash
set -e -o pipefail

##########################################
# load environment variables for RecombineX
source ./../../env.sh

###########################################
# set project-specific variables
batch_id="Batch_S288C-SK1" # The batch id for this analysis. This batch id will also be used as the output directory name. Default = "Batch_S288C-SK1".
master_sample_table="Master_Sample_Table.${batch_id}.txt" # The master sample table for this batch. Default = "Master_Sample_Table.${batch_id}.txt".
window_size=250 # The window size for sliding-window-based CNV profiling. This setting should be kept consistent with the same setting in ./../01.Reference_Genome_Preprocessing. Default = "250" (i.e. 250 bp).
threads=4 # The number of threads to use. Default = "4".
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
###########################################


# process the pipeline
###########################################
# Normally no need to change the following parameters
reference_genome_assembly_dir="./../01.Reference_Genome_Preprocessing" # The relative path to 01.Reference_Genome_Preprocessing. Default = "./../01.Reference_Genome_Preprocessing".
gamete_reads_dir="./../00.Gamete_Reads" # The relative path to 00.Gamete_Reads. Default = "./../00.Gamete_Reads".
output_dir="$batch_id" # The output directory. Default = "$batch_id".
mapping_quality_cutoff_for_mpileup=30 # The mapping quality cutoff for filtering the resulting bam file. Default = "30".
min_mappability=0.85 # The minimal mappability for sliding-window-based CNV profiling. Default = "0.85".
excluded_chr_list_for_cnv_profiling="" # The relative path to the list for specifying chromosomes/scaffolds/contigs in relabeled parental genomes (in ./../01.Reference_Genome_Preprocessing) to be exclued for CNV profiling. Default = "".
##################################

perl $RECOMBINEX_HOME/scripts/batch_read_mapping_to_reference_genome.pl \
    -i $master_sample_table \
    -t $threads \
    -q $mapping_quality_cutoff_for_mpileup \
    -reference_genome_assembly_dir $reference_genome_assembly_dir \
    -gamete_reads_dir $gamete_reads_dir \
    -output_dir $output_dir \
    -min_mappability $min_mappability \
    -window $window_size \
    -step $window_size \
    -ploidy 1 \
    -excluded_chr_list_for_cnv_profiling $excluded_chr_list_for_cnv_profiling \
    -debug $debug

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

