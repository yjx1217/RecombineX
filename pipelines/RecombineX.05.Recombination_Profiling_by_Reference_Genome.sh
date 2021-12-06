#!/bin/bash
set -e -o pipefail

##########################################
# load environment variables for RecombineX
source ./../../env.sh

###########################################
# set project-specific variables
batch_id="Batch_S288C-SK1" # The batch id used for the gamete read mapping analysis. Default = "Batch_S288C-SK1".
master_sample_table="Master_Sample_Table.$batch_id.txt" # The master sample table for this batch. Default = "Master_Sample_Table.${batch_id}.txt".
merging_range=5000 # The distance range (bp) for merging nearby COs. Default = "5000" (i.e. 5000 bp). 
net_quality_cutoff=50 # The net quality difference cutoff used in tetrad genotyping. Default = "50".
color_scheme="$RECOMBINEX_HOME/data/Saccharomyces_cerevisiae.color_scheme.txt" # The color scheme to use for plotting genotypes. Default = "$RECOMBINEX_HOME/data/Saccharomyces_cerevisiae.color_scheme.txt".
plot_individual_recombination_event="no" # Whether to plot individual recombination event: "yes" or "no". Default = "no".
flanking="4000" # The recombination event flanking region (bp) for plotting. Default = "4000".
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
###########################################


# process the pipeline
###########################################
# Normally no need to change the following parameters
genome_dir="./../01.Reference_Genome_Preprocessing" # The relative path to the 01.Reference_Genome_Preprocessing.
genotype_dir="./../04.Gamete_Genotyping_by_Reference_Genome" # The relative path to the 04.Gamete_Genotyping_by_Reference_Genome directory.
output_dir=$batch_id # output directory to create within this current directory
min_marker_number=1 # The minimal number of markers to be considered for trustful linkage blocks. Default = "1".
min_block_size=1 # The minimal marker-bounded block size (bp) to be considered for trustful linkage blocks. Default = "1" (i.e. 1 bp).
###########################################


# filtering tetrads by spore mapping depth
#perl $RECOMBINEX_HOME/scripts/filter_tetrads_by_sequencing_depth.pl  -i $master_sample_table -o $master_sample_table.filtered \
#    -d $spore_mapping_depth -c $mapping_depth_cutoff  

perl $RECOMBINEX_HOME/scripts/batch_recombination_profiling_by_reference_genome.pl -s $master_sample_table -n $min_marker_number -l $min_block_size \
    -d $merging_range -q $net_quality_cutoff -b $batch_id -g $genotype_dir -o $output_dir

if [[ $plot_individual_recombination_event == "yes" ]]
then
    perl $RECOMBINEX_HOME/scripts/batch_recombination_plotting_by_reference_genome.pl \
	--sample_table $master_sample_table \
	--genome_dir $genome_dir \
	--genotype_dir $genotype_dir \
	--recombination_profiling_dir $(pwd) \
	--batch_id $batch_id \
	--qual_diff $net_quality_cutoff \
	--color_scheme $color_scheme \
	--plot_centromere "no" \
	--flanking $flanking
fi

# clean up intermediate files
# if [[ $debug = "no" ]]
# then

# fi

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
