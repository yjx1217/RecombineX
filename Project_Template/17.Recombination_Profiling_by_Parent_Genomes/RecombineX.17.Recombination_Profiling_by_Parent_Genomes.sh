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
net_quality_cutoff=30 # The net quality difference cutoff used in tetrad genotyping. Default = "30".
color_scheme="$RECOMBINEX_HOME/data/Saccharomyces_cerevisiae.color_scheme.txt" # The color scheme to use for plotting genotypes. Default = "$RECOMBINEX_HOME/data/Saccharomyces_cerevisiae.color_scheme.txt".
plot_individual_recombination_event="yes" # Whether to plot individual recombination event, "yes" by default. Default = "yes".
flanking="4000" # The recombination event flanking region (bp) for plotting. Default = "4000".
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
###########################################

# process the pipeline
###########################################
# Normally no need to change the following parameters
genome_dir="./../11.Parent_Genome_Preprocessing" # The relative path to the 11.Parent_Genome_Preprocessing directory.
genotype_dir="./../16.Tetrad_Genotyping_by_Parent_Genomes" # The relative path to the 16.Tetrad_Genotyping_by_Parent_Genomes directory.
output_dir=$batch_id # output directory to create within this current directory
min_marker_number=1 # The minimal number of markers to be considered for trustful linkage blocks. Default = "1".
min_block_size=5 # The minimal marker-bounded block size (bp) to be considered for trustful linkage blocks. Default = "5" (i.e. 5 bp).
###########################################

# filtering tetrads by spore mapping depth
#perl $RECOMBINEX_HOME/scripts/filter_tetrads_by_sequencing_depth.pl  -i $master_sample_table -o $master_sample_table.filtered \
#    -d $spore_mapping_depth -c $mapping_depth_cutoff  

perl $RECOMBINEX_HOME/scripts/batch_recombination_profiling_by_parent_genomes.pl -s $master_sample_table -n $min_marker_number -l $min_block_size \
    -d $merging_range -q $net_quality_cutoff -b $batch_id -g $genotype_dir -o $output_dir

if [[ $plot_individual_recombination_event == "yes" ]]
then
    perl $RECOMBINEX_HOME/scripts/batch_recombination_plotting_by_parent_genomes.pl \
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
