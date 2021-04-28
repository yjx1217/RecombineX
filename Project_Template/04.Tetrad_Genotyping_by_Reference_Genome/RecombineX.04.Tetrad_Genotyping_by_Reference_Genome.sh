#!/bin/bash
set -e -o pipefail

##########################################
# load environment variables for RecombineX
source ./../../env.sh

###########################################
# set project-specific variables
batch_id="Batch_S288C-SK1" # The batch id used for the gamete read mapping analysis. Default = "Batch_S288C-SK1" 
master_sample_table="Master_Sample_Table.${batch_id}.txt" # The master sample table for this batch. Default = "Master_Sample_Table.${batch_id}.txt".
net_quality_cutoff=20 # The net quality cutoff for genotyping. Default = "20".
apply_cnv_filter="yes" # Whether to set gamete genotype to NA for potential CNV regions in gametes. Set this option to "no" if the gamete sequencing depth is very low (e.g. <= 1). Default = "yes".
allow_heteroduplex="no" # Whether to consider the possibility of heteroduplex formation. Default = "no".
chr_list="$RECOMBINEX_HOME/data/Saccharomyces_cerevisiae.chr_list.txt" # The included chromosome list for the analyzed genome. Default = "$RECOMBINEX_HOME/data/Saccharomyces_cerevisiae.chr_list.txt".
color_scheme="$RECOMBINEX_HOME/data/Saccharomyces_cerevisiae.color_scheme.txt" # The color scheme to use for plotting genotypes. This file is a tab-delimited two column list file in which the first column is parent_id and the second column is the hex color code. Default = "$RECOMBINEX_HOME/data/Saccharomyces_cerevisiae.color_scheme.txt". 
plot_centromere="yes" # Whether to plot centromere in the generated genotyping plots. Please note that enable this option requires that you have the ref.centromere.relabel.gff file ready in the "./../01.Reference_Genome_Preprocessing" directory. Default = "yes". 
same_cross_combination_for_the_batch="yes" # Wether all the samples in the current batch come from the same cross combination (i.e. shared the same parents): "yes" or "no". When "yes", RecombineX will automatically profle and plot the parental allele frequency for every markers. Default = "yes".
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
###########################################


# process the pipeline
############################################
# No need to change the following parameters in general
genome_dir="./../01.Reference_Genome_Preprocessing" # The relative path to the "01.Reference_Genome_Preprocessing directory. Default = "./../01.Reference_Genome_Preprocessing".
marker_dir="./../02.Polymorphic_Markers_by_Reference_based_Read_Mapping" # The relative path to the 02.Polymorphic_Markers_by_Reference_based_Read_Mapping directory. Default = "./../02.Polymorphic_Markers_by_Reference_based_Read_Mapping".
gamete_read_mapping_dir="./../03.Gamete_Read_Mapping_to_Reference_Genome" # The relative path to the 03.Gamete_Read_Mapping_to_Reference_Genome directory. Default = "./../03.Gamete_Read_Mapping_to_Reference_Genome".
output_dir="${batch_id}" # The output directory for this batch
marker_type="SNP" # The types of markers to use: "BOTH" or "SNP". Default = "SNP".
basecall_purity_cutoff=0.9 # The basecall purity cutoff for genotyping. Default = "0.9".

#############################################


test_directory_existence () {
    dir=$1
    if [[ ! -d $dir ]]
    then
        echo "The specified directory $dir does not exists! Process terminated!"
        exit 1
    fi
}


echo "Testing the existence of the reference genome directory: $genome_dir ..."
test_directory_existence $genome_dir
echo "Testing the existence of the maker directory: $marker_dir ..."
test_directory_existence $marker_dir
echo "Testing the existence of the gamete read mapping genome directory for the specified batch: $gamete_read_mapping_dir/$batch_id ..."
test_directory_existence $gamete_read_mapping_dir/$batch_id


# generate genotypes
echo "generate genotypes ..."
perl $RECOMBINEX_HOME/scripts/batch_tetrad_genotyping_by_reference_genome.pl \
    -s $master_sample_table \
    -q $net_quality_cutoff \
    -p $basecall_purity_cutoff \
    -apply_cnv_filter $apply_cnv_filter \
    -allow_heteroduplex $allow_heteroduplex \
    -c $chr_list \
    -b $batch_id \
    -marker_dir $marker_dir \
    -marker_type $marker_type \
    -gamete_read_mapping_dir $gamete_read_mapping_dir \
    -output_dir $output_dir


echo "generating Rqtl genotype inputs ..."
perl $RECOMBINEX_HOME/scripts/batch_recombinex_genotype2rqtl_genotype_by_reference_genome.pl \
    -s $master_sample_table \
    -q $net_quality_cutoff \
    -b $batch_id \
    -genotype_dir $(pwd) 

echo "summarizing genotype segregation ratio ..."
perl $RECOMBINEX_HOME/scripts/batch_genotype_segregation_summarizing_by_reference_genome.pl \
    -s $master_sample_table \
    -q $net_quality_cutoff \
    -genotype_dir $(pwd) \
    -b $batch_id \
    -output $output_dir/all_samples.segregation_summary.txt

echo "plot genotypes ..."
perl $RECOMBINEX_HOME/scripts/batch_genotype_plotting_by_reference_genome.pl \
    -s $master_sample_table \
    -m raw \
    -genome_dir $genome_dir \
    -q $net_quality_cutoff \
    -genotype_dir $(pwd) \
    -batch_id $batch_id \
    -output_dir $output_dir \
    -color_scheme $color_scheme \
    -plot_centromere $plot_centromere

perl $RECOMBINEX_HOME/scripts/batch_genotype_plotting_by_reference_genome.pl \
    -s $master_sample_table \
    -m inferred \
    -genome_dir $genome_dir \
    -q $net_quality_cutoff \
    -genotype_dir $(pwd) \
    -batch_id $batch_id \
    -output_dir $output_dir \
    -color_scheme $color_scheme \
    -plot_centromere $plot_centromere

# check parental allele frequency
if [[ "$same_cross_combination_for_the_batch" == "yes" ]]
then

    perl $RECOMBINEX_HOME/scripts/batch_parental_allele_frequency_in_tetrads_profiling_by_reference_genome.pl \
	-s $master_sample_table \
	-batch_id $batch_id \
	-qual_diff $net_quality_cutoff \
	-m raw \
	-o $output_dir/$batch_id.parental_allele_frequency.raw.txt
    
    Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_parental_allele_frequency_in_tetrads.R \
	--input $output_dir/$batch_id.parental_allele_frequency.raw.txt \
	--output $output_dir/$batch_id.parental_allele_frequency.raw.plot.pdf \
	--color_scheme $color_scheme
    
    perl $RECOMBINEX_HOME/scripts/batch_parental_allele_frequency_in_tetrads_profiling_by_reference_genome.pl \
	-s $master_sample_table \
	-batch_id $batch_id \
	-qual_diff $net_quality_cutoff \
	-m inferred \
	-o $output_dir/$batch_id.parental_allele_frequency.inferred.txt
    
    Rscript  --vanilla --slave $RECOMBINEX_HOME/scripts/plot_parental_allele_frequency_in_tetrads.R \
	--input $output_dir/$batch_id.parental_allele_frequency.inferred.txt \
	--output $output_dir/$batch_id.parental_allele_frequency.inferred.plot.pdf \
	--color_scheme $color_scheme
    
fi

if [[ -f Rplots.pdf ]]
then
    rm Rplots.pdf
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
    echo "#########################################################################"
    echo ""
    echo "RecombineX message: This bash script has been successfully processed! :)"
    echo ""
    echo "#########################################################################"
    echo ""
    exit 0
fi
############################

