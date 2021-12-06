#!/bin/bash
set -e -o pipefail

##########################################
# load environment variables for RecombineX
source ./../../env.sh

###########################################
# set project-specific variables
batch_id="Batch_S288C-SK1" # The batch id used for gamete read mapping. Default = "Batch_S288C-SK1". 
master_sample_table="Master_Sample_Table.${batch_id}.txt" # The master sample table for this batch. Default = "Master_Sample_Table.${batch_id}.txt".
marker_dir="./../14.Polymorphic_Markers_by_Consensus" # The relative path to the "12.Polymorphic_Markers_by_Cross_Parent_Genome_Alignment" or "14.Polymorphic_Markers_by_Consensus" directory (for parental-genome-based analysis). Default = ./../14.Polymorphic_Markers_by_Consensus".
net_quality_cutoff=50 # The net quality cutoff for genotyping. Default = "50".
apply_cnv_filter="yes" # Whether to set gamete genotype to NA for potential CNV regions in gametes. Set this option to "no" if the gamete sequencing depth is very low (e.g. <= 1). Default = "yes".
allow_heteroduplex="no" # Whether to consider the possibility of heteroduplex formation. Default = "no". 
chr_list="$RECOMBINEX_HOME/data/Saccharomyces_cerevisiae.chr_list.txt" # The chromosome list for the analyzed genome. Default = "$RECOMBINEX_HOME/data/Saccharomyces_cerevisiae.chr_list.txt".
color_scheme="$RECOMBINEX_HOME/data/Saccharomyces_cerevisiae.color_scheme.txt" # The color scheme to use for plotting genotypes. Default = "$RECOMBINEX_HOME/data/Saccharomyces_cerevisiae.color_scheme.txt".
plot_centromere="yes" # Whether to plot centromere in the generated genotyping plots. Please note that enable this option requires that you have the parent1_tag.centromere.relabel.gff and parent2_tag.centromere.relabel.gff files ready in the "./../11.Parent_Genome_Preprocessing" directory (for parental-genome-based analysis). Default = "yes".
same_cross_combination_for_the_batch="yes" # Wether all the samples in the current batch come from the same cross combination (i.e. shared the same parents): "yes" or "no". When setting as "yes", RecombineX will automatically profle and plot the parental allele frequency for every markers. Default = "yes".
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
###########################################


# process the pipeline
############################################
# Normally no need to change the following parameters
genome_dir="./../11.Parent_Genome_Preprocessing" # The relative path to the "11.Parent_Genome_Preprocessing" directory (for parental-genome-based analysis). Default = "./../11.Parent_Genome_Preprocessing".
gamete_read_mapping_dir="./../15.Gamete_Read_Mapping_to_Parent_Genomes" # The relative path to the "15.Gamete_Read_Mapping_to_Parent_Genomes" directory (for parent-genome-based analysis). Default = "./../15.Gamete_Read_Mapping_to_Parent_Genomes".
basecall_purity_cutoff=0.9 # The basecall purity cutoff for genotyping. Default = "0.9".
############################################

test_directory_existence () {
    dir=$1
    if [[ ! -d $dir ]]
    then
        echo "The specified directory $dir does not exists! Process terminated!"
        exit 1
    fi
}


echo "Testing the existence of the parental genome directory: $genome_dir ..."
test_directory_existence $genome_dir
echo "Testing the existence of the maker directory: $marker_dir ..."
test_directory_existence $marker_dir
echo "Testing the existence of the gamete read mapping genome directory for the specified batch: $gamete_read_mapping_dir/$batch_id ..."
test_directory_existence $gamete_read_mapping_dir/$batch_id


output_dir="${batch_id}" # The output directory for this batch
marker_type="SNP" # The types of markers to use: "SNP" or "BOTH". Default = "SNP".

# generate genotypes
echo "generate genotypes ..."
perl $RECOMBINEX_HOME/scripts/batch_tetrad_genotyping_by_parent_genomes.pl \
    -s $master_sample_table \
    -q $net_quality_cutoff \
    -p $basecall_purity_cutoff \
    -c $chr_list \
    -b $batch_id \
    -apply_cnv_filter $apply_cnv_filter \
    -allow_heteroduplex $allow_heteroduplex \
    -marker_dir $marker_dir \
    -marker_type $marker_type \
    -gamete_read_mapping_dir $gamete_read_mapping_dir \
    -output_dir $output_dir

echo "generating Rqtl genotype inputs ..."
perl $RECOMBINEX_HOME/scripts/batch_recombinex_genotype2rqtl_genotype_by_parent_genomes.pl \
    -s $master_sample_table \
    -q $net_quality_cutoff \
    -b $batch_id \
    -genotype_dir $(pwd)

echo "summarizing genotype segregation ratio ..."
perl $RECOMBINEX_HOME/scripts/batch_genotype_segregation_summarizing_by_parent_genomes.pl \
    -s $master_sample_table \
    -q $net_quality_cutoff \
    -genotype_dir $(pwd) \
    -b $batch_id \
    -output $output_dir/all_samples.segregation_summary.txt

echo "plot genotypes ..."
perl $RECOMBINEX_HOME/scripts/batch_genotype_plotting_by_parent_genomes.pl \
    -s $master_sample_table \
    -m raw \
    -genome_dir $genome_dir \
    -q $net_quality_cutoff \
    -genotype_dir $(pwd) \
    -batch_id $batch_id \
    -output_dir $output_dir \
    -color_scheme $color_scheme \
    -plot_centromere $plot_centromere

perl $RECOMBINEX_HOME/scripts/batch_genotype_plotting_by_parent_genomes.pl \
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
# check parental allele frequency
if [[ "$same_cross_combination_for_the_batch" == "yes" ]]
then

    perl $RECOMBINEX_HOME/scripts/batch_parental_allele_frequency_in_tetrads_profiling_by_parent_genomes.pl \
	-s $master_sample_table \
	-batch_id $batch_id \
	-qual_diff $net_quality_cutoff \
	-m raw \
	-o $output_dir/$batch_id.parental_allele_frequency.raw.txt
    
    Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_parental_allele_frequency_in_tetrads.R \
	--input $output_dir/$batch_id.parental_allele_frequency.raw.txt \
	--output $output_dir/$batch_id.parental_allele_frequency.raw.plot.pdf \
	--color_scheme $color_scheme
    
    perl $RECOMBINEX_HOME/scripts/batch_parental_allele_frequency_in_tetrads_profiling_by_parent_genomes.pl \
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

