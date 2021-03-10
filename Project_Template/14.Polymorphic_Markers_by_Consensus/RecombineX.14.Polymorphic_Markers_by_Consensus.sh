#!/bin/bash
set -e -o pipefail

##########################################
# load environment variables for RecombineX
source ./../../env.sh

###########################################
# set project-specific variables
parent1_tag="S288C" # The tag for parental genome 1 used in 11.Parent_Genome_Preprocessing. Default = "S288C".
parent2_tag="SK1" # The tag for parental genome 2 used in 11.Parent_Genome_Preprocessing. Default = "SK1".
use_centromere_annotation="yes" # Whether to use centromere annotation. Please note that enable this option requires that you have the parent1_tag.centromere.relabel.gff and parent2_tag.centromere.relabel.gff files ready in the "./../11.Parent_Genome_Preprocessing" directory. Default = "yes".
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
###########################################



# process the pipeline
##########################################
# Normally no need to change the following parameters
raw_parent1_assembly="./../11.Parent_Genome_Preprocessing/$parent1_tag.genome.raw.relabel.fa" # relabeled raw genome 1 assembly file generated in the 01.Parent_Genome_Preprocessing directory.
raw_parent2_assembly="./../11.Parent_Genome_Preprocessing/$parent2_tag.genome.raw.relabel.fa" # relabeled raw genome 2 assembly file generated in the 01.Parent_Genome_Preprocessing directory.
##########################################


test_file_existence () {
    filename=$1
    if [[ ! -f $filename ]]
    then
        echo "The file $filename does not exists! process terminated!"
        exit
    else
	echo "Test passed!"
    fi
}

if [[ $use_centromere_annotation == "yes" ]]
then
    parent1_centromere_gff="./../11.Parent_Genome_Preprocessing/$parent1_tag.centromere.relabel.gff"
    parent2_centromere_gff="./../11.Parent_Genome_Preprocessing/$parent2_tag.centromere.relabel.gff"
fi

echo ""
echo "Testing the existence of raw_parent1_assembly: $raw_parent1_assembly"
test_file_existence $raw_parent1_assembly
echo ""
echo "Testing the existence of raw_parent2_assembly: $raw_parent2_assembly"
test_file_existence $raw_parent2_assembly

if [[ $use_centromere_annotation == "yes" ]] 
then
    echo ""
    echo "Testing the existence of parent1_centromere_gff: $parent1_centromere_gff"
    test_file_existence $parent1_centromere_gff
    echo ""
    echo "Testing the existence of parent2_centromere_gff: $parent2_centromere_gff"
    test_file_existence $parent2_centromere_gff
fi

echo ""
parent1_based_prefix="${parent1_tag}-${parent2_tag}.$parent1_tag"
parent2_based_prefix="${parent1_tag}-${parent2_tag}.$parent2_tag"

parent1_based_snp_vcf_by_cross_parent_genome_alignment="./../12.Polymorphic_Markers_by_Cross_Parent_Genome_Alignment/$parent1_based_prefix.final.SNP.markers.vcf.gz"
# parent1_based_indel_vcf_by_cross_parent_genome_alignment="./../12.Polymorphic_Markers_by_Cross_Parent_Genome_Alignment/$parent1_based_prefix.final.INDEL.markers.vcf.gz"
parent2_based_snp_vcf_by_cross_parent_genome_alignment="./../12.Polymorphic_Markers_by_Cross_Parent_Genome_Alignment/$parent2_based_prefix.final.SNP.markers.vcf.gz"
# parent2_based_indel_vcf_by_cross_parent_genome_alignment="./../12.Polymorphic_Markers_by_Cross_Parent_Genome_Alignment/$parent2_based_prefix.final.INDEL.markers.vcf.gz"

parent1_based_snp_vcf_by_cross_parent_read_mapping="./../13.Polymorphic_Markers_by_Cross_Parent_Read_Mapping/$parent1_based_prefix.read_mapping.SNP.filter.vcf.gz"
# parent1_based_indel_vcf_by_cross_parent_read_mapping="./../13.Polymorphic_Markers_by_Cross_Parent_Read_Mapping/$parent1_based_prefix.read_mapping.INDEL.filter.vcf.gz"
parent2_based_snp_vcf_by_cross_parent_read_mapping="./../13.Polymorphic_Markers_by_Cross_Parent_Read_Mapping/$parent2_based_prefix.read_mapping.SNP.filter.vcf.gz"
# parent2_based_indel_vcf_by_cross_parent_read_mapping="./../13.Polymorphic_Markers_by_Cross_Parent_Read_Mapping/$parent2_based_prefix.read_mapping.INDEL.filter.vcf.gz"

$bedtools_dir/bedtools intersect -wa -f 1 -r -a $parent1_based_snp_vcf_by_cross_parent_genome_alignment -b $parent1_based_snp_vcf_by_cross_parent_read_mapping >$parent1_based_prefix.intersect.SNP.markers.vcf
# $bedtools_dir/bedtools intersect -wa -f 1 -r -a $parent1_based_indel_vcf_by_cross_parent_genome_alignment -b $parent1_based_indel_vcf_by_cross_parent_read_mapping >$parent1_based_prefix.intersect.INDEL.markers.vcf

$bedtools_dir/bedtools intersect -wa -f 1 -r -a $parent2_based_snp_vcf_by_cross_parent_genome_alignment -b $parent2_based_snp_vcf_by_cross_parent_read_mapping >$parent2_based_prefix.intersect.SNP.markers.vcf
# $bedtools_dir/bedtools intersect -wa -f 1 -r -a $parent2_based_indel_vcf_by_cross_parent_genome_alignment -b $parent2_based_indel_vcf_by_cross_parent_read_mapping >$parent2_based_prefix.intersect.INDEL.markers.vcf

gzip *.vcf

cp $raw_parent1_assembly $parent1_tag.genome.fa
cp $raw_parent2_assembly $parent2_tag.genome.fa

$samtools_dir/samtools faidx $parent1_tag.genome.fa
$samtools_dir/samtools faidx $parent2_tag.genome.fa

parent1_fai="$parent1_tag.genome.fa.fai"
parent2_fai="$parent2_tag.genome.fa.fai"


# apply depth filter
perl $RECOMBINEX_HOME/scripts/filter_parent_based_markers_by_depth.pl \
     -i $parent1_based_prefix.intersect.SNP.markers.vcf.gz \
     -depth_summary ./../13.Polymorphic_Markers_by_Cross_Parent_Read_Mapping/${parent1_tag}-${parent2_tag}_${parent1_tag}_based/$parent1_based_prefix.coverage_summary.txt \
     -depth_detail ./../13.Polymorphic_Markers_by_Cross_Parent_Read_Mapping/${parent1_tag}-${parent2_tag}_${parent1_tag}_based/$parent1_based_prefix.depth.txt.gz \
     -o $parent1_based_prefix.depth_filter.SNP.markers.vcf.gz

# perl $RECOMBINEX_HOME/scripts/filter_parent_based_markers_by_depth.pl \
#      -i $parent1_based_prefix.intersect.INDEL.markers.vcf.gz \
#      -depth_summary ./../13.Polymorphic_Markers_by_Cross_Parent_Read_Mapping/${parent1_tag}-${parent2_tag}_${parent1_tag}_based/$parent1_based_prefix.coverage_summary.txt \
#      -depth_detail ./../13.Polymorphic_Markers_by_Cross_Parent_Read_Mapping/${parent1_tag}-${parent2_tag}_${parent1_tag}_based/$parent1_based_prefix.depth.txt.gz \
#      -o $parent1_based_prefix.depth_filter.INDEL.markers.vcf.gz

perl $RECOMBINEX_HOME/scripts/filter_parent_based_markers_by_depth.pl \
     -i $parent2_based_prefix.intersect.SNP.markers.vcf.gz \
     -depth_summary ./../13.Polymorphic_Markers_by_Cross_Parent_Read_Mapping/${parent1_tag}-${parent2_tag}_${parent2_tag}_based/$parent2_based_prefix.coverage_summary.txt \
     -depth_detail ./../13.Polymorphic_Markers_by_Cross_Parent_Read_Mapping/${parent1_tag}-${parent2_tag}_${parent2_tag}_based/$parent2_based_prefix.depth.txt.gz \
     -o $parent2_based_prefix.depth_filter.SNP.markers.vcf.gz
# perl $RECOMBINEX_HOME/scripts/filter_parent_based_markers_by_depth.pl \
#      -i $parent2_based_prefix.intersect.INDEL.markers.vcf.gz \
#      -depth_summary ./../13.Polymorphic_Markers_by_Cross_Parent_Read_Mapping/${parent1_tag}-${parent2_tag}_${parent2_tag}_based/$parent2_based_prefix.coverage_summary.txt \
#      -depth_detail ./../13.Polymorphic_Markers_by_Cross_Parent_Read_Mapping/${parent1_tag}-${parent2_tag}_${parent2_tag}_based/$parent2_based_prefix.depth.txt.gz \
#      -o $parent2_based_prefix.depth_filter.INDEL.markers.vcf.gz

# apply the reciprocal filter
perl $RECOMBINEX_HOME/scripts/filter_nonreciprocal_markers.pl \
     -g1_based_vcf $parent1_based_prefix.depth_filter.SNP.markers.vcf.gz \
     -g2_based_vcf $parent2_based_prefix.depth_filter.SNP.markers.vcf.gz \
     -g1 $parent1_tag.genome.fa \
     -g2 $parent2_tag.genome.fa \
     -g1_based_prefix $parent1_based_prefix.final.SNP \
     -g2_based_prefix $parent2_based_prefix.final.SNP
# perl $RECOMBINEX_HOME/scripts/filter_nonreciprocal_markers.pl \
#      -g1_based_vcf $parent1_based_prefix.intersect.INDEL.markers.vcf.gz \
#      -g2_based_vcf $parent2_based_prefix.intersect.INDEL.markers.vcf.gz \
#      -g1 $parent1_tag.genome.fa \
#      -g2 $parent2_tag.genome.fa \
#      -g1_based_prefix $parent1_based_prefix.final.INDEL \
#      -g2_based_prefix $parent2_based_prefix.final.INDEL




gzip *.markers.vcf
gzip *.markers.txt


# calculate intermarker distance

perl $RECOMBINEX_HOME/scripts/intermarker_distance.pl -i $parent1_based_prefix.final.SNP.markers.txt.gz \
     -p $parent1_based_prefix.final.SNP.markers
# perl $RECOMBINEX_HOME/scripts/intermarker_distance.pl -i $parent1_based_prefix.final.INDEL.markers.txt.gz \
#      -p $parent1_based_prefix.final.INDEL.markers
# zcat $parent1_based_prefix.final.INDEL.markers.txt.gz | tail -n +2 | gzip > $parent1_based_prefix.final.INDEL.markers.content.txt.gz
# cat $parent1_based_prefix.final.SNP.markers.txt.gz $parent1_based_prefix.final.INDEL.markers.content.txt.gz > $parent1_based_prefix.final.BOTH.markers.txt.gz
# rm $parent1_based_prefix.final.INDEL.markers.content.txt.gz
# perl $RECOMBINEX_HOME/scripts/intermarker_distance.pl -i $parent1_based_prefix.final.BOTH.markers.txt.gz \
#      -p $parent1_based_prefix.final.BOTH.markers

perl $RECOMBINEX_HOME/scripts/intermarker_distance.pl -i $parent2_based_prefix.final.SNP.markers.txt.gz \
     -p $parent2_based_prefix.final.SNP.markers
# perl $RECOMBINEX_HOME/scripts/intermarker_distance.pl -i $parent2_based_prefix.final.INDEL.markers.txt.gz \
#      -p $parent2_based_prefix.final.INDEL.markers

# zcat $parent2_based_prefix.final.INDEL.markers.txt.gz | tail -n +2 | gzip > $parent2_based_prefix.final.INDEL.markers.content.txt.gz
# cat $parent2_based_prefix.final.SNP.markers.txt.gz $parent2_based_prefix.final.INDEL.markers.content.txt.gz > $parent2_based_prefix.final.BOTH.markers.txt.gz
# rm $parent2_based_prefix.final.INDEL.markers.content.txt.gz
# perl $RECOMBINEX_HOME/scripts/intermarker_distance.pl -i $parent2_based_prefix.final.BOTH.markers.txt.gz \
#      -p $parent2_based_prefix.final.BOTH.markers


# plot markers
if [[ $use_centromere_annotation == "yes" ]]
then
    Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent1_fai \
	    --marker_table $parent1_based_prefix.final.SNP.markers.txt.gz \
	    --centromere_gff $parent1_centromere_gff --output_prefix $parent1_based_prefix.final.SNP.markers
    rm Rplots.pdf
    # Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent1_fai \
    # 	    --marker_table $parent1_based_prefix.final.INDEL.markers.txt.gz \
    # 	    --centromere_gff $parent1_centromere_gff --output_prefix $parent1_based_prefix.final.INDEL.markers
    # rm Rplots.pdf
    # Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent1_fai \
    # 	    --marker_table $parent1_based_prefix.final.BOTH.markers.txt.gz \
    # 	    --centromere_gff $parent1_centromere_gff --output_prefix $parent1_based_prefix.final.BOTH.markers
    # rm Rplots.pdf

    Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent2_fai \
	    --marker_table $parent2_based_prefix.final.SNP.markers.txt.gz \
	    --centromere_gff $parent2_centromere_gff --output_prefix $parent2_based_prefix.final.SNP.markers
    rm Rplots.pdf
    # Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent2_fai \
    # 	    --marker_table $parent2_based_prefix.final.INDEL.markers.txt.gz \
    # 	    --centromere_gff $parent2_centromere_gff --output_prefix $parent2_based_prefix.final.INDEL.markers
    # rm Rplots.pdf
    # Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent2_fai \
    # 	    --marker_table $parent2_based_prefix.final.BOTH.markers.txt.gz \
    # 	    --centromere_gff $parent2_centromere_gff --output_prefix $parent2_based_prefix.final.BOTH.markers
    # rm Rplots.pdf

else
    Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent1_fai \
	    --marker_table $parent1_based_prefix.final.SNP.markers.txt.gz \
	    --output_prefix $parent1_based_prefix.final.SNP.markers
    rm Rplots.pdf
    # Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent1_fai \
    # 	    --marker_table $parent1_based_prefix.final.INDEL.markers.txt.gz \
    # 	    --output_prefix $parent1_based_prefix.final.INDEL.markers
    # rm Rplots.pdf
    # Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent1_fai \
    # 	    --marker_table $parent1_based_prefix.final.BOTH.markers.txt.gz \
    # 	    --output_prefix $parent1_based_prefix.final.BOTH.markers
    # rm Rplots.pdf

    Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent2_fai \
	    --marker_table $parent2_based_prefix.final.SNP.markers.txt.gz \
	    --output_prefix $parent2_based_prefix.final.SNP.markers
    rm Rplots.pdf
    # Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent2_fai \
    # 	    --marker_table $parent2_based_prefix.final.INDEL.markers.txt.gz \
    # 	    --output_prefix $parent2_based_prefix.final.INDEL.markers
    # rm Rplots.pdf
    # Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent2_fai \
    # 	    --marker_table $parent2_based_prefix.final.BOTH.markers.txt.gz \
    # 	    --output_prefix $parent2_based_prefix.final.BOTH.markers
    # rm Rplots.pdf
    
fi




# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm $parent1_tag.genome.fa
    rm $parent2_tag.genome.fa
    rm $parent1_fai
    rm $parent2_fai
    # rm $parent1_based_prefix.final.BOTH.markers.txt.gz
    # rm $parent2_based_prefix.final.BOTH.markers.txt.gz
    # rm $parent1_based_prefix.intersect.*.markers.vcf.gz
    # rm $parent2_based_prefix.intersect.*.markers.vcf.gz
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


