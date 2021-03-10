#!/bin/bash
set -e -o pipefail

##########################################
# load environment variables for RecombineX
source ./../../env.sh

###########################################
# set project-specific variables
parent_tag="S288C" # The unique tag (preferably containing only letters and numbers) for this specific parental genome. This tag will be used throughout this project. Default = "S288C".
parent_genome_assembly="./../00.Parent_Genomes/$parent_tag.genome.fa" # Input parental genome assembly in FASTA format with or without gz compression (e.g. *.fa, *.fasta, *.fa.gz, *.fasta.gz). Default = "./../00.Parent_Genomes/$parent_tag.genome.fa".
use_centromere_annotation="yes" # Whether to use the centromere annotation information. When set to "yes" (default), you need to also provide the path to the centromere GFF file in the "centromere_gff=" option below. Default = "yes".
centromere_gff="./../00.Parent_Genomes/$parent_tag.centromere.gff" # Path to the centromere annotation GFF3 file. Required when use_centromere_annotation="yes". Default = "./../00.Parent_Genomes/$parent_tag.centromere.gff".
window_size=250 # The window size for the non-overlapping sliding-window-based CNV profiling. Default = 250 (i.e. 250 bp). 
threads=4 # The number of threads to use. Default = "4".
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
###########################################


# process the pipeline
##########################################
# Normally no need to change the following parameters 
# for CNV profiling
lower_quantile=15
upper_quantile=85
min_mappability=0.85 # The minimal mappability for sliding-window-based CNV profiling. Default = "0.85".
excluded_chr_list_for_cnv_profiling="" # The relative path to the list for specifying chromosomes/scaffolds/contigs to be exclued for CNV profiling. We strongly recommend to exclude the organelle (e.g. Mitochondria and Choloraplast) genomes and plasmids if exists. Use "" if there is no chromosome/scaffold/contig for exclusion. Default = "". 
raw_read_length=100 # RecombineX will fix this value for simplicity. There is no need to adjust it for the actual lengths of your Illumina reads.

test_file_existence () {
    filename=$1
    if [[ ! -f $filename ]]
    then
        echo "the file $filename does not exists! process terminated!"
        exit
    fi
}

echo ""
echo "check the existence of genome_assembly"
test_file_existence $parent_genome_assembly
if [[ $use_centromere_annotation == "yes" ]] 
then
    echo ""
    echo "check the existence of centromere_gff"
    test_file_existence $centromere_gff
fi

if [[ $centromere_gff != "" ]] 
then
    echo ""
    echo "check the existence of centromere_gff"
    test_file_existence $centromere_gff
fi

echo "covert all the bases in the input genome assembly into upercases first .."
perl $RECOMBINEX_HOME/scripts/switch_letter_cases_in_fasta.pl -i $parent_genome_assembly -o $parent_tag.genome.raw.fa -c upper

# relabel the genome assembly file
echo ""
echo "relabel sequences in the genome assembly with the genome_tag prefix .."
cat $parent_tag.genome.raw.fa |sed "s/>/>${parent_tag}_/g" > $parent_tag.genome.raw.relabel.fa
$samtools_dir/samtools faidx $parent_tag.genome.raw.relabel.fa

$windowmasker_dir/windowmasker -mk_counts -in $parent_tag.genome.raw.relabel.fa -out $parent_tag.genome.raw.relabel.masking_library.ustat 
$windowmasker_dir/windowmasker -ustat $parent_tag.genome.raw.relabel.masking_library.ustat -dust true \
			       -in $parent_tag.genome.raw.relabel.fa -out $parent_tag.genome.softmask.relabel.fa -outfmt fasta
perl $RECOMBINEX_HOME/scripts/softmask2hardmask.pl -i $parent_tag.genome.softmask.relabel.fa -o $parent_tag.genome.hardmask.relabel.fa
$samtools_dir/samtools faidx $parent_tag.genome.hardmask.relabel.fa
### slower solution ###
# perl $RECOMBINEX_HOME/scripts/find_motif_in_genome.pl -i $parent_tag.genome.hardmask.relabel.fa -m $RECOMBINEX_HOME/data/hardmask.motif.txt -p $parent_tag.genome.hardmask.relabel
# cat $parent_tag.genome.hardmask.relabel.masking_details.txt | tail -n+2 | awk {'print $3, $4-1, $5'} OFS='\t' > $parent_tag.genome.hardmask.relabel.masking_details.bed
### faster solution ###
$ucsc_dir/faToTwoBit -long -noMask $parent_tag.genome.hardmask.relabel.fa $parent_tag.genome.hardmask.relabel.2bit
$ucsc_dir/twoBitInfo -nBed -nBed $parent_tag.genome.hardmask.relabel.2bit $parent_tag.genome.hardmask.relabel.masking_details.bed

# determine GC range for FREEC
perl $RECOMBINEX_HOME/scripts/prepare_genome_for_FREEC.pl -i $parent_tag.genome.raw.relabel.fa -p $parent_tag -e $excluded_chr_list_for_cnv_profiling
$samtools_dir/samtools faidx $parent_tag.FREEC.fa
$bedtools_dir/bedtools makewindows -g $parent_tag.FREEC.fa.fai -w $window_size > $parent_tag.FREEC.window.$window_size.bed
$bedtools_dir/bedtools nuc -fi $parent_tag.FREEC.fa -bed $parent_tag.FREEC.window.$window_size.bed > $parent_tag.FREEC.GC_content.txt
perl $RECOMBINEX_HOME/scripts/cal_GC_range_for_FREEC.pl -i $parent_tag.FREEC.GC_content.txt -w $window_size -lower_quantile $lower_quantile -upper_quantile $upper_quantile -min_mappability $min_mappability -o $parent_tag.FREEC.GC_range.txt 

# mappability calculation by gemtools
$gemtools_dir/gemtools index -t $threads -i $parent_tag.FREEC.fa -o $parent_tag.FREEC.gem
$gemtools_dir/gem-mappability -T $threads -I $parent_tag.FREEC.gem -l $raw_read_length -m 0.02 -e 0.02 -o $parent_tag.FREEC



# relabel GFF files
echo ""
echo "relabel the sequence field in the GFF files with the genome_tag prefix .."

if [[ $centromere_gff != "" ]]
then
    perl $RECOMBINEX_HOME/scripts/relabel_chr_in_gff.pl -i $centromere_gff -t $parent_tag -o $parent_tag.centromere.relabel.gff
fi

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm $parent_tag.genome.raw.fa
    rm $parent_tag.genome.softmask.relabel.fa
    rm $parent_tag.genome.raw.relabel.masking_library.ustat
    rm $parent_tag.genome.hardmask.relabel.2bit
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

