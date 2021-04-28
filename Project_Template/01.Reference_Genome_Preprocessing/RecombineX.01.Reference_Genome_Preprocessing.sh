#!/bin/bash
set -e -o pipefail

##########################################
# load environment variables for RecombineX
source ./../../env.sh

###########################################
# set project-specific variables
reference_genome_assembly="./../00.Reference_Genome/SGDref.genome.fa.gz" # The input reference genome assembly in FASTA format with or without gz compression (e.g. *.fa, *.fasta, *.fa.gz, *.fasta.gz). Default = "./../00.Reference_Genome/SGDref.genome.fa.gz".
use_centromere_annotation="yes" # Whether to use the centromere annotation information. When set to "yes" (default), you need to also provide the path to the centromere GFF file in the "centromere_gff=" option below. Default = "yes".
centromere_gff="./../00.Reference_Genome/SGDref.centromere.gff" # Path to the centromere annotation GFF3 file. Required when "use_centromere_annotation="yes". Otherwise, leave it empty. Default = "./../00.Reference_Genome/SGDref.centromere.gff".
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
check_duplicates_by_windowmasker="no" # Whether to mark duplicated regions using windowmasker. Default = "no". 
ram_for_windowmasker="1536" # Accessible RAM (in MB) for windowmasker. Default = 1536.

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
test_file_existence $reference_genome_assembly
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
perl $RECOMBINEX_HOME/scripts/switch_letter_cases_in_fasta.pl -i $reference_genome_assembly -o ref.genome.raw.fa -c upper

# relabel the genome assembly file
echo ""
echo "relabel sequences in the genome assembly with the genome_tag prefix .."
cat ref.genome.raw.fa |sed "s/>/>ref_/g" > ref.genome.raw.relabel.fa
$samtools_dir/samtools faidx ref.genome.raw.relabel.fa

if [[ $check_duplicates_by_windowmasker == "yes" ]]
then
    $windowmasker_dir/windowmasker -checkdup true -mk_counts -in ref.genome.raw.relabel.fa -out ref.genome.raw.relabel.masking_library.ustat -mem $ram_for_windowmasker
else 
    $windowmasker_dir/windowmasker -mk_counts -in ref.genome.raw.relabel.fa -out ref.genome.raw.relabel.masking_library.ustat -mem $ram_for_windowmasker
fi
    
$windowmasker_dir/windowmasker -ustat ref.genome.raw.relabel.masking_library.ustat -dust true \
                               -in ref.genome.raw.relabel.fa -out ref.genome.softmask.relabel.fa -outfmt fasta
perl $RECOMBINEX_HOME/scripts/softmask2hardmask.pl -i ref.genome.softmask.relabel.fa -o ref.genome.hardmask.relabel.fa
$samtools_dir/samtools faidx ref.genome.hardmask.relabel.fa
### slower solution ###
# perl $RECOMBINEX_HOME/scripts/find_motif_in_genome.pl -i ref.genome.hardmask.relabel.fa -m $RECOMBINEX_HOME/data/hardmask.motif.txt -p ref.genome.hardmask.relabel
# cat ref.genome.hardmask.relabel.masking_details.txt | tail -n+2 | awk {'print $3, $4-1, $5'} OFS='\t' > ref.genome.hardmask.relabel.masking_details.bed
### faster solution ###
$ucsc_dir/faToTwoBit -long -noMask ref.genome.hardmask.relabel.fa ref.genome.hardmask.relabel.2bit
$ucsc_dir/twoBitInfo -nBed -nBed ref.genome.hardmask.relabel.2bit ref.genome.hardmask.relabel.masking_details.bed

# determine GC range for FREEC
perl $RECOMBINEX_HOME/scripts/prepare_genome_for_FREEC.pl -i ref.genome.raw.relabel.fa -p ref -e $excluded_chr_list_for_cnv_profiling
$samtools_dir/samtools faidx ref.FREEC.fa
$bedtools_dir/bedtools makewindows -g ref.FREEC.fa.fai -w $window_size > ref.FREEC.window.$window_size.bed
$bedtools_dir/bedtools nuc -fi ref.FREEC.fa -bed ref.FREEC.window.$window_size.bed > ref.FREEC.GC_content.txt
perl $RECOMBINEX_HOME/scripts/cal_GC_range_for_FREEC.pl -i ref.FREEC.GC_content.txt -w $window_size -lower_quantile $lower_quantile -upper_quantile $upper_quantile -min_mappability $min_mappability -o ref.FREEC.GC_range.txt 

# mappability calculation by gemtools
$gemtools_dir/gemtools index -t $threads -i ref.FREEC.fa -o ref.FREEC.gem
$gemtools_dir/gem-mappability -T $threads -I ref.FREEC.gem -l $raw_read_length -m 0.02 -e 0.02 -o ref.FREEC



# relabel GFF files
echo ""
echo "relabel the sequence field in the GFF files with the genome_tag prefix .."

if [[ $centromere_gff != "" ]]
then
    perl $RECOMBINEX_HOME/scripts/relabel_chr_in_gff.pl -i $centromere_gff -t ref -o ref.centromere.relabel.gff
fi

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm ref.genome.raw.fa
    rm ref.genome.softmask.relabel.fa
    rm ref.genome.raw.relabel.masking_library.ustat
    rm ref.genome.hardmask.relabel.2bit
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

