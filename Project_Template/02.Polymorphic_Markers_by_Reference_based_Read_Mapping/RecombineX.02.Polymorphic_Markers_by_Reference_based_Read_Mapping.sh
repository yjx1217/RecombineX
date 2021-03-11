#!/bin/bash
set -e -o pipefail

##########################################
# load environment variables for RecombineX
source ./../../env.sh

###########################################
# set project-specific variables
parent1_tag="S288C" # The unique tag of parental genome 1. This tag will be used throughout the project. Default = "S288C".
parent2_tag="SK1" # The unique tag of parental genome 2. This tag will be used throughout the project. Default = "SK1".
parent1_reads_R1="./../00.Parent_Reads/$parent1_tag.R1.fq.gz" # The path to the paired-end R1 reads of the parental genome 1. Default = "./../00.Parent_Reads/$parent1_tag.R1.fq.gz".
parent1_reads_R2="./../00.Parent_Reads/$parent1_tag.R2.fq.gz" # The path to the paired-end R2 reads of the parental genome 1. Default = "./../00.Parent_Reads/$parent1_tag.R2.fq.gz".
parent2_reads_R1="./../00.Parent_Reads/$parent2_tag.R1.fq.gz" # The path to the paired-end R1 reads of the parental genome 2. Default = "./../00.Parent_Reads/$parent2_tag.R1.fq.gz".
parent2_reads_R2="./../00.Parent_Reads/$parent2_tag.R2.fq.gz" # The path to the paired-end R2 reads of the parental genome 2. Default = "./../00.Parent_Reads/$parent2_tag.R2.fq.gz".
use_centromere_annotation="yes" # Whether to use the centromere annotation information. Please note that enabling this option requires that you have the parent1_tag.centromere.relabel.gff and parent2_tag.centromere.relabel.gff files ready in the "./../01.Reference_Genome_Preprocessing" directory. Default = "yes".
ploidy=1 # The ploidy of the parental genome. (e.g. "1" for haploid and "2" for diploid). For diploid parents, only homozygous SNPs will be used as markers. If the parental genome is purely homozygous, it is recommended to set "ploidy=1" to maximize the power of CNV profiling. Default = "1".
window_size=250 # The window size for the non-overlapping sliding-window-based CNV profiling. Default = 250 (i.e. 250 bp). 
threads=4 # The number of threads to use. Default = "4".
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
##########################################


# process the pipeline
##########################################
# Normally no need to change the following parameters 
reference_genome_preprocessing_dir="./../01.Reference_Genome_Preprocessing" # The path to the 01.Reference_Genome_Preprocessing directory.
reference_raw_assembly="$reference_genome_preprocessing_dir/ref.genome.raw.relabel.fa" # The relabeled reference genome raw assembly generated in the 01.Reference_Genome_Preprocessing directory.
reference_hardmask_bed="$reference_genome_preprocessing_dir/ref.genome.hardmask.relabel.masking_details.bed" # The masking details bed file for the relabeled and hardmasked reference genome assembly file generated in the 01.Reference_Genome_Preprocessing directory.
excluded_chr_list_for_cnv_profiling="" # The relative path to the list for specifying chromosomes/scaffolds/contigs to be exclued for CNV profiling. We strongly recommend to exclude the organelle (e.g. Mitochondria and Choloraplast) genomes and plasmids if exists. Use "" if there is no chromosome/scaffold/contig for exclusion. Default = "". 
mapping_quality_cutoff=30 # The minimal mapping quality to be considered. Default = "30".
variant_calling_quality_cutoff=30 # The minimal variant calling quality to be considered. Default = "30".

min_mappability=0.85 # The minimal mappability for sliding-window-based CNV profiling. Default = "0.85".
cluster_window_size=10 # Adjacent variants within the specified window (unit: bp) will be both filtered out if any of them is INDEL. Default = "10".
#######################################


test_file_existence () {
    filename=$1
    if [[ ! -f $filename ]]
    then
        echo "the file $filename does not exists! process terminated!"
        exit
    else
	echo "test pass."
    fi
}

echo ""
echo "check the existence of reference_raw_assembly: $reference_raw_assembly"
test_file_existence $reference_raw_assembly
echo ""
echo "check the existence of reference_hardmask_bed: $reference_hardmask_bed"
test_file_existence $reference_hardmask_bed

if [[ $use_centromere_annotation == "yes" ]]
then
    reference_centromere_gff="$reference_genome_preprocessing_dir/ref.centromere.relabel.gff"
    echo ""
    echo "check the existence of reference_centromere_gff: $reference_centromere_gff"
    test_file_existence $reference_centromere_gff
fi

echo ""
echo "check the existence of genome1_reads_R1: $parent1_reads_R1"
test_file_existence $parent1_reads_R1
echo ""
echo "check the existence of parent1_reads_R2: $parent1_reads_R2"
test_file_existence $parent1_reads_R2
echo ""
echo "check the existence of parent2_reads_R1: $parent2_reads_R1"
test_file_existence $parent2_reads_R1
echo ""
echo "check the existence of parent2_reads_R2: $parent2_reads_R2"
test_file_existence $parent2_reads_R2
echo ""

reference_based_output_dir="${parent1_tag}-${parent2_tag}_Reference_based"
mkdir $reference_based_output_dir

cd $reference_based_output_dir

ln -s ./../$reference_raw_assembly ref.genome.raw.fa
ln -s ./../$parent1_reads_R1 $parent1_tag.R1.raw.fq.gz
ln -s ./../$parent1_reads_R2 $parent1_tag.R2.raw.fq.gz
ln -s ./../$parent2_reads_R1 $parent2_tag.R1.raw.fq.gz
ln -s ./../$parent2_reads_R2 $parent2_tag.R2.raw.fq.gz

adapter="$trimmomatic_dir/adapters/TruSeq3-PE-2.fa"
cp $adapter adapter.fa

mkdir tmp

# index reference sequence
$samtools_dir/samtools faidx ref.genome.raw.fa
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar CreateSequenceDictionary  \
    -REFERENCE ref.genome.raw.fa \
    -OUTPUT ref.genome.raw.dict
$bwa_dir/bwa index ref.genome.raw.fa

# trim the reads
$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $trimmomatic_dir/trimmomatic.jar PE -threads $threads -phred33 $parent1_tag.R1.raw.fq.gz $parent1_tag.R2.raw.fq.gz $parent1_tag.R1.trimmed.PE.fq.gz $parent1_tag.R1.trimmed.SE.fq.gz $parent1_tag.R2.trimmed.PE.fq.gz $parent1_tag.R2.trimmed.SE.fq.gz ILLUMINACLIP:adapter.fa:2:30:10  SLIDINGWINDOW:5:20 MINLEN:36

if [[ $debug == "no" ]]
then
    rm $parent1_tag.R1.trimmed.SE.fq.gz
    rm $parent1_tag.R2.trimmed.SE.fq.gz
fi

# map reads to the reference genome
$bwa_dir/bwa mem -t $threads -M ref.genome.raw.fa  $parent1_tag.R1.trimmed.PE.fq.gz $parent1_tag.R2.trimmed.PE.fq.gz | $samtools_dir/samtools view -bS -q $mapping_quality_cutoff - >${parent1_tag}-ref.ref.bam

if [[ $debug == "no" ]]
then
    rm $parent1_tag.R1.trimmed.PE.fq.gz
    rm $parent1_tag.R2.trimmed.PE.fq.gz
fi

# sort bam file by picard-tools SortSam
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar SortSam \
    -INPUT ${parent1_tag}-ref.ref.bam \
    -OUTPUT ${parent1_tag}-ref.ref.sort.bam \
    -SORT_ORDER coordinate

if [[ $debug == "no" ]]
then
    rm ${parent1_tag}-ref.ref.bam
fi

# fixmate
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar FixMateInformation \
    -INPUT ${parent1_tag}-ref.ref.sort.bam \
    -OUTPUT ${parent1_tag}-ref.ref.fixmate.bam

if [[ $debug == "no" ]]
then
    rm ${parent1_tag}-ref.ref.sort.bam
fi

# add or replace read groups and sort
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -jar $picard_dir/picard.jar AddOrReplaceReadGroups \
    -INPUT ${parent1_tag}-ref.ref.fixmate.bam \
    -OUTPUT ${parent1_tag}-ref.ref.rdgrp.bam \
    -SORT_ORDER coordinate \
    -RGID "$parent1_tag" \
    -RGLB "$parent1_tag" \
    -RGPL "Illumina" \
    -RGPU "$parent1_tag" \
    -RGSM "$parent1_tag" \
    -RGCN "RGCN" 

if [[ $debug == "no" ]]
then
    rm ${parent1_tag}-ref.ref.fixmate.bam
fi

# remove duplicates
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar MarkDuplicates \
    -INPUT ${parent1_tag}-ref.ref.rdgrp.bam \
    -REMOVE_DUPLICATES true  \
    -METRICS_FILE ${parent1_tag}-ref.ref.dedup.matrics  \
    -OUTPUT ${parent1_tag}-ref.ref.dedup.bam 

# index the dedup.bam file
$samtools_dir/samtools index ${parent1_tag}-ref.ref.dedup.bam

if [[ $debug == "no" ]]
then
    rm ${parent1_tag}-ref.ref.rdgrp.bam
fi

# GATK local realign
# find realigner targets
$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $gatk3_dir/gatk3.jar \
    -nt $threads \
    -R ref.genome.raw.fa \
    -T RealignerTargetCreator \
    -I ${parent1_tag}-ref.ref.dedup.bam  \
    -o ${parent1_tag}-ref.ref.realn.intervals 
# run realigner
$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $gatk3_dir/gatk3.jar  \
    -R ref.genome.raw.fa -T IndelRealigner \
    -I ${parent1_tag}-ref.ref.dedup.bam -targetIntervals ${parent1_tag}-ref.ref.realn.intervals  \
    -o ${parent1_tag}-ref.ref.realn.bam

if [[ $debug == "no" ]]
then
    rm ${parent1_tag}-ref.ref.dedup.bam
    rm ${parent1_tag}-ref.ref.dedup.bam.bai
    rm ${parent1_tag}-ref.ref.dedup.matrics
    rm ${parent1_tag}-ref.ref.realn.intervals
fi

# generate samtools mpileup
$samtools_dir/samtools mpileup -C 0 -q $mapping_quality_cutoff -f ref.genome.raw.fa ${parent1_tag}-ref.ref.realn.bam |gzip -c >${parent1_tag}-ref.ref.mpileup.gz

# calculate per-base depth
$samtools_dir/samtools depth -aa ${parent1_tag}-ref.ref.realn.bam |gzip -c >${parent1_tag}-ref.ref.depth.txt.gz

# compute basic alignment statistics by samtools
$samtools_dir/samtools flagstat ${parent1_tag}-ref.ref.realn.bam >${parent1_tag}-ref.ref.samstat

# compute insert size statistics
$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar CollectInsertSizeMetrics \
    I=${parent1_tag}-ref.ref.realn.bam \
    O=${parent1_tag}-ref.ref.insert_size_metrics.txt \
    H=${parent1_tag}-ref.ref.insert_size_histogram.pdf \
    M=0.5

# calculate read mapping coverage statistics
perl $RECOMBINEX_HOME/scripts/summarize_mapping_coverage.pl \
    -r ref.genome.raw.fa \
    -s ${parent1_tag}-ref.ref.samstat \
    -d ${parent1_tag}-ref.ref.depth.txt.gz \
    -c 5 \
    -t $parent1_tag -o ${parent1_tag}-ref.ref.coverage_summary.txt


# scan for CNV by FREEC
#parent1_raw_read_length=$(gunzip -c $parent1_tag.R1.raw.fq.gz |awk 'NR%4==2{print length($0)}' | head -1 || true)
# here we fix read length to 100bp for simplicity
parent1_raw_read_length=100
step_size=$window_size
min_expected_gc=$(cat ./../$reference_genome_preprocessing_dir/ref.FREEC.GC_range.txt|egrep -v "#"|cut -f 1)
max_expected_gc=$(cat ./../$reference_genome_preprocessing_dir/ref.FREEC.GC_range.txt|egrep -v "#"|cut -f 2)
echo "min_expected_gc=$min_expected_gc, max_expected_gc=$max_expected_gc";

if [ -z "$excluded_chr_list_for_cnv_profiling" ]
then
    perl $RECOMBINEX_HOME/scripts/run_FREEC_wrapper_lite_for_RecombineX.pl \
	-r ref \
	-bam ${parent1_tag}-ref.ref.realn.bam \
	-prefix ${parent1_tag}-ref.ref \
	-ploidy $ploidy \
	-bedtools $bedtools_dir/bedtools \
	-samtools $samtools_dir/samtools \
	-freec $freec_dir/freec \
	-window $window_size \
	-step $step_size \
	-read_length_for_mappability $parent1_raw_read_length \
	-min_mappability $min_mappability \
	-min_expected_gc $min_expected_gc \
	-max_expected_gc $max_expected_gc \
	-mates_orientation 0 \
	-refseq_genome_preprocessing_dir "./../$reference_genome_preprocessing_dir" \
	-threads $threads \
then
    perl $RECOMBINEX_HOME/scripts/run_FREEC_wrapper_lite_for_RecombineX.pl \
	-r ref \
	-bam ${parent1_tag}-ref.ref.realn.bam \
	-prefix ${parent1_tag}-ref.ref \
	-ploidy $ploidy \
	-bedtools $bedtools_dir/bedtools \
	-samtools $samtools_dir/samtools \
	-freec $freec_dir/freec \
	-window $window_size \
	-step $step_size \
	-read_length_for_mappability $parent1_raw_read_length \
	-min_mappability $min_mappability \
	-min_expected_gc $min_expected_gc \
	-max_expected_gc $max_expected_gc \
	-mates_orientation 0 \
	-refseq_genome_preprocessing_dir "./../$reference_genome_preprocessing_dir" \
	-threads $threads \
	-excluded_chr_list ./../$excluded_chr_list_for_cnv_profiling
fi
mv FREEC.config.txt ${parent1_tag}-ref.ref.FREEC.config.txt

Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/CNV_segmentation_by_DNAcopy.R \
    --input ${parent1_tag}-ref.ref.FREEC.bam_ratio.txt \
    --prefix ${parent1_tag}-ref.ref \
    --window $window_size \
    --ploidy $ploidy \
    --genome_fai ref.genome.raw.fa.fai

perl $RECOMBINEX_HOME/scripts/adjust_FREEC_copynumber_by_DNAcopy_copynumber.pl \
    -i ${parent1_tag}-ref.ref.FREEC.bam_ratio.sorted.txt \
    -a ${parent1_tag}-ref.ref.FREEC.bam_ratio.sorted.resegmented.lite.txt \
    -o ${parent1_tag}-ref.ref.FREEC.bam_ratio.sorted.adjusted.txt

if [[ -s ${parent1_tag}-ref.ref.FREEC.bam_ratio.sorted.adjusted.txt ]]
then
    Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_CNV_for_FREEC.R \
	--ploidy $ploidy \
	--genome_fai ref.genome.raw.fa.fai \
	--input ${parent1_tag}-ref.ref.FREEC.bam_ratio.sorted.adjusted.txt \
	--output ${parent1_tag}-ref.ref.CNV_plot.pdf
    rm Rplots.pdf
    if [[ -s ${parent1_tag}-ref.ref.FREEC.bam_ratio.sorted.resegmented.CNVs.txt ]]
    then
        Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/assess_CNV_significance_for_FREEC.R \
	    --cnv ${parent1_tag}-ref.ref.FREEC.bam_ratio.sorted.resegmented.CNVs.txt \
	    --ratio ${parent1_tag}-ref.ref.FREEC.bam_ratio.sorted.adjusted.txt \
	    --genome_fai ref.genome.raw.fa.fai \
	    --output ${parent1_tag}-ref.ref.CNV_significance_test.txt
	cat ${parent1_tag}-ref.ref.CNV_significance_test.txt | tail -n +2 | awk '{ if ($6 < 0.05) print $0 }' | awk '{print $1, $2-1, $3}' OFS='\t' > ${parent1_tag}-ref.ref.significant_CNV.bed
    else
        echo -e "chr\tstart\tend\tcopy_number\tstatus\tMWU_test_p_value\tKS_test_p_value" > ${parent1_tag}-ref.ref.CNV_significance_test.txt
        echo -e "\t\t" > ${parent1_tag}-ref.ref.significant_CNV.bed
        echo 
    fi
else
    echo "Exception encountered for FREEC! Exit! ..." > ${parent1_tag}-ref.ref.realn.bam.no_FREEC.txt
    exit
fi


mkdir ${parent1_tag}-ref.ref.FREEC_intermediate_files
mv for_CNV.bam ${parent1_tag}-ref.ref.FREEC_intermediate_files
mv GC_profile.${window_size}bp.cnp ${parent1_tag}-ref.ref.FREEC_intermediate_files
mv FREEC.bam* ${parent1_tag}-ref.ref.FREEC_intermediate_files

if [[ $debug == "no" ]]
then
    rm -r ${parent1_tag}-ref.ref.FREEC_intermediate_files
    rm $parent1_tag.R1.raw.fq.gz
    rm $parent1_tag.R2.raw.fq.gz
fi

# SNP and INDEL calling

python3 $freebayes_dir/../scripts/fasta_generate_regions.py ref.genome.raw.fa.fai 100000 > ref.genome.region.txt
cat ref.genome.region.txt | $parallel_dir/parallel -k -j $threads $freebayes_dir/freebayes -f ref.genome.raw.fa  -p $ploidy ${parent1_tag}-ref.ref.realn.bam --region {} | python3 $vcflib_dir/../scripts/vcffirstheader | $vcflib_dir/vcfstreamsort -w 1000 | $vcflib_dir/vcfuniq > ${parent1_tag}-ref.ref.caller.raw.vcf 

# cat ref.genome.raw.fa.fai |cut -f 1 | $parallel_dir/parallel -I% --max-args 1 -k -j $threads $java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=1 -jar $gatk4_dir/gatk4.jar HaplotypeCaller -R ref.genome.raw.fa -I ${parent1_tag}-ref.ref.realn.bam  -ploidy $ploidy -L % -O ${parent1_tag}-ref.ref.caller.raw.by_chr.%.vcf

# cat ${parent1_tag}-ref.ref.caller.raw.by_chr.*.vcf | python3 $vcflib_dir/../scripts/vcffirstheader | $vcflib_dir/vcfstreamsort -w 1000 | $vcflib_dir/vcfuniq > ${parent1_tag}-ref.ref.caller.raw.unsort.vcf

# java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar SortVcf  -SEQUENCE_DICTIONARY ref.genome.raw.dict -I ${parent1_tag}-ref.ref.caller.raw.unsort.vcf  -O  ${parent1_tag}-ref.ref.caller.raw.vcf

# if [[ $debug == "no" ]]
# then
#     rm ${parent1_tag}-ref.ref.caller.raw.by_chr.*.vcf 
#     rm ${parent1_tag}-ref.ref.caller.raw.unsort.vcf
#     rm ${parent1_tag}-ref.ref.caller.raw.*.idx
# fi

$vt_dir/vt decompose_blocksub ${parent1_tag}-ref.ref.caller.raw.vcf -a -o ${parent1_tag}-ref.ref.caller.decompose.vcf
$vt_dir/vt normalize ${parent1_tag}-ref.ref.caller.decompose.vcf -r ref.genome.raw.fa -f "~VARIANT_CONTAINS_N"| $vt_dir/vt uniq - -o ${parent1_tag}-ref.ref.caller.normalize.vcf
$vt_dir/vt annotate_variants ${parent1_tag}-ref.ref.caller.normalize.vcf -r ref.genome.raw.fa -o ${parent1_tag}-ref.ref.caller.annotate.vcf

cat ${parent1_tag}-ref.ref.caller.annotate.vcf | egrep "^#" > ${parent1_tag}-ref.ref.caller.annotate.vcf.header
$bedtools_dir/bedtools subtract -a ${parent1_tag}-ref.ref.caller.annotate.vcf -b ./../$reference_hardmask_bed > ${parent1_tag}-ref.ref.caller.annotate_hardmask.vcf.content
cat ${parent1_tag}-ref.ref.caller.annotate.vcf.header ${parent1_tag}-ref.ref.caller.annotate_hardmask.vcf.content > ${parent1_tag}-ref.ref.caller.annotate_hardmask.vcf

parent1_based_CNV_bed_line_count=$(cat ${parent1_tag}-ref.ref.significant_CNV.bed |sed '/^\s*$/d' | wc -l)
if [[ "$parent1_based_CNV_bed_line_count" > 0 ]]
then
    $bedtools_dir/bedtools subtract -a ${parent1_tag}-ref.ref.caller.annotate_hardmask.vcf -b ${parent1_tag}-ref.ref.significant_CNV.bed > ${parent1_tag}-ref.ref.caller.annotate_hardmask_CNVmask.vcf.content
    cat ${parent1_tag}-ref.ref.caller.annotate.vcf.header ${parent1_tag}-ref.ref.caller.annotate_hardmask_CNVmask.vcf.content > ${parent1_tag}-ref.ref.caller.annotate_hardmask_CNVmask.vcf
    rm ${parent1_tag}-ref.ref.caller.annotate.vcf.header
    rm ${parent1_tag}-ref.ref.caller.annotate_hardmask.vcf.content
    rm ${parent1_tag}-ref.ref.caller.annotate_hardmask_CNVmask.vcf.content
else
    cp ${parent1_tag}-ref.ref.caller.annotate_hardmask.vcf ${parent1_tag}-ref.ref.caller.annotate_hardmask_CNVmask.vcf
    rm ${parent1_tag}-ref.ref.caller.annotate.vcf.header
    rm ${parent1_tag}-ref.ref.caller.annotate_hardmask.vcf.content
fi

perl $RECOMBINEX_HOME/scripts/filter_vcf_by_window.pl \
    -i ${parent1_tag}-ref.ref.caller.annotate_hardmask_CNVmask.vcf \
    -w $cluster_window_size \
    -o ${parent1_tag}-ref.ref.caller.annotate_hardmask_CNVmask.thin.vcf

$vt_dir/vt view ${parent1_tag}-ref.ref.caller.annotate_hardmask_CNVmask.thin.vcf -f "VTYPE==SNP" -o ${parent1_tag}-ref.ref.caller.annotate_hardmask_CNVmask.thin.SNP.vcf
$vt_dir/vt view ${parent1_tag}-ref.ref.caller.annotate_hardmask_CNVmask.thin.vcf -f "VTYPE==INDEL" -o ${parent1_tag}-ref.ref.caller.annotate_hardmask_CNVmask.thin.INDEL.vcf

$vcflib_dir/vcffilter -f "QUAL > $variant_calling_quality_cutoff & QUAL / AO > 1 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ${parent1_tag}-ref.ref.caller.annotate_hardmask_CNVmask.thin.SNP.vcf |egrep "^#|AC=1" \
    > ${parent1_tag}-ref.ref.read_mapping.SNP.filter.vcf
$vcflib_dir/vcffilter -f "QUAL > $variant_calling_quality_cutoff & QUAL / AO > 1 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ${parent1_tag}-ref.ref.caller.annotate_hardmask_CNVmask.thin.INDEL.vcf |egrep "^#|AC=1" \
    > ${parent1_tag}-ref.ref.read_mapping.INDEL.filter.vcf

# compress vcf files
gzip *.vcf 

if [[ $debug == "no" ]]
then
    rm ${parent1_tag}-ref.ref.caller.raw.vcf.gz
    rm ${parent1_tag}-ref.ref.caller.decompose.vcf.gz
    rm ${parent1_tag}-ref.ref.caller.normalize.vcf.gz
    # rm ${parent1_tag}-ref.ref.caller.annotate.vcf.gz
    # rm ${parent1_tag}-ref.ref.caller.annotate_hardmask.vcf.gz
    # rm ${parent1_tag}-ref.ref.caller.annotate_hardmask_CNVmask.vcf.gz
    # rm ${parent1_tag}-ref.ref.caller.annotate_hardmask_CNVmask.thin.vcf.gz
fi

# trim the reads
$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $trimmomatic_dir/trimmomatic.jar PE -threads $threads -phred33 $parent2_tag.R1.raw.fq.gz $parent2_tag.R2.raw.fq.gz $parent2_tag.R1.trimmed.PE.fq.gz $parent2_tag.R1.trimmed.SE.fq.gz $parent2_tag.R2.trimmed.PE.fq.gz $parent2_tag.R2.trimmed.SE.fq.gz ILLUMINACLIP:adapter.fa:2:30:10  SLIDINGWINDOW:5:20 MINLEN:36

if [[ $debug == "no" ]]
then
    rm $parent2_tag.R1.trimmed.SE.fq.gz
    rm $parent2_tag.R2.trimmed.SE.fq.gz
    rm adapter.fa
fi

# map reads to the reference genome
$bwa_dir/bwa mem -t $threads -M ref.genome.raw.fa  $parent2_tag.R1.trimmed.PE.fq.gz $parent2_tag.R2.trimmed.PE.fq.gz | $samtools_dir/samtools view -bS -q $mapping_quality_cutoff - >${parent2_tag}-ref.ref.bam

if [[ $debug == "no" ]]
then
    rm ref.genome.raw.fa.bwt
    rm ref.genome.raw.fa.pac
    rm ref.genome.raw.fa.ann
    rm ref.genome.raw.fa.amb
    rm ref.genome.raw.fa.sa
    rm $parent2_tag.R1.trimmed.PE.fq.gz
    rm $parent2_tag.R2.trimmed.PE.fq.gz
fi

# sort bam file by picard-tools SortSam
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar SortSam \
    -INPUT ${parent2_tag}-ref.ref.bam \
    -OUTPUT ${parent2_tag}-ref.ref.sort.bam \
    -SORT_ORDER coordinate

if [[ $debug == "no" ]]
then
    rm ${parent2_tag}-ref.ref.bam
fi

# fixmate
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar FixMateInformation \
    -INPUT ${parent2_tag}-ref.ref.sort.bam \
    -OUTPUT ${parent2_tag}-ref.ref.fixmate.bam

if [[ $debug == "no" ]]
then
    rm ${parent2_tag}-ref.ref.sort.bam
fi

# add or replace read groups and sort
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar AddOrReplaceReadGroups \
    -INPUT ${parent2_tag}-ref.ref.fixmate.bam \
    -OUTPUT ${parent2_tag}-ref.ref.rdgrp.bam \
    -SORT_ORDER coordinate \
    -RGID "$parent2_tag" \
    -RGLB "$parent2_tag" \
    -RGPL "Illumina" \
    -RGPU "$parent2_tag" \
    -RGSM "$parent2_tag" \
    -RGCN "RGCN" 

if [[ $debug == "no" ]]
then
    rm ${parent2_tag}-ref.ref.fixmate.bam
fi

# remove duplicates
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar MarkDuplicates \
    -INPUT ${parent2_tag}-ref.ref.rdgrp.bam \
    -REMOVE_DUPLICATES true  \
    -METRICS_FILE ${parent2_tag}-ref.ref.dedup.matrics \
    -OUTPUT ${parent2_tag}-ref.ref.dedup.bam 

if [[ $debug == "no" ]]
then
    rm ${parent2_tag}-ref.ref.rdgrp.bam
fi

# index the dedup.bam file
$samtools_dir/samtools index ${parent2_tag}-ref.ref.dedup.bam

# GATK local realign
# find realigner targets
$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $gatk3_dir/gatk3.jar \
    -nt $threads \
    -R ref.genome.raw.fa \
    -T RealignerTargetCreator \
    -I ${parent2_tag}-ref.ref.dedup.bam  \
    -o ${parent2_tag}-ref.ref.realn.intervals 
# run realigner
$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $gatk3_dir/gatk3.jar  \
    -R ref.genome.raw.fa \
    -T IndelRealigner \
    -I ${parent2_tag}-ref.ref.dedup.bam \
    -targetIntervals ${parent2_tag}-ref.ref.realn.intervals \
    -o ${parent2_tag}-ref.ref.realn.bam

if [[ $debug == "no" ]]
then
    rm ${parent2_tag}-ref.ref.dedup.bam
    rm ${parent2_tag}-ref.ref.dedup.bam.bai
    rm ${parent2_tag}-ref.ref.dedup.matrics
    rm ${parent2_tag}-ref.ref.realn.intervals
fi

# generate samtools mpileup
$samtools_dir/samtools mpileup -C 0 -q $mapping_quality_cutoff -f ref.genome.raw.fa  ${parent2_tag}-ref.ref.realn.bam |gzip -c >${parent2_tag}-ref.ref.mpileup.gz

# compute basic alignment statistics by samtools
$samtools_dir/samtools flagstat ${parent2_tag}-ref.ref.realn.bam >${parent2_tag}-ref.ref.samstat

# calculate per-base depth
$samtools_dir/samtools depth -aa ${parent2_tag}-ref.ref.realn.bam |gzip -c >${parent2_tag}-ref.ref.depth.txt.gz

# compute insert size statistics
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar CollectInsertSizeMetrics \
    -I ${parent2_tag}-ref.ref.realn.bam \
    -O ${parent2_tag}-ref.ref.insert_size_metrics.txt \
    -H ${parent2_tag}-ref.ref.insert_size_histogram.pdf \
    -M 0.5

# calculate read mapping coverage statistics
perl $RECOMBINEX_HOME/scripts/summarize_mapping_coverage.pl \
    -r ref.genome.raw.fa \
    -s ${parent2_tag}-ref.ref.samstat \
    -d ${parent2_tag}-ref.ref.depth.txt.gz \
    -c 5 \
    -t $parent2_tag \
    -o ${parent2_tag}-ref.ref.coverage_summary.txt

# scan for CNV by FREEC
#parent2_raw_read_length=$(gunzip -c $parent2_tag.R1.raw.fq.gz |awk 'NR%4==2{print length($0)}' | head -1 || true)
parent2_raw_read_length=100
if [ -z "$excluded_chr_list_for_cnv_profiling" ]
then
    perl $RECOMBINEX_HOME/scripts/run_FREEC_wrapper_lite_for_RecombineX.pl \
	-r ref \
	-bam ${parent2_tag}-ref.ref.realn.bam \
	-prefix ${parent2_tag}-ref.ref \
	-ploidy $ploidy \
	-bedtools $bedtools_dir/bedtools \
	-samtools $samtools_dir/samtools \
	-freec $freec_dir/freec \
	-window $window_size \
	-step $step_size \
	-read_length_for_mappability $parent2_raw_read_length \
	-min_mappability $min_mappability \
	-min_expected_gc $min_expected_gc \
	-max_expected_gc $max_expected_gc \
	-mates_orientation 0 \
	-refseq_genome_preprocessing_dir "./../$reference_genome_preprocessing_dir" \
	-threads $threads 
else
    perl $RECOMBINEX_HOME/scripts/run_FREEC_wrapper_lite_for_RecombineX.pl \
	-r ref.genome.raw.fa \
	-bam ${parent2_tag}-ref.ref.realn.bam \
	-prefix ${parent2_tag}-ref.ref \
	-ploidy $ploidy \
	-bedtools $bedtools_dir/bedtools \
	-samtools $samtools_dir/samtools \
	-freec $freec_dir/freec \
	-window $window_size \
	-step $step_size \
	-read_length_for_mappability $parent2_raw_read_length \
	-min_mappability $min_mappability \
	-min_expected_gc $min_expected_gc \
	-max_expected_gc $max_expected_gc \
	-mates_orientation 0 \
	-threads $threads \
	-refseq_genome_preprocessing_dir "./../$reference_genome_preprocessing_dir" \
	-excluded_chr_list ./../$excluded_chr_list_for_cnv_profiling
fi
mv FREEC.config.txt ${parent2_tag}-ref.ref.FREEC.config.txt

Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/CNV_segmentation_by_DNAcopy.R \
    --input ${parent2_tag}-ref.ref.FREEC.bam_ratio.txt \
    --prefix ${parent2_tag}-ref.ref \
    --window $window_size \
    --ploidy $ploidy \
    --genome_fai ref.genome.raw.fa.fai

perl $RECOMBINEX_HOME/scripts/adjust_FREEC_copynumber_by_DNAcopy_copynumber.pl \
    -i ${parent2_tag}-ref.ref.FREEC.bam_ratio.sorted.txt \
    -a ${parent2_tag}-ref.ref.FREEC.bam_ratio.sorted.resegmented.lite.txt \
    -o ${parent2_tag}-ref.ref.FREEC.bam_ratio.sorted.adjusted.txt

if [[ -s ${parent2_tag}-ref.ref.FREEC.bam_ratio.sorted.adjusted.txt ]]
then
    Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_CNV_for_FREEC.R \
	--ploidy $ploidy \
	--genome_fai ref.genome.raw.fa.fai \
	--input ${parent2_tag}-ref.ref.FREEC.bam_ratio.sorted.adjusted.txt \
	--output ${parent2_tag}-ref.ref.CNV_plot.pdf
    rm Rplots.pdf
    if [[ -s ${parent2_tag}-ref.ref.FREEC.bam_ratio.sorted.resegmented.CNVs.txt ]]
    then
        Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/assess_CNV_significance_for_FREEC.R \
	    --cnv ${parent2_tag}-ref.ref.FREEC.bam_ratio.sorted.resegmented.CNVs.txt \
	    --ratio ${parent2_tag}-ref.ref.FREEC.bam_ratio.sorted.adjusted.txt \
	    --genome_fai ref.genome.raw.fa.fai \
	    --output ${parent2_tag}-ref.ref.CNV_significance_test.txt
	cat ${parent2_tag}-ref.ref.CNV_significance_test.txt | tail -n +2 | awk '{ if ($6 < 0.05) print $0 }' | awk '{print $1, $2-1, $3}' OFS='\t' > ${parent2_tag}-ref.ref.significant_CNV.bed
    else
        echo -e "chr\tstart\tend\tcopy_number\tstatus\tMWU_test_p_value\tKS_test_p_value" > ${parent2_tag}-ref.ref.CNV_significance_test.txt
        echo -e "\t\t" > ${parent2_tag}-ref.ref.significant_CNV.bed
        echo 
    fi
else
    echo "Exception encountered for FREEC! Exit! ..." > ${parent2_tag}-ref.ref.realn.bam.no_FREEC.txt
    exit
fi

mkdir ${parent2_tag}-ref.ref.FREEC_intermediate_files
mv for_CNV.bam ${parent2_tag}-ref.ref.FREEC_intermediate_files
mv GC_profile.${window_size}bp.cnp ${parent2_tag}-ref.ref.FREEC_intermediate_files
mv FREEC.bam* ${parent2_tag}-ref.ref.FREEC_intermediate_files

if [[ $debug == "no" ]]
then
    rm -r ${parent2_tag}-ref.ref.FREEC_intermediate_files
    rm $parent2_tag.R1.raw.fq.gz
    rm $parent2_tag.R2.raw.fq.gz
fi

# SNP and INDEL calling

python3 $freebayes_dir/../scripts/fasta_generate_regions.py ref.genome.raw.fa.fai 100000 > ref.genome.region.txt
cat ref.genome.region.txt | $parallel_dir/parallel -k -j $threads $freebayes_dir/freebayes -f ref.genome.raw.fa  -p $ploidy ${parent2_tag}-ref.ref.realn.bam --region {} | python3 $vcflib_dir/../scripts/vcffirstheader | $vcflib_dir/vcfstreamsort -w 1000 | $vcflib_dir/vcfuniq > ${parent2_tag}-ref.ref.caller.raw.vcf 

# cat ref.genome.raw.fa.fai |cut -f 1 | $parallel_dir/parallel -I% --max-args 1 -k -j $threads $java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=1 -jar $gatk4_dir/gatk4.jar HaplotypeCaller -R ref.genome.raw.fa -I ${parent2_tag}-ref.ref.realn.bam  -ploidy $ploidy -L % -O ${parent2_tag}-ref.ref.caller.raw.by_chr.%.vcf
# cat ${parent2_tag}-ref.ref.caller.raw.by_chr.*.vcf | python3 $vcflib_dir/../scripts/vcffirstheader | $vcflib_dir/vcfstreamsort -w 1000 | $vcflib_dir/vcfuniq > ${parent2_tag}-ref.ref.caller.raw.unsort.vcf
# java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar SortVcf  -SEQUENCE_DICTIONARY ref.genome.raw.dict -I ${parent2_tag}-ref.ref.caller.raw.unsort.vcf  -O  ${parent2_tag}-ref.ref.caller.raw.vcf

# if [[ $debug == "no" ]]
# then
#     rm ${parent2_tag}-ref.ref.caller.raw.by_chr.*.vcf 
#     rm ${parent2_tag}-ref.ref.caller.raw.unsort.vcf
#     rm ${parent2_tag}-ref.ref.caller.raw.*.idx
# fi

$vt_dir/vt decompose_blocksub ${parent2_tag}-ref.ref.caller.raw.vcf -a -o ${parent2_tag}-ref.ref.caller.decompose.vcf

$vt_dir/vt normalize ${parent2_tag}-ref.ref.caller.decompose.vcf -r ref.genome.raw.fa -f "~VARIANT_CONTAINS_N"| $vt_dir/vt uniq - -o ${parent2_tag}-ref.ref.caller.normalize.vcf
$vt_dir/vt annotate_variants ${parent2_tag}-ref.ref.caller.normalize.vcf -r ref.genome.raw.fa -o ${parent2_tag}-ref.ref.caller.annotate.vcf

cat ${parent2_tag}-ref.ref.caller.annotate.vcf | egrep "^#" > ${parent2_tag}-ref.ref.caller.annotate.vcf.header
$bedtools_dir/bedtools subtract -a ${parent2_tag}-ref.ref.caller.annotate.vcf -b ./../$reference_hardmask_bed > ${parent2_tag}-ref.ref.caller.annotate_hardmask.vcf.content
cat ${parent2_tag}-ref.ref.caller.annotate.vcf.header ${parent2_tag}-ref.ref.caller.annotate_hardmask.vcf.content > ${parent2_tag}-ref.ref.caller.annotate_hardmask.vcf

parent2_based_CNV_bed_line_count=$(cat ${parent2_tag}-ref.ref.significant_CNV.bed |sed '/^\s*$/d' | wc -l)
if [[ "$parent2_based_CNV_bed_line_count" > 0 ]]
then
    $bedtools_dir/bedtools subtract -a ${parent2_tag}-ref.ref.caller.annotate_hardmask.vcf -b ${parent2_tag}-ref.ref.significant_CNV.bed > ${parent2_tag}-ref.ref.caller.annotate_hardmask_CNVmask.vcf.content
    cat ${parent2_tag}-ref.ref.caller.annotate.vcf.header ${parent2_tag}-ref.ref.caller.annotate_hardmask_CNVmask.vcf.content >${parent2_tag}-ref.ref.caller.annotate_hardmask_CNVmask.vcf
    rm ${parent2_tag}-ref.ref.caller.annotate.vcf.header
    rm ${parent2_tag}-ref.ref.caller.annotate_hardmask.vcf.content
    rm ${parent2_tag}-ref.ref.caller.annotate_hardmask_CNVmask.vcf.content
else
    cp ${parent2_tag}-ref.ref.caller.annotate_hardmask.vcf ${parent2_tag}-ref.ref.caller.annotate_hardmask_CNVmask.vcf
    rm ${parent2_tag}-ref.ref.caller.annotate.vcf.header
    rm ${parent2_tag}-ref.ref.caller.annotate_hardmask.vcf.content
fi

perl $RECOMBINEX_HOME/scripts/filter_vcf_by_window.pl \
    -i ${parent2_tag}-ref.ref.caller.annotate_hardmask_CNVmask.vcf \
    -w $cluster_window_size \
    -o ${parent2_tag}-ref.ref.caller.annotate_hardmask_CNVmask.thin.vcf

$vt_dir/vt view ${parent2_tag}-ref.ref.caller.annotate_hardmask_CNVmask.thin.vcf -f "VTYPE==SNP" -o ${parent2_tag}-ref.ref.caller.annotate_hardmask_CNVmask.thin.SNP.vcf
$vt_dir/vt view ${parent2_tag}-ref.ref.caller.annotate_hardmask_CNVmask.thin.vcf -f "VTYPE==INDEL" -o ${parent2_tag}-ref.ref.caller.annotate_hardmask_CNVmask.thin.INDEL.vcf

$vcflib_dir/vcffilter -f "QUAL > $variant_calling_quality_cutoff & QUAL / AO > 1 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ${parent2_tag}-ref.ref.caller.annotate_hardmask_CNVmask.thin.SNP.vcf |egrep "^#|AC=1" \
    > ${parent2_tag}-ref.ref.read_mapping.SNP.filter.vcf
$vcflib_dir/vcffilter -f "QUAL > $variant_calling_quality_cutoff & QUAL / AO > 1 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ${parent2_tag}-ref.ref.caller.annotate_hardmask_CNVmask.thin.INDEL.vcf |egrep "^#|AC=1" \
    > ${parent2_tag}-ref.ref.read_mapping.INDEL.filter.vcf

# compress vcf files
gzip *.vcf 

rm -r tmp

if [[ $debug == "no" ]]
then
    rm ${parent2_tag}-ref.ref.caller.raw.vcf.gz
    rm ${parent2_tag}-ref.ref.caller.decompose.vcf.gz
    rm ${parent2_tag}-ref.ref.caller.normalize.vcf.gz
    # rm ${parent2_tag}-ref.ref.caller.annotate.vcf.gz
    # rm ${parent2_tag}-ref.ref.caller.annotate_hardmask.vcf.gz
    # rm ${parent2_tag}-ref.ref.caller.annotate_hardmask_CNVmask.vcf.gz
    # rm ${parent2_tag}-ref.ref.caller.annotate_hardmask_CNVmask.thin.vcf.gz
fi

# extract markers
reference_based_prefix="${parent1_tag}-${parent2_tag}.ref"

perl $RECOMBINEX_HOME/scripts/extract_polymorphic_markers_by_reference_based_vcf.pl \
    -r ref.genome.raw.fa \
    -v1 ${parent1_tag}-ref.ref.read_mapping.SNP.filter.vcf.gz \
    -v2 ${parent2_tag}-ref.ref.read_mapping.SNP.filter.vcf.gz \
    -b0 ./../$reference_hardmask_bed \
    -b1 ${parent1_tag}-ref.ref.significant_CNV.bed \
    -b2 ${parent2_tag}-ref.ref.significant_CNV.bed \
    -q $variant_calling_quality_cutoff \
    -o $reference_based_prefix.depth_filter.SNP.markers.txt.gz


# perl $RECOMBINEX_HOME/scripts/extract_polymorphic_markers_by_reference_based_vcf.pl \
#     -r ref.genome.raw.fa \
#     -v1 ${parent1_tag}-ref.ref.read_mapping.INDEL.filter.vcf.gz \
#     -v2 ${parent2_tag}-ref.ref.read_mapping.INDEL.filter.vcf.gz \
#     -b0 ./../$reference_hardmask_bed \
#     -b1 ${parent1_tag}-ref.ref.significant_CNV.bed \
#     -b2 ${parent2_tag}-ref.ref.significant_CNV.bed \
#     -q $variant_calling_quality_cutoff \
#     -o $reference_based_prefix.depth_filter.INDEL.markers.txt.gz


perl $RECOMBINEX_HOME/scripts/filter_reference_based_markers_by_mpileup_and_depth.pl \
     -i $reference_based_prefix.depth_filter.SNP.markers.txt.gz \
     -ref_depth_summary ${parent1_tag}-ref.ref.coverage_summary.txt \
     -query_depth_summary ${parent2_tag}-ref.ref.coverage_summary.txt \
     -ref_mpileup ${parent1_tag}-ref.ref.mpileup.gz \
     -query_mpileup ${parent2_tag}-ref.ref.mpileup.gz \
     -o $reference_based_prefix.final.SNP.markers.txt.gz

# perl $RECOMBINEX_HOME/scripts/filter_reference_based_markers_by_mpileup_and_depth.pl \
#      -i $reference_based_prefix.depth_filter.INDEL.markers.txt.gz \
#      -ref_depth_summary ${parent1_tag}-ref.ref.coverage_summary.txt \
#      -query_depth_summary ${parent2_tag}-ref.ref.coverage_summary.txt \
#      -ref_mpileup ${parent1_tag}-ref.ref.mpileup.gz \
#      -query_mpileup ${parent2_tag}-ref.ref.mpileup.gz \
#      -o $reference_based_prefix.final.INDEL.markers.txt.gz


# copy the final marker vcf files back to the task-specific directory
cp $reference_based_prefix.final.SNP.markers.txt.gz ./../
# cp $reference_based_prefix.final.INDEL.markers.txt.gz ./../

cd ..

# calculate intermarker distance

perl $RECOMBINEX_HOME/scripts/intermarker_distance.pl \
    -i $reference_based_prefix.final.SNP.markers.txt.gz \
    -p $reference_based_prefix.final.SNP.markers

# perl $RECOMBINEX_HOME/scripts/intermarker_distance.pl \
#     -i $reference_based_prefix.final.INDEL.markers.txt.gz \
#     -p $reference_based_prefix.final.INDEL.markers
# gunzip -c $reference_based_prefix.final.INDEL.markers.txt.gz | tail -n +2 | gzip > $reference_based_prefix.final.INDEL.markers.content.txt.gz
# cat $reference_based_prefix.final.SNP.markers.txt.gz $reference_based_prefix.final.INDEL.markers.content.txt.gz > $reference_based_prefix.final.BOTH.markers.txt.gz
# rm $reference_based_prefix.final.INDEL.markers.content.txt.gz
# perl $RECOMBINEX_HOME/scripts/intermarker_distance.pl \
#     -i $reference_based_prefix.final.BOTH.markers.txt.gz \
#     -p $reference_based_prefix.final.BOTH.markers

reference_fai="./$reference_based_output_dir/ref.genome.raw.fa.fai"

# plot markers
if [[ $use_centromere_annotation == "yes" ]]
then
    Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R \
	--genome_fai $reference_fai \
	--marker_table $reference_based_prefix.final.SNP.markers.txt.gz \
	--centromere_gff $reference_centromere_gff \
	--output_prefix $reference_based_prefix.final.SNP.markers
    rm Rplots.pdf
    # Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R \
    # 	--genome_fai $reference_fai \
    # 	--marker_table $reference_based_prefix.final.INDEL.markers.txt.gz \
    # 	--centromere_gff $reference_centromere_gff \
    # 	--output_prefix $reference_based_prefix.final.INDEL.markers
    # rm Rplots.pdf
    # Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R \
    # 	--genome_fai $reference_fai \
    # 	--marker_table $reference_based_prefix.final.BOTH.markers.txt.gz \
    # 	--centromere_gff $reference_centromere_gff \
    # 	--output_prefix $reference_based_prefix.final.BOTH.markers
    # rm Rplots.pdf
else
    Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R \
	--genome_fai $reference_fai \
	--marker_table $reference_based_prefix.final.SNP.markers.txt.gz \
	--output_prefix $reference_based_prefix.final.SNP.markers
    rm Rplots.pdf
    # Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $reference_fai \
    # 	    --marker_table $reference_based_prefix.final.INDEL.markers.txt.gz \
    # 	    --output_prefix $reference_based_prefix.final.INDEL.markers
    # rm Rplots.pdf
    # Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $reference_fai \
    # 	    --marker_table $reference_based_prefix.final.BOTH.markers.txt.gz \
    # 	    --output_prefix $reference_based_prefix.final.BOTH.markers
    # rm Rplots.pdf

fi




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


