#!/bin/bash
set -e -o pipefail

##########################################
# load environment variables for RecombineX
source ./../../env.sh

###########################################
# set project-specific variables
parent1_tag="S288C" # The unique tag of parental genome 1 used in 11.Parent_Genome_Preprocessing. Default = "S288C".
parent2_tag="SK1" # The unique tag of parental genome 2 used in 11.Parent_Genome_Preprocessing. Default = "SK1".
parent1_reads_R1="./../00.Parent_Reads/$parent1_tag.R1.fq.gz" # The path to the paired-end R1 reads of the parental genome 1. Default = "./../00.Parent_Reads/$parent1_tag.R1.fq.gz".
parent1_reads_R2="./../00.Parent_Reads/$parent1_tag.R2.fq.gz" # The path to the paired-end R2 reads of the parental genome 1. Default = "./../00.Parent_Reads/$parent1_tag.R2.fq.gz".
parent2_reads_R1="./../00.Parent_Reads/$parent2_tag.R1.fq.gz" # The path to the paired-end R1 reads of the parental genome 2. Default = "./../00.Parent_Reads/$parent2_tag.R1.fq.gz".
parent2_reads_R2="./../00.Parent_Reads/$parent2_tag.R2.fq.gz" # The path to the paired-end R2 reads of the parental genome 2. Default = "./../00.Parent_Reads/$parent2_tag.R2.fq.gz".
ploidy=1 # The ploidy of the parental genome: "1" for haploid and "2" for diploid. Default = "2". If the parental genome is purely homozygous, it is recommended to set "ploidy=1" to maximize the power of CNV profiling. Default= "1".
window_size=250 # The window size for non-overlapping sliding-window-based CNV profiling. Default = "250" (i.e. 250 bp). 
apply_cnv_filter="yes". # Whether to apply the CNV filter for marker candidates. Set this option to "yes" when running RecombineX for the mitochondrial genome. Default = "yes". 
threads=4 # The number of threads to use. Default = "4".
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
##########################################




# process the pipeline
##########################################
# Normally no need to change the following parameters
parent_genome_preprocessing_dir="./../11.Parent_Genome_Preprocessing" # The relative path to 11.Parent_Genome_Preprocessing directory.
parent1_raw_assembly="$parent_genome_preprocessing_dir/$parent1_tag.genome.raw.relabel.fa" # The relabeled genome 1 raw assembly file generated in 11.Parent_Genome_Preprocessing directory.
parent2_raw_assembly="$parent_genome_preprocessing_dir/$parent2_tag.genome.raw.relabel.fa" # The relabeled genome 2 raw assembly file generated in 11.Parent_Genome_Preprocessing directory.
parent1_hardmask_bed="$parent_genome_preprocessing_dir/$parent1_tag.genome.hardmask.relabel.masking_details.bed" # The masking details bed file for the relabeled and hardmasked genome 1 assembly file generated in 11.Parent_Genome_Preprocessing directory.
parent2_hardmask_bed="$parent_genome_preprocessing_dir/$parent2_tag.genome.hardmask.relabel.masking_details.bed" # The masking details bed file for the relabeled and hardmasked genome 2 assembly file generated in 11.Parent_Genome_Preprocessing directory.
excluded_chr_list_for_cnv_profiling="" # The relative path to the list for specifying chromosomes/scaffolds/contigs in relabeled parental genomes (in ./../11.Parent_Genome_Preprocessing) to be exclued for CNV profiling. Default = "".
mapping_quality_cutoff=30 # The minimal mapping quality to be considered. Default = "30".
variant_calling_quality_cutoff=30 # The minimal variant calling quality to be considered. Default = "30".

min_mappability=0.85 # The minimal mappability for sliding-window-based CNV profiling. Default = "0.85". 
cluster_window_size=10 # Adjacent variants within the specified window (unit: bp) will be both filtered out if any of them is INDEL. Default = "10".
###############################################


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
echo "test the existence of parent1_raw_assembly: $parent1_raw_assembly"
test_file_existence $parent1_raw_assembly
echo ""
echo "test the existence of parent1_hardmask_bed: $parent1_hardmask_bed"
test_file_existence $parent1_hardmask_bed
echo ""
echo "test the existence of parent2_raw_assembly: $parent2_raw_assembly"
test_file_existence $parent2_raw_assembly
echo ""
echo "test the existence of parent2_hardmask_bed: $parent2_hardmask_bed"
test_file_existence $parent2_hardmask_bed

echo ""
echo "test the existence of parent1_reads_R1: $parent1_reads_R1"
test_file_existence $parent1_reads_R1
echo ""
echo "test the existence of parent1_reads_R2: $parent1_reads_R2"
test_file_existence $parent1_reads_R2
echo ""
echo "test the existence of parent2_reads_R1: $parent2_reads_R1"
test_file_existence $parent2_reads_R1
echo ""
echo "test the existence of parent2_reads_R2: $parent2_reads_R2"
test_file_existence $parent2_reads_R2
echo ""



adapter="$trimmomatic_dir/adapters/TruSeq3-PE-2.fa"
parent1_based_output_dir="${parent1_tag}-${parent2_tag}_${parent1_tag}_based"
parent2_based_output_dir="${parent1_tag}-${parent2_tag}_${parent2_tag}_based"
mkdir $parent1_based_output_dir
mkdir $parent2_based_output_dir

parent1_based_prefix="${parent1_tag}-${parent2_tag}.${parent1_tag}"
parent2_based_prefix="${parent1_tag}-${parent2_tag}.${parent2_tag}"

cd $parent1_based_output_dir

ln -s ./../$parent1_raw_assembly $parent1_tag.genome.raw.fa
ln -s ./../$parent2_reads_R1 $parent2_tag.R1.raw.fq.gz
ln -s ./../$parent2_reads_R2 $parent2_tag.R2.raw.fq.gz

cp $adapter adapter.fa

mkdir tmp
$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $trimmomatic_dir/trimmomatic.jar PE -threads $threads -phred33 $parent2_tag.R1.raw.fq.gz $parent2_tag.R2.raw.fq.gz $parent2_tag.R1.trimmed.PE.fq.gz $parent2_tag.R1.trimmed.SE.fq.gz $parent2_tag.R2.trimmed.PE.fq.gz $parent2_tag.R2.trimmed.SE.fq.gz ILLUMINACLIP:adapter.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:36

if [[ $debug == "no" ]]
then
    rm $parent2_tag.R1.trimmed.SE.fq.gz
    rm $parent2_tag.R2.trimmed.SE.fq.gz
    rm adapter.fa
fi

$bwa_dir/bwa index $parent1_tag.genome.raw.fa
$bwa_dir/bwa mem -t $threads -M $parent1_tag.genome.raw.fa  $parent2_tag.R1.trimmed.PE.fq.gz $parent2_tag.R2.trimmed.PE.fq.gz |$samtools_dir/samtools view -bS -q $mapping_quality_cutoff -F 3340 -f 2 - >$parent1_based_prefix.bam

if [[ $debug == "no" ]]
then
    rm $parent1_tag.genome.raw.fa.bwt
    rm $parent1_tag.genome.raw.fa.pac
    rm $parent1_tag.genome.raw.fa.ann
    rm $parent1_tag.genome.raw.fa.amb
    rm $parent1_tag.genome.raw.fa.sa
    rm $parent2_tag.R1.trimmed.PE.fq.gz
    rm $parent2_tag.R2.trimmed.PE.fq.gz
fi

# index reference sequence
$samtools_dir/samtools faidx $parent1_tag.genome.raw.fa
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar CreateSequenceDictionary \
    -REFERENCE $parent1_tag.genome.raw.fa \
    -OUTPUT $parent1_tag.genome.raw.dict


# sort bam file by picard-tools SortSam
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar SortSam \
    -INPUT $parent1_based_prefix.bam \
    -OUTPUT $parent1_based_prefix.sort.bam \
    -SORT_ORDER coordinate

if [[ $debug == "no" ]]
then
    rm $parent1_based_prefix.bam
fi

# fixmate
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar FixMateInformation \
    -INPUT $parent1_based_prefix.sort.bam \
    -OUTPUT $parent1_based_prefix.fixmate.bam

if [[ $debug == "no" ]]
then
    rm $parent1_based_prefix.sort.bam
fi

# add or replace read groups and sort
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar AddOrReplaceReadGroups \
    -INPUT $parent1_based_prefix.fixmate.bam \
    -OUTPUT $parent1_based_prefix.rdgrp.bam \
    -SORT_ORDER coordinate \
    -RGID "$parent1_based_prefix" \
    -RGLB "$parent1_based_prefix" \
    -RGPL "Illumina" \
    -RGPU "$parent1_based_prefix" \
    -RGSM "$parent1_based_prefix" \
    -RGCN "RGCN" 

if [[ $debug == "no" ]]
then
    rm $parent1_based_prefix.fixmate.bam
fi

# remove duplicates
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar MarkDuplicates \
    -INPUT $parent1_based_prefix.rdgrp.bam \
    -REMOVE_DUPLICATES true  \
    -METRICS_FILE $parent1_based_prefix.dedup.matrics \
    -OUTPUT $parent1_based_prefix.dedup.bam 

if [[ $debug == "no" ]]
then
    rm $parent1_based_prefix.rdgrp.bam
fi

# index the dedup.bam file
$samtools_dir/samtools index $parent1_based_prefix.dedup.bam

# GATK local realign
# find realigner targets
$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $gatk3_dir/gatk3.jar \
    -nt $threads \
    -R $parent1_tag.genome.raw.fa \
    -T RealignerTargetCreator \
    -I $parent1_based_prefix.dedup.bam \
    -o $parent1_based_prefix.realn.intervals 
# run realigner
$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $gatk3_dir/gatk3.jar  \
    -R $parent1_tag.genome.raw.fa \
    -T IndelRealigner \
    -I $parent1_based_prefix.dedup.bam \
    -targetIntervals $parent1_based_prefix.realn.intervals  \
    -o $parent1_based_prefix.realn.bam

if [[ $debug == "no" ]]
then
    rm $parent1_based_prefix.dedup.bam
    rm $parent1_based_prefix.dedup.bam.bai
    rm $parent1_based_prefix.dedup.matrics
    rm $parent1_based_prefix.realn.intervals
fi

# generate samtools mpileup
$samtools_dir/samtools mpileup -C 0 -q $mapping_quality_cutoff -f $parent1_tag.genome.raw.fa $parent1_based_prefix.realn.bam |gzip -c >$parent1_based_prefix.mpileup.gz
# calculate per-base depth
$samtools_dir/samtools depth -aa $parent1_based_prefix.realn.bam |gzip -c >$parent1_based_prefix.depth.txt.gz

# compute basic alignment statistics by samtools
$samtools_dir/samtools flagstat $parent1_based_prefix.realn.bam >$parent1_based_prefix.samstat

# compute insert size statistics
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar CollectInsertSizeMetrics \
    -I $parent1_based_prefix.realn.bam \
    -O $parent1_based_prefix.insert_size_metrics.txt \
    -H $parent1_based_prefix.insert_size_histogram.pdf \
    -M 0.5

# calculate read mapping coverage statistics
perl $RECOMBINEX_HOME/scripts/summarize_mapping_coverage.pl \
    -r $parent1_tag.genome.raw.fa \
    -s $parent1_based_prefix.samstat \
    -d $parent1_based_prefix.depth.txt.gz \
    -c 5 \
    -t $parent2_tag \
    -o $parent1_based_prefix.coverage_summary.txt


echo "mapping finished!"

# scan for CNV by FREEC
echo ""
echo "scan CNV .."
#parent2_raw_read_length=$(gunzip -c $parent2_tag.R1.raw.fq.gz |awk 'NR%4==2{print length($0)}' | head -1 || true)
parent2_raw_read_length=100
step_size=$window_size
parent1_min_expected_gc=$(cat ./../$parent_genome_preprocessing_dir/$parent1_tag.FREEC.GC_range.txt|egrep -v "#"|cut -f 1)
parent1_max_expected_gc=$(cat ./../$parent_genome_preprocessing_dir/$parent1_tag.FREEC.GC_range.txt|egrep -v "#"|cut -f 2)

echo ""
echo "window_size=$window_size, step_size=$step_size, read_length=$parent2_raw_read_length"
echo "parent1_min_expected_gc=$parent1_min_expected_gc, parent1_max_expected_gc=$parent1_max_expected_gc"
echo ""
if [ -z "$excluded_chr_list_for_cnv_profiling" ]
then
    perl $RECOMBINEX_HOME/scripts/run_FREEC_wrapper_lite_for_RecombineX.pl \
	-r $parent1_tag \
	-bam $parent1_based_prefix.realn.bam \
	-prefix $parent1_based_prefix \
	-ploidy $ploidy \
	-bedtools $bedtools_dir/bedtools \
	-samtools $samtools_dir/samtools \
	-freec $freec_dir/freec \
	-window $window_size \
	-step $step_size \
	-read_length_for_mappability $parent2_raw_read_length  \
	-min_mappability $min_mappability \
	-min_expected_gc $parent1_min_expected_gc \
	-max_expected_gc $parent1_max_expected_gc \
	-refseq_genome_preprocessing_dir "./../$parent_genome_preprocessing_dir" \
	-mates_orientation 0 \
	-threads $threads
else
    echo "excluded chr list: $excluded_chr_list_for_cnv_profiling"
    perl $RECOMBINEX_HOME/scripts/run_FREEC_wrapper_lite_for_RecombineX.pl \
	-r $parent1_tag \
	-bam $parent1_based_prefix.realn.bam \
	-prefix $parent1_based_prefix \
	-ploidy $ploidy \
	-bedtools $bedtools_dir/bedtools \
	-samtools $samtools_dir/samtools \
	-freec $freec_dir/freec \
	-window $window_size \
	-step $step_size \
	-read_length_for_mappability $parent2_raw_read_length  \
	-min_mappability $min_mappability \
	-mates_orientation 0 \
	-min_expected_gc $parent1_min_expected_gc \
	-max_expected_gc $parent1_max_expected_gc \
	-refseq_genome_preprocessing_dir "./../$parent_genome_preprocessing_dir" \
	-threads $threads \
	-excluded_chr_list ./../$excluded_chr_list_for_cnv_profiling 
fi

Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/CNV_segmentation_by_DNAcopy.R \
    --input $parent1_based_prefix.FREEC.bam_ratio.txt \
    --prefix $parent1_based_prefix \
    --window $window_size \
    --ploidy $ploidy \
    --genome_fai $parent1_tag.genome.raw.fa.fai

perl $RECOMBINEX_HOME/scripts/adjust_FREEC_copynumber_by_DNAcopy_copynumber.pl \
    -i $parent1_based_prefix.FREEC.bam_ratio.sorted.txt \
    -a $parent1_based_prefix.FREEC.bam_ratio.sorted.resegmented.lite.txt \
    -o $parent1_based_prefix.FREEC.bam_ratio.sorted.adjusted.txt

if [[ -s $parent1_based_prefix.FREEC.bam_ratio.sorted.adjusted.txt ]]
then
    Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_CNV_for_FREEC.R \
	--ploidy $ploidy \
	--genome_fai $parent1_tag.genome.raw.fa.fai \
	--input $parent1_based_prefix.FREEC.bam_ratio.sorted.adjusted.txt \
	--output $parent1_based_prefix.CNV_plot.pdf
    rm Rplots.pdf
    if [[ -s $parent1_based_prefix.FREEC.bam_ratio.sorted.resegmented.CNVs.txt ]]
    then
	Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/assess_CNV_significance_for_FREEC.R \
	    --cnv $parent1_based_prefix.FREEC.bam_ratio.sorted.resegmented.CNVs.txt  \
	    --ratio $parent1_based_prefix.FREEC.bam_ratio.sorted.adjusted.txt \
	    --genome_fai $parent1_tag.genome.raw.fa.fai \
	    --output $parent1_based_prefix.CNV_significance_test.txt
	cat $parent1_based_prefix.CNV_significance_test.txt | tail -n +2 | awk '{ if ($6 < 0.05) print $0 }' | awk '{print $1, $2-1, $3}' OFS='\t' > $parent1_based_prefix.significant_CNV.bed
    else
        echo -e "chr\tstart\tend\tcopy_number\tstatus\tMWU_test_p_value\tKS_test_p_value" > $parent1_based_prefix.CNV_significance_test.txt
	echo -e "\t\t" > $parent1_based_prefix.significant_CNV.bed
	echo 
    fi
else
    echo "Exception encountered for FREEC! Exit! ..." > $parent1_based_prefix.realn.bam.no_FREEC.txt
    exit
fi

if [[ $debug == "no" ]]
then
    rm for_CNV.bam
    rm FREEC.bam
    rm FREEC.bam.header.old.sam
    rm FREEC.bam.header.new.sam
    rm GC_profile.${window_size}bp.cnp
    rm $parent2_tag.R1.raw.fq.gz
    rm $parent2_tag.R2.raw.fq.gz
fi

# SNP and INDEL calling

python3 $freebayes_dir/../scripts/fasta_generate_regions.py $parent1_tag.genome.raw.fa.fai 100000 > $parent1_tag.genome.region.txt
cat $parent1_tag.genome.region.txt | $parallel_dir/parallel -k -j $threads $freebayes_dir/freebayes -f $parent1_tag.genome.raw.fa  -p $ploidy $parent1_based_prefix.realn.bam --region {} | python3 $vcflib_dir/../scripts/vcffirstheader | $vcflib_dir/vcfstreamsort -w 1000 | $vcflib_dir/vcfuniq > $parent1_based_prefix.caller.raw.vcf 

# cat $parent1_tag.genome.raw.fa.fai |cut -f 1 | $parallel_dir/parallel -I% --max-args 1 -k -j $threads $java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=1 -jar $gatk4_dir/gatk4.jar HaplotypeCaller -R $parent1_tag.genome.raw.fa -I $parent1_based_prefix.realn.bam  -ploidy $ploidy -L % -O $parent1_based_prefix.caller.raw.by_chr.%.vcf

# cat $parent1_based_prefix.caller.raw.by_chr.*.vcf | python3 $vcflib_dir/../scripts/vcffirstheader | $vcflib_dir/vcfstreamsort -w 1000 | $vcflib_dir/vcfuniq > $parent1_based_prefix.caller.raw.unsort.vcf

# java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar SortVcf  -SEQUENCE_DICTIONARY $parent1_tag.genome.raw.dict -I $parent1_based_prefix.caller.raw.unsort.vcf -O $parent1_based_prefix.caller.raw.vcf

# if [[ $debug == "no" ]]
# then
#     rm $parent1_based_prefix.caller.raw.by_chr.*.vcf 
#     rm $parent1_based_prefix.caller.raw.unsort.vcf
#     rm $parent1_based_prefix.caller.raw.*.idx
# fi

$vt_dir/vt decompose_blocksub $parent1_based_prefix.caller.raw.vcf -a -o $parent1_based_prefix.caller.decompose.vcf
$vt_dir/vt normalize $parent1_based_prefix.caller.decompose.vcf -r $parent1_tag.genome.raw.fa -f "~VARIANT_CONTAINS_N"| $vt_dir/vt uniq - -o  $parent1_based_prefix.caller.normalize.vcf
$vt_dir/vt annotate_variants $parent1_based_prefix.caller.normalize.vcf -r $parent1_tag.genome.raw.fa -o $parent1_based_prefix.caller.annotate.vcf

cat $parent1_based_prefix.caller.annotate.vcf | egrep "^#" > $parent1_based_prefix.caller.annotate.vcf.header
$bedtools_dir/bedtools subtract -a $parent1_based_prefix.caller.annotate.vcf -b ./../$parent1_hardmask_bed > $parent1_based_prefix.caller.annotate_hardmask.vcf.content
cat $parent1_based_prefix.caller.annotate.vcf.header $parent1_based_prefix.caller.annotate_hardmask.vcf.content > $parent1_based_prefix.caller.annotate_hardmask.vcf

parent1_based_CNV_bed_line_count=$(cat $parent1_based_prefix.significant_CNV.bed |sed '/^\s*$/d' | wc -l)
if [[ "$apply_cnv_filter" == "yes" && "$parent1_based_CNV_bed_line_count" > 0 ]]
then
    $bedtools_dir/bedtools subtract -a $parent1_based_prefix.caller.annotate_hardmask.vcf -b $parent1_based_prefix.significant_CNV.bed > $parent1_based_prefix.caller.annotate_hardmask_CNVmask.vcf.content
    cat $parent1_based_prefix.caller.annotate.vcf.header $parent1_based_prefix.caller.annotate_hardmask_CNVmask.vcf.content > $parent1_based_prefix.caller.annotate_hardmask_CNVmask.vcf
    rm $parent1_based_prefix.caller.annotate.vcf.header
    rm $parent1_based_prefix.caller.annotate_hardmask.vcf.content
    rm $parent1_based_prefix.caller.annotate_hardmask_CNVmask.vcf.content
else
    cp $parent1_based_prefix.caller.annotate_hardmask.vcf $parent1_based_prefix.caller.annotate_hardmask_CNVmask.vcf
    rm $parent1_based_prefix.caller.annotate.vcf.header
    rm $parent1_based_prefix.caller.annotate_hardmask.vcf.content
fi

perl $RECOMBINEX_HOME/scripts/filter_vcf_by_window.pl \
    -i $parent1_based_prefix.caller.annotate_hardmask_CNVmask.vcf \
    -w $cluster_window_size \
    -o $parent1_based_prefix.caller.annotate_hardmask_CNVmask.thin.vcf

$vt_dir/vt view $parent1_based_prefix.caller.annotate_hardmask_CNVmask.thin.vcf -f "VTYPE==SNP" -o $parent1_based_prefix.caller.annotate_hardmask_CNVmask.thin.SNP.vcf
$vt_dir/vt view $parent1_based_prefix.caller.annotate_hardmask_CNVmask.thin.vcf -f "VTYPE==INDEL" -o $parent1_based_prefix.caller.annotate_hardmask_CNVmask.thin.INDEL.vcf

$vcflib_dir/vcffilter -f "QUAL > $variant_calling_quality_cutoff & QUAL / AO > 1 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" $parent1_based_prefix.caller.annotate_hardmask_CNVmask.thin.SNP.vcf |egrep "^#|AC=1" \
    > $parent1_based_prefix.read_mapping.SNP.filter.vcf
$vcflib_dir/vcffilter -f "QUAL > $variant_calling_quality_cutoff & QUAL / AO > 1 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" $parent1_based_prefix.caller.annotate_hardmask_CNVmask.thin.INDEL.vcf |egrep "^#|AC=1" \
    > $parent1_based_prefix.read_mapping.INDEL.filter.vcf

rm -r tmp
# compress vcf files
gzip *.vcf 

# copy the filtered vcf files back to the task-specific directory
cp $parent1_based_prefix.read_mapping.SNP.filter.vcf.gz ./../
# cp $parent1_based_prefix.read_mapping.INDEL.filter.vcf.gz ./../

if [[ $debug == "no" ]]
then
    rm $parent1_based_prefix.caller.raw.vcf.gz
    rm $parent1_based_prefix.caller.decompose.vcf.gz
    rm $parent1_based_prefix.caller.normalize.vcf.gz
    # rm $parent1_based_prefix.caller.annotate.vcf.gz
    # rm $parent1_based_prefix.caller.annotate_hardmask.vcf.gz
    # rm $parent1_based_prefix.caller.annotate_hardmask_CNVmask.vcf.gz
    # rm $parent1_based_prefix.caller.annotate_hardmask_CNVmask.thin.vcf.gz
    # rm $parent1_based_prefix.caller.annotate_hardmask_CNVmask.thin.SNP.vcf.gz
    # rm $parent1_based_prefix.caller.annotate_hardmask_CNVmask.thin.INDEL.vcf.gz
fi

cd ..

#######

cd $parent2_based_output_dir
ln -s ./../$parent2_raw_assembly $parent2_tag.genome.raw.fa
ln -s ./../$parent1_reads_R1 $parent1_tag.R1.raw.fq.gz
ln -s ./../$parent1_reads_R2 $parent1_tag.R2.raw.fq.gz

cp $adapter adapter.fa

mkdir tmp

$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $trimmomatic_dir/trimmomatic.jar PE -threads $threads -phred33 $parent1_tag.R1.raw.fq.gz $parent1_tag.R2.raw.fq.gz  $parent1_tag.R1.trimmed.PE.fq.gz $parent1_tag.R1.trimmed.SE.fq.gz $parent1_tag.R2.trimmed.PE.fq.gz $parent1_tag.R2.trimmed.SE.fq.gz ILLUMINACLIP:adapter.fa:2:30:10  SLIDINGWINDOW:5:20 MINLEN:36

if [[ $debug == "no" ]]
then
    rm $parent1_tag.R1.trimmed.SE.fq.gz
    rm $parent1_tag.R2.trimmed.SE.fq.gz
    rm adapter.fa
fi

$bwa_dir/bwa index $parent2_tag.genome.raw.fa
$bwa_dir/bwa mem -t $threads -M $parent2_tag.genome.raw.fa $parent1_tag.R1.trimmed.PE.fq.gz $parent1_tag.R2.trimmed.PE.fq.gz | $samtools_dir/samtools view -bS -q $mapping_quality_cutoff -F 3340 -f 2 - >$parent2_based_prefix.bam

if [[ $debug == "no" ]]
then
    rm $parent2_tag.genome.raw.fa.bwt
    rm $parent2_tag.genome.raw.fa.pac
    rm $parent2_tag.genome.raw.fa.ann
    rm $parent2_tag.genome.raw.fa.amb
    rm $parent2_tag.genome.raw.fa.sa
    rm $parent1_tag.R1.trimmed.PE.fq.gz
    rm $parent1_tag.R2.trimmed.PE.fq.gz
fi

# index reference sequence
$samtools_dir/samtools faidx $parent2_tag.genome.raw.fa
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar CreateSequenceDictionary \
    -REFERENCE $parent2_tag.genome.raw.fa \
    -OUTPUT $parent2_tag.genome.raw.dict

# sort bam file by picard-tools SortSam
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar SortSam \
    -INPUT $parent2_based_prefix.bam \
    -OUTPUT $parent2_based_prefix.sort.bam \
    -SORT_ORDER coordinate

if [[ $debug == "no" ]]
then
    rm $parent2_based_prefix.bam
fi

# fixmate
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar FixMateInformation \
    -INPUT $parent2_based_prefix.sort.bam \
    -OUTPUT $parent2_based_prefix.fixmate.bam

if [[ $debug == "no" ]]
then
    rm $parent2_based_prefix.sort.bam
fi

# add or replace read groups and sort
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar AddOrReplaceReadGroups \
    -INPUT $parent2_based_prefix.fixmate.bam \
    -OUTPUT $parent2_based_prefix.rdgrp.bam \
    -SORT_ORDER coordinate \
    -RGID "$parent2_based_prefix" \
    -RGLB "$parent2_based_prefix" \
    -RGPL "Illumina" \
    -RGPU "$parent2_based_prefix" \
    -RGSM "$parent2_based_prefix" \
    -RGCN "RGCN"

if [[ $debug == "no" ]]
then
    rm $parent2_based_prefix.fixmate.bam
fi

# remove duplicates
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar MarkDuplicates \
    -INPUT $parent2_based_prefix.rdgrp.bam \
    -REMOVE_DUPLICATES true  \
    -METRICS_FILE $parent2_based_prefix.dedup.matrics \
    -OUTPUT $parent2_based_prefix.dedup.bam 

# index the dedup.bam file
$samtools_dir/samtools index $parent2_based_prefix.dedup.bam

if [[ $debug == "no" ]]
then
    rm $parent2_based_prefix.rdgrp.bam
fi

# GATK local realign
# find realigner targets
$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $gatk3_dir/gatk3.jar \
    -nt $threads \
    -R $parent2_tag.genome.raw.fa \
    -T RealignerTargetCreator \
    -I $parent2_based_prefix.dedup.bam \
    -o $parent2_based_prefix.realn.intervals
# run realigner
$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $gatk3_dir/gatk3.jar \
    -R $parent2_tag.genome.raw.fa \
    -T IndelRealigner \
    -I $parent2_based_prefix.dedup.bam \
    -targetIntervals $parent2_based_prefix.realn.intervals  \
    -o $parent2_based_prefix.realn.bam

if [[ $debug == "no" ]]
then
    rm $parent2_based_prefix.dedup.bam
    rm $parent2_based_prefix.dedup.bam.bai
    rm $parent2_based_prefix.dedup.matrics
    rm $parent2_based_prefix.realn.intervals
fi

# generate samtools mpileup
$samtools_dir/samtools mpileup -C 0 -q $mapping_quality_cutoff -f $parent2_tag.genome.raw.fa $parent2_based_prefix.realn.bam |gzip -c >$parent2_based_prefix.mpileup.gz

# calculate per-base depth
$samtools_dir/samtools depth -aa $parent2_based_prefix.realn.bam |gzip -c >$parent2_based_prefix.depth.txt.gz

# compute basic alignment statistics by samtools
$samtools_dir/samtools flagstat $parent2_based_prefix.realn.bam >$parent2_based_prefix.samstat

# compute insert size statistics
$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar CollectInsertSizeMetrics \
    -I $parent2_based_prefix.realn.bam \
    -O $parent2_based_prefix.insert_size_metrics.txt \
    -H $parent2_based_prefix.insert_size_histogram.pdf \
    -M 0.5

# calculate read mapping coverage statistics
perl $RECOMBINEX_HOME/scripts/summarize_mapping_coverage.pl \
    -r $parent2_tag.genome.raw.fa \
    -s $parent2_based_prefix.samstat \
    -d $parent2_based_prefix.depth.txt.gz \
    -c 5 \
    -t $parent2_tag \
    -o $parent2_based_prefix.coverage_summary.txt


# scan for CNV by FREEC

#parent1_raw_read_length=$(gunzip -c $parent1_tag.R1.raw.fq.gz |awk 'NR%4==2{print length($0)}' | head -1 || true)
parent1_raw_read_length=100
parent2_min_expected_gc=$(cat ./../$parent_genome_preprocessing_dir/$parent2_tag.FREEC.GC_range.txt|egrep -v "#"|cut -f 1)
parent2_max_expected_gc=$(cat ./../$parent_genome_preprocessing_dir/$parent2_tag.FREEC.GC_range.txt|egrep -v "#"|cut -f 2)
echo ""
echo "window_size=$window_size, step_size=$step_size, read_length=$parent1_raw_read_length"
echo "parent2_min_expected_gc=$parent2_min_expected_gc, parent2_max_expected_gc=$parent2_max_expected_gc"
echo ""


if [ -z "$excluded_chr_list_for_cnv_profiling" ]
then
perl $RECOMBINEX_HOME/scripts/run_FREEC_wrapper_lite_for_RecombineX.pl \
    -r $parent2_tag \
    -bam $parent2_based_prefix.realn.bam \
    -prefix $parent2_based_prefix \
    -ploidy $ploidy \
    -bedtools $bedtools_dir/bedtools \
    -samtools $samtools_dir/samtools \
    -freec $freec_dir/freec \
    -window $window_size \
    -step $step_size \
    -read_length_for_mappability $parent1_raw_read_length  \
    -min_mappability $min_mappability \
    -min_expected_gc $parent2_min_expected_gc \
    -max_expected_gc $parent2_max_expected_gc \
    -refseq_genome_preprocessing_dir "./../$parent_genome_preprocessing_dir" \
    -mates_orientation 0 \
    -threads $threads 
else
    echo "excluded chr list: $excluded_chr_list_for_cnv_profiling"
    perl $RECOMBINEX_HOME/scripts/run_FREEC_wrapper_lite_for_RecombineX.pl \
	-r $parent2_tag \
	-bam $parent2_based_prefix.realn.bam \
	-prefix $parent2_based_prefix \
	-ploidy $ploidy \
	-bedtools $bedtools_dir/bedtools \
	-samtools $samtools_dir/samtools \
	-freec $freec_dir/freec \
	-window $window_size \
	-step $step_size \
	-read_length_for_mappability $parent1_raw_read_length  \
	-min_mappability $min_mappability \
	-min_expected_gc $parent2_min_expected_gc \
	-max_expected_gc $parent2_max_expected_gc \
	-refseq_genome_preprocessing_dir "./../$parent_genome_preprocessing_dir" \
	-mates_orientation 0 \
	-threads $threads \
	-excluded_chr_list ./../$excluded_chr_list_for_cnv_profiling
fi
Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/CNV_segmentation_by_DNAcopy.R \
    --input $parent2_based_prefix.FREEC.bam_ratio.txt \
    --prefix $parent2_based_prefix \
    --window $window_size \
    --ploidy $ploidy \
    --genome_fai $parent2_tag.genome.raw.fa.fai

perl $RECOMBINEX_HOME/scripts/adjust_FREEC_copynumber_by_DNAcopy_copynumber.pl \
    -i $parent2_based_prefix.FREEC.bam_ratio.sorted.txt \
    -a $parent2_based_prefix.FREEC.bam_ratio.sorted.resegmented.lite.txt \
    -o $parent2_based_prefix.FREEC.bam_ratio.sorted.adjusted.txt

if [[ -s $parent2_based_prefix.FREEC.bam_ratio.sorted.adjusted.txt ]]
then
    Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_CNV_for_FREEC.R \
	--ploidy $ploidy \
	--genome_fai $parent2_tag.genome.raw.fa.fai \
	--input $parent2_based_prefix.FREEC.bam_ratio.sorted.adjusted.txt \
	--output $parent2_based_prefix.CNV_plot.pdf
    rm Rplots.pdf
    if [[ -s $parent2_based_prefix.FREEC.bam_ratio.sorted.resegmented.CNVs.txt ]]
    then
	Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/assess_CNV_significance_for_FREEC.R \
	    --cnv $parent2_based_prefix.FREEC.bam_ratio.sorted.resegmented.CNVs.txt \
	    --ratio $parent2_based_prefix.FREEC.bam_ratio.sorted.adjusted.txt \
	    --genome_fai $parent2_tag.genome.raw.fa.fai \
	    --output $parent2_based_prefix.CNV_significance_test.txt
	cat $parent2_based_prefix.CNV_significance_test.txt | tail -n +2 | awk '{ if ($6 < 0.05) print $0 }' | awk '{print $1, $2-1, $3}' OFS='\t' > $parent2_based_prefix.significant_CNV.bed
    else
        echo -e "chr\tstart\tend\tcopy_number\tstatus\tMWU_test_p_value\tKS_test_p_value" > $parent2_based_prefix.CNV_significance_test.txt
	echo -e "\t\t" > $parent2_based_prefix.significant_CNV.bed
    fi
else
    echo "Exception encountered for FREEC! Exit! ..." > $parent2_based_prefix.realn.bam.no_FREEC.txt
    exit
fi


if [[ $debug == "no" ]]
then
    rm for_CNV.bam
    rm FREEC.bam
    rm FREEC.bam.header.old.sam
    rm FREEC.bam.header.new.sam
    rm GC_profile.${window_size}bp.cnp
    rm $parent1_tag.R1.raw.fq.gz
    rm $parent1_tag.R2.raw.fq.gz
fi

# SNP and INDEL calling

python3 $freebayes_dir/../scripts/fasta_generate_regions.py $parent2_tag.genome.raw.fa.fai 100000 > $parent2_tag.genome.region.txt
cat $parent2_tag.genome.region.txt | $parallel_dir/parallel -k -j $threads $freebayes_dir/freebayes -f $parent2_tag.genome.raw.fa  -p $ploidy $parent2_based_prefix.realn.bam --region {} | python3 $vcflib_dir/../scripts/vcffirstheader | $vcflib_dir/vcfstreamsort -w 1000 | $vcflib_dir/vcfuniq > $parent2_based_prefix.caller.raw.vcf

# cat $parent2_tag.genome.raw.fa.fai |cut -f 1 | $parallel_dir/parallel -I% --max-args 1 -k -j $threads $java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=1 -jar $gatk4_dir/gatk4.jar HaplotypeCaller -R $parent2_tag.genome.raw.fa -I $parent2_based_prefix.realn.bam  -ploidy $ploidy -L % -O $parent2_based_prefix.caller.raw.by_chr.%.vcf
# cat $parent2_based_prefix.caller.raw.by_chr.*.vcf | python3 $vcflib_dir/../scripts/vcffirstheader | $vcflib_dir/vcfstreamsort -w 1000 | $vcflib_dir/vcfuniq > $parent2_based_prefix.caller.raw.unsort.vcf
# java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar SortVcf  -SEQUENCE_DICTIONARY $parent2_tag.genome.raw.dict -I $parent2_based_prefix.caller.raw.unsort.vcf -O $parent2_based_prefix.caller.raw.vcf

# if [[ $debug == "no" ]]
# then
#     rm $parent2_based_prefix.caller.raw.by_chr.*.vcf 
#     rm $parent2_based_prefix.caller.raw.unsort.vcf
#     rm $parent2_based_prefix.caller.raw.*.idx

# fi

$vt_dir/vt decompose_blocksub $parent2_based_prefix.caller.raw.vcf -a -o $parent2_based_prefix.caller.decompose.vcf

$vt_dir/vt normalize $parent2_based_prefix.caller.decompose.vcf -r $parent2_tag.genome.raw.fa -f "~VARIANT_CONTAINS_N"| $vt_dir/vt uniq - -o  $parent2_based_prefix.caller.normalize.vcf
$vt_dir/vt annotate_variants $parent2_based_prefix.caller.normalize.vcf -r $parent2_tag.genome.raw.fa -o $parent2_based_prefix.caller.annotate.vcf

cat $parent2_based_prefix.caller.annotate.vcf | egrep "^#" > $parent2_based_prefix.caller.annotate.vcf.header
$bedtools_dir/bedtools subtract -a $parent2_based_prefix.caller.annotate.vcf -b ./../$parent2_hardmask_bed > $parent2_based_prefix.caller.annotate_hardmask.vcf.content
cat $parent2_based_prefix.caller.annotate.vcf.header $parent2_based_prefix.caller.annotate_hardmask.vcf.content > $parent2_based_prefix.caller.annotate_hardmask.vcf

parent2_based_CNV_bed_line_count=$(cat $parent2_based_prefix.significant_CNV.bed |sed '/^\s*$/d' | wc -l)
if [[ "$apply_cnv_filter" == "yes" && "$parent2_based_CNV_bed_line_count" > 0 ]]
then
    $bedtools_dir/bedtools subtract -a $parent2_based_prefix.caller.annotate_hardmask.vcf -b $parent2_based_prefix.significant_CNV.bed > $parent2_based_prefix.caller.annotate_hardmask_CNVmask.vcf.content
    cat $parent2_based_prefix.caller.annotate.vcf.header $parent2_based_prefix.caller.annotate_hardmask_CNVmask.vcf.content > $parent2_based_prefix.caller.annotate_hardmask_CNVmask.vcf
    rm $parent2_based_prefix.caller.annotate.vcf.header
    rm $parent2_based_prefix.caller.annotate_hardmask.vcf.content
    rm $parent2_based_prefix.caller.annotate_hardmask_CNVmask.vcf.content
else
    cp $parent2_based_prefix.caller.annotate_hardmask.vcf $parent2_based_prefix.caller.annotate_hardmask_CNVmask.vcf
    rm $parent2_based_prefix.caller.annotate.vcf.header
    rm $parent2_based_prefix.caller.annotate_hardmask.vcf.content
fi

perl $RECOMBINEX_HOME/scripts/filter_vcf_by_window.pl \
    -i $parent2_based_prefix.caller.annotate_hardmask_CNVmask.vcf \
    -w $cluster_window_size \
    -o $parent2_based_prefix.caller.annotate_hardmask_CNVmask.thin.vcf

$vt_dir/vt view $parent2_based_prefix.caller.annotate_hardmask_CNVmask.thin.vcf -f "VTYPE==SNP" -o $parent2_based_prefix.caller.annotate_hardmask_CNVmask.thin.SNP.vcf
$vt_dir/vt view $parent2_based_prefix.caller.annotate_hardmask_CNVmask.thin.vcf -f "VTYPE==INDEL" -o $parent2_based_prefix.caller.annotate_hardmask_CNVmask.thin.INDEL.vcf

$vcflib_dir/vcffilter -f "QUAL > $variant_calling_quality_cutoff & QUAL / AO > 1 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" $parent2_based_prefix.caller.annotate_hardmask_CNVmask.thin.SNP.vcf |egrep "^#|AC=1" \
    > $parent2_based_prefix.read_mapping.SNP.filter.vcf
$vcflib_dir/vcffilter -f "QUAL > $variant_calling_quality_cutoff & QUAL / AO > 1 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" $parent2_based_prefix.caller.annotate_hardmask_CNVmask.thin.INDEL.vcf |egrep "^#|AC=1" \
    > $parent2_based_prefix.read_mapping.INDEL.filter.vcf

rm -r tmp

# compress the vcf files
gzip *.vcf

# copy the filtered vcf files back to the task-specific directory
cp $parent2_based_prefix.read_mapping.SNP.filter.vcf.gz ./../
# cp $parent2_based_prefix.read_mapping.INDEL.filter.vcf.gz ./../

if [[ $debug == "no" ]]
then
    rm $parent2_based_prefix.caller.raw.vcf.gz
    rm $parent2_based_prefix.caller.decompose.vcf.gz
    rm $parent2_based_prefix.caller.normalize.vcf.gz
    # rm $parent2_based_prefix.caller.annotate.vcf.gz
    # rm $parent2_based_prefix.caller.annotate_hardmask.vcf.gz
    # rm $parent2_based_prefix.caller.annotate_hardmask_CNVmask.vcf.gz
    # rm $parent2_based_prefix.caller.annotate_hardmask_CNVmask.thin.vcf.gz
    # rm $parent2_based_prefix.caller.annotate_hardmask_CNVmask.thin.SNP.vcf.gz
    # rm $parent2_based_prefix.caller.annotate_hardmask_CNVmask.thin.INDEL.vcf.gz
fi

cd ..


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


