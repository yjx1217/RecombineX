#!/bin/bash
set -e -o pipefail

##########################################
# load environment variables for RecombineX
source ./../../env.sh

###########################################
# set project-specific variables
parent1_tag="S288C" # The relabeling tag of genome 1 used in 11.Parent_Genome_Preprocessing. Default = "S288C".
parent2_tag="SK1" # The relabeling tag of genome 2 used in 11.Parent_Genome_Preprocessing. Default = "SK1".
genome_aln_mode="by_chromosome" # The mode for running whole genome alignment, can be "by_chromosome" (default) or "by_genome" (when interchromosomal rearrangements are involved, e.g. translocations). Default = "by_chromosome".
chr_list="$RECOMBINEX_HOME/data/Saccharomyces_cerevisiae.chr_list.txt" # The chromosome list for whole-genome alignment by chromosome, only needed when genome_aln_mode="by_chromosome". Set this option to "" when genome_aln_mode="by_genome". Default = "$RECOMBINEX_HOME/data/Saccharomyces_cerevisiae.chr_list.txt".
use_centromere_annotation="yes" # whether to use the centromere annotation information. Please note that enabling this option requires that you have the parent1_tag.centromere.relabel.gff and parent2_tag.centromere.relabel.gff files ready in the "./../11.Parent_Genome_Preprocessing" directory. Set this option to "no" when running RecombineX for the mitochondrial genome. Default = "yes".
threads=4 # The number of threads to use. Default = "4".
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
###########################################



# process the pipeline

#############################################
# Normally no need to change the following parameters
mincluster=65 # The minimum length of a cluster of matches: the "-c|mincluster" option for mummer's nucmer program. Default = "65".
minmatch=20 # The the minimum length of a single match: the "-l|minmatch" option for mummer's nucmer program. Default = "20".
maxgap=90 # The the maximum gap between two adjacent matches in a cluster: the "-g|maxgap" option for mummer's nucmer progam. Default = "90".
cluster_window_size=10 # Adjacent variants within specified window (unit: bp) will be all filtered out if any of them is INDEL. Default = "10".
############################################

test_file_existence () {
    filename=$1
    if [[ ! -f $filename ]]
    then
        echo "the file $filename does not exists! process terminated!"
        exit
    fi
}

# check option settings

parent1_raw_assembly="./../11.Parent_Genome_Preprocessing/$parent1_tag.genome.raw.relabel.fa" # relabeled genome 1 raw assembly file generated in 01.Parent_Genome_Preprocessing
parent2_raw_assembly="./../11.Parent_Genome_Preprocessing/$parent2_tag.genome.raw.relabel.fa" # relabeled genome 2 raw assembly file generated in 01.Parent_Genome_Preprocessing

parent1_hardmask_bed="./../11.Parent_Genome_Preprocessing/$parent1_tag.genome.hardmask.relabel.masking_details.bed" # the masking details bed file for the relabeled and hardmasked genome 1 assembly file generated in 01.Parent_Genome_Preprocessing
parent2_hardmask_bed="./../11.Parent_Genome_Preprocessing/$parent2_tag.genome.hardmask.relabel.masking_details.bed" # the masking details bed file for the relabeled and hardmasked genome 2 assembly file generated in 01.Parent_Genome_Preprocessing

echo ""
echo "check the existence of parent1_raw_assembly .."
test_file_existence $parent1_raw_assembly
echo ""
echo "check the existence of parent2_raw_assembly .."
test_file_existence $parent2_raw_assembly
echo ""
echo "test the existence of parent1_hardmask_bed"
test_file_existence $parent1_hardmask_bed
echo ""
echo "test the existence of parent2_hardmask_bed"
test_file_existence $parent2_hardmask_bed


if [[ $use_centromere_annotation == "yes" ]]
then
    parent1_centromere_gff="./../11.Parent_Genome_Preprocessing/$parent1_tag.centromere.relabel.gff" # relabeled centromere gff file generated in 01.Parent_Genome_Preprocessing; only needed if plot_centromere="yes"
    parent2_centromere_gff="./../11.Parent_Genome_Preprocessing/$parent2_tag.centromere.relabel.gff" # relabeled centromere gff file generated in 01.Parent_Genome_Preprocessing; only needed if plot_centromere="yes"
echo ""
echo "check the existence of parent1_centromere .."
test_file_existence $parent1_centromere_gff
echo ""
echo "check the existence of parent2_centromere .."
test_file_existence $parent2_centromere_gff
fi

parent1_based_output_dir="${parent1_tag}-${parent2_tag}_${parent1_tag}_based"
parent2_based_output_dir="${parent1_tag}-${parent2_tag}_${parent2_tag}_based"
mkdir -p ./${parent1_based_output_dir}/${parent1_tag}
mkdir -p ./${parent1_based_output_dir}/${parent2_tag}
mkdir -p ./${parent2_based_output_dir}/${parent1_tag}
mkdir -p ./${parent2_based_output_dir}/${parent2_tag}

cp $parent1_raw_assembly $parent1_tag.genome.fa
cp $parent2_raw_assembly $parent2_tag.genome.fa

$samtools_dir/samtools faidx $parent1_tag.genome.fa
$samtools_dir/samtools faidx $parent2_tag.genome.fa

parent1_fai="$parent1_tag.genome.fa.fai" 
parent2_fai="$parent2_tag.genome.fa.fai"

parent1_based_prefix="${parent1_tag}-${parent2_tag}.${parent1_tag}"
parent2_based_prefix="${parent1_tag}-${parent2_tag}.${parent2_tag}"

parent1_based_vcf_file_list=""
parent2_based_vcf_file_list=""

if [[ $genome_aln_mode == "by_chromosome" ]]
then
    echo ""
    echo "check the existence of chr_list .."
    test_file_existence $chr_list
    sed -i '/^$/d' $chr_list
    $ucsc_dir/faSplit byname $parent1_tag.genome.fa ./${parent1_based_output_dir}/${parent1_tag}/
    $ucsc_dir/faSplit byname $parent2_tag.genome.fa ./${parent1_based_output_dir}/${parent2_tag}/
    $ucsc_dir/faSplit byname $parent1_tag.genome.fa ./${parent2_based_output_dir}/${parent1_tag}/
    $ucsc_dir/faSplit byname $parent2_tag.genome.fa ./${parent2_based_output_dir}/${parent2_tag}/
    
    while read chr
    do
	parent1_based_prefix="${parent1_tag}-${parent2_tag}.${parent1_tag}.${chr}"
	parent2_based_prefix="${parent1_tag}-${parent2_tag}.${parent2_tag}.${chr}"
	$mummer3_dir/nucmer -g $maxgap -l $minmatch -c $mincluster -p $parent1_based_prefix ./${parent1_based_output_dir}/${parent1_tag}/${parent1_tag}_${chr}.fa ./${parent1_based_output_dir}/${parent2_tag}/${parent2_tag}_${chr}.fa
	$mummer3_dir/nucmer -g $maxgap -l $minmatch -c $mincluster -p $parent2_based_prefix ./${parent2_based_output_dir}/${parent2_tag}/${parent2_tag}_${chr}.fa ./${parent2_based_output_dir}/${parent1_tag}/${parent1_tag}_${chr}.fa

	$mummer3_dir/delta-filter -1 $parent1_based_prefix.delta > $parent1_based_prefix.1delta
	$mummer3_dir/delta-filter -1 $parent2_based_prefix.delta > $parent2_based_prefix.1delta                                                                                                                                                
	#unique_length=2000 # The minimal size of unique anchors. Default = "2000". This is the "--unique-length" option for Assemblytics; In general, a minimal value of "1000" is recommended. Here we want to be more conservative in marker identification, so a value of "2000" is chosen. If you are working with mammalian-size genomes, a value of "10000" is recommended. 
	# $assemblytics_dir/Assemblytics_uniq_anchor.py --delta $parent1_based_prefix.delta --unique-length $unique_length --keep-small-uniques --out $parent1_based_prefix.1delta
	# $assemblytics_dir/Assemblytics_uniq_anchor.py --delta $parent2_based_prefix.delta --unique-length $unique_length --keep-small-uniques --out $parent2_based_prefix.1delta

	# gunzip -c  $parent1_based_prefix.1delta.Assemblytics.unique_length_filtered_l${unique_length}.delta.gz >$parent1_based_prefix.1delta
	# gunzip -c  $parent2_based_prefix.1delta.Assemblytics.unique_length_filtered_l${unique_length}.delta.gz >$parent2_based_prefix.1delta

	$mummer3_dir/show-coords -T -r -c -l -d   $parent1_based_prefix.1delta > $parent1_based_prefix.1coords
	$mummer3_dir/show-coords -T -r -c -l -d   $parent2_based_prefix.1delta > $parent2_based_prefix.1coords

	$mummer3_dir/show-snps -Clr -T $parent1_based_prefix.1delta > $parent1_based_prefix.snps
	$mummer3_dir/show-snps -Clr -T $parent2_based_prefix.1delta > $parent2_based_prefix.snps

	perl $RECOMBINEX_HOME/scripts/mummer2vcf.pl -r ${parent1_tag}.genome.fa -i $parent1_based_prefix.snps -t BOTH -p $parent1_based_prefix
	perl $RECOMBINEX_HOME/scripts/mummer2vcf.pl -r ${parent2_tag}.genome.fa -i $parent2_based_prefix.snps -t BOTH -p $parent2_based_prefix

	awk '{printf("##contig=<ID=%s,length=%d>\n",$1,$2);}' $parent1_tag.genome.fa.fai > $parent1_based_prefix.vcf_header.sequences.txt
	awk '{printf("##contig=<ID=%s,length=%d>\n",$1,$2);}' $parent2_tag.genome.fa.fai > $parent2_based_prefix.vcf_header.sequences.txt
	sed -i -e "/##reference/r $parent1_based_prefix.vcf_header.sequences.txt" $parent1_based_prefix.mummer2vcf.BOTH.vcf
	sed -i -e "/##reference/r $parent2_based_prefix.vcf_header.sequences.txt" $parent2_based_prefix.mummer2vcf.BOTH.vcf
	mv  $parent1_based_prefix.vcf_header.sequences.txt $parent1_based_output_dir
	mv  $parent2_based_prefix.vcf_header.sequences.txt $parent2_based_output_dir

	$vt_dir/vt sort $parent1_based_prefix.mummer2vcf.BOTH.vcf -o $parent1_based_prefix.mummer2vcf.BOTH.sort.vcf
	$vt_dir/vt sort $parent2_based_prefix.mummer2vcf.BOTH.vcf -o $parent2_based_prefix.mummer2vcf.BOTH.sort.vcf

	parent1_based_vcf_file_list+="$parent1_based_prefix.mummer2vcf.BOTH.sort.vcf "
	parent2_based_vcf_file_list+="$parent2_based_prefix.mummer2vcf.BOTH.sort.vcf "

	mv $parent1_based_prefix.delta $parent1_based_output_dir
	mv $parent1_based_prefix.1delta $parent1_based_output_dir
	#mv $parent1_based_prefix.1delta.* $parent1_based_output_dir
	mv $parent1_based_prefix.1coords $parent1_based_output_dir
	mv $parent1_based_prefix.snps $parent1_based_output_dir
	mv $parent1_based_prefix.mummer2vcf.BOTH.vcf $parent1_based_output_dir

	mv $parent2_based_prefix.delta $parent2_based_output_dir
	mv $parent2_based_prefix.1delta $parent2_based_output_dir
	# mv $parent2_based_prefix.1delta.* $parent2_based_output_dir
	mv $parent2_based_prefix.1coords $parent2_based_output_dir
	mv $parent2_based_prefix.snps $parent2_based_output_dir
	mv $parent2_based_prefix.mummer2vcf.BOTH.vcf $parent2_based_output_dir
    done < $chr_list
    cat $parent1_based_prefix.mummer2vcf.BOTH.sort.vcf | egrep "^#"  > ${parent1_tag}-${parent2_tag}.${parent1_tag}.vcf_header.txt
    cat $parent2_based_prefix.mummer2vcf.BOTH.sort.vcf | egrep "^#"  > ${parent1_tag}-${parent2_tag}.${parent2_tag}.vcf_header.txt
else
    $mummer3_dir/nucmer -g $maxgap -l $minmatch -c $mincluster -p $parent1_based_prefix $parent1_tag.genome.fa $parent2_tag.genome.fa
    $mummer3_dir/nucmer -g $maxgap -l $minmatch -c $mincluster -p $parent2_based_prefix $parent2_tag.genome.fa $parent1_tag.genome.fa

    $mummer3_dir/delta-filter -1 $parent1_based_prefix.delta > $parent1_based_prefix.1delta
    $mummer3_dir/delta-filter -1 $parent2_based_prefix.delta > $parent2_based_prefix.1delta                                                                                                                                                
    #unique_length=2000 # The minimal size of unique anchors. Default = "2000". This is the "--unique-length" option for Assemblytics; In general, a minimal value of "1000" is recommended. Here we want to be more conservative in marker identification, so a value of "2000" is chosen. If you are working with mammalian-size genomes, a value of "10000" is recommended. 
    #$assemblytics_dir/Assemblytics_uniq_anchor.py --delta $parent1_based_prefix.delta --unique-length $unique_length --keep-small-uniques --out $parent1_based_prefix.1delta
    #$assemblytics_dir/Assemblytics_uniq_anchor.py --delta $parent2_based_prefix.delta --unique-length $unique_length --keep-small-uniques --out $parent2_based_prefix.1delta
    #gunzip -c $parent1_based_prefix.1delta.Assemblytics.unique_length_filtered_l${unique_length}.delta.gz >$parent1_based_prefix.1delta
    #gunzip -c $parent2_based_prefix.1delta.Assemblytics.unique_length_filtered_l${unique_length}.delta.gz >$parent2_based_prefix.1delta

    $mummer3_dir/show-coords -T -r -c -l -d $parent1_based_prefix.1delta > $parent1_based_prefix.1coords
    $mummer3_dir/show-coords -T -r -c -l -d $parent2_based_prefix.1delta > $parent2_based_prefix.1coords

    $mummer3_dir/show-snps -Clr -T $parent1_based_prefix.1delta > $parent1_based_prefix.snps
    $mummer3_dir/show-snps -Clr -T $parent2_based_prefix.1delta > $parent2_based_prefix.snps

    perl $RECOMBINEX_HOME/scripts/mummer2vcf.pl -r ${parent1_tag}.genome.fa -i $parent1_based_prefix.snps -t BOTH -p $parent1_based_prefix
    perl $RECOMBINEX_HOME/scripts/mummer2vcf.pl -r ${parent2_tag}.genome.fa -i $parent2_based_prefix.snps -t BOTH -p $parent2_based_prefix

    awk '{printf("##contig=<ID=%s,length=%d>\n",$1,$2);}' $parent1_tag.genome.fa.fai > $parent1_based_prefix.vcf_header.sequences.txt
    awk '{printf("##contig=<ID=%s,length=%d>\n",$1,$2);}' $parent2_tag.genome.fa.fai > $parent2_based_prefix.vcf_header.sequences.txt
    sed -i -e "/##reference/r $parent1_based_prefix.vcf_header.sequences.txt" $parent1_based_prefix.mummer2vcf.BOTH.vcf
    sed -i -e "/##reference/r $parent2_based_prefix.vcf_header.sequences.txt" $parent2_based_prefix.mummer2vcf.BOTH.vcf
    mv  $parent1_based_prefix.vcf_header.sequences.txt $parent1_based_output_dir
    mv  $parent2_based_prefix.vcf_header.sequences.txt $parent2_based_output_dir

    $vt_dir/vt sort $parent1_based_prefix.mummer2vcf.BOTH.vcf -o $parent1_based_prefix.mummer2vcf.BOTH.sort.vcf
    $vt_dir/vt sort $parent2_based_prefix.mummer2vcf.BOTH.vcf -o $parent2_based_prefix.mummer2vcf.BOTH.sort.vcf

    cat $parent1_based_prefix.mummer2vcf.BOTH.sort.vcf | egrep "^#"  > $parent1_based_prefix.vcf_header.txt
    cat $parent2_based_prefix.mummer2vcf.BOTH.sort.vcf | egrep "^#"  > $parent2_based_prefix.vcf_header.txt

    parent1_based_vcf_file_list+="$parent1_based_prefix.mummer2vcf.BOTH.sort.vcf "
    parent2_based_vcf_file_list+="$parent2_based_prefix.mummer2vcf.BOTH.sort.vcf "

    mv $parent1_based_prefix.delta $parent1_based_output_dir
    mv $parent1_based_prefix.1delta $parent1_based_output_dir
    #mv $parent1_based_prefix.1delta.* $parent1_based_output_dir
    mv $parent1_based_prefix.1coords $parent1_based_output_dir
    gzip $parent1_based_prefix.snps
    mv $parent1_based_prefix.snps.gz $parent1_based_output_dir
    gzip $parent1_based_prefix.mummer2vcf.BOTH.vcf
    mv $parent1_based_prefix.mummer2vcf.BOTH.vcf.gz $parent1_based_output_dir
    
    mv $parent2_based_prefix.delta $parent2_based_output_dir
    mv $parent2_based_prefix.1delta $parent2_based_output_dir
    #mv $parent2_based_prefix.1delta.* $parent2_based_output_dir
    mv $parent2_based_prefix.1coords $parent2_based_output_dir
    gzip $parent2_based_prefix.snps
    mv $parent2_based_prefix.snps.gz $parent2_based_output_dir
    gzip $parent2_based_prefix.mummer2vcf.BOTH.vcf
    mv $parent2_based_prefix.mummer2vcf.BOTH.vcf.gz $parent2_based_output_dir
fi

parent1_based_prefix="${parent1_tag}-${parent2_tag}.${parent1_tag}"
parent2_based_prefix="${parent1_tag}-${parent2_tag}.${parent2_tag}"

cat $parent1_based_vcf_file_list | egrep -v "^#" > $parent1_based_prefix.mummer2vcf.BOTH.sort.noheader.vcf;
cat $parent2_based_vcf_file_list | egrep -v "^#" > $parent2_based_prefix.mummer2vcf.BOTH.sort.noheader.vcf;

cat $parent1_based_prefix.vcf_header.txt $parent1_based_prefix.mummer2vcf.BOTH.sort.noheader.vcf > $parent1_based_prefix.mummer2vcf.BOTH.sort.withheader.vcf
cat $parent2_based_prefix.vcf_header.txt $parent2_based_prefix.mummer2vcf.BOTH.sort.noheader.vcf > $parent2_based_prefix.mummer2vcf.BOTH.sort.withheader.vcf

mv $parent1_based_vcf_file_list $parent1_based_output_dir
mv $parent2_based_vcf_file_list $parent2_based_output_dir

$vt_dir/vt decompose_blocksub $parent1_based_prefix.mummer2vcf.BOTH.sort.withheader.vcf -a -o $parent1_based_prefix.mummer2vcf.BOTH.decompose.vcf
$vt_dir/vt decompose_blocksub $parent2_based_prefix.mummer2vcf.BOTH.sort.withheader.vcf -a -o $parent2_based_prefix.mummer2vcf.BOTH.decompose.vcf

$vt_dir/vt normalize $parent1_based_prefix.mummer2vcf.BOTH.decompose.vcf -r $parent1_tag.genome.fa -f "~VARIANT_CONTAINS_N"| $vt_dir/vt uniq - -o  $parent1_based_prefix.mummer2vcf.BOTH.normalize.vcf
$vt_dir/vt normalize $parent2_based_prefix.mummer2vcf.BOTH.decompose.vcf -r $parent2_tag.genome.fa -f "~VARIANT_CONTAINS_N"| $vt_dir/vt uniq - -o  $parent2_based_prefix.mummer2vcf.BOTH.normalize.vcf

$vt_dir/vt annotate_variants $parent1_based_prefix.mummer2vcf.BOTH.normalize.vcf -r $parent1_tag.genome.fa -o $parent1_based_prefix.mummer2vcf.BOTH.annotate.vcf
$vt_dir/vt annotate_variants $parent2_based_prefix.mummer2vcf.BOTH.normalize.vcf -r $parent2_tag.genome.fa -o $parent2_based_prefix.mummer2vcf.BOTH.annotate.vcf

cat $parent1_based_prefix.mummer2vcf.BOTH.annotate.vcf | egrep "^#" > $parent1_based_prefix.mummer2vcf.BOTH.annotate.vcf.header
$bedtools_dir/bedtools subtract -a $parent1_based_prefix.mummer2vcf.BOTH.annotate.vcf -b $parent1_hardmask_bed > $parent1_based_prefix.mummer2vcf.BOTH.annotate_hardmask.vcf.content
cat $parent1_based_prefix.mummer2vcf.BOTH.annotate.vcf.header $parent1_based_prefix.mummer2vcf.BOTH.annotate_hardmask.vcf.content > $parent1_based_prefix.mummer2vcf.BOTH.annotate_hardmask.vcf
rm $parent1_based_prefix.mummer2vcf.BOTH.annotate.vcf.header
rm $parent1_based_prefix.mummer2vcf.BOTH.annotate_hardmask.vcf.content

cat $parent2_based_prefix.mummer2vcf.BOTH.annotate.vcf | egrep "^#" > $parent2_based_prefix.mummer2vcf.BOTH.annotate.vcf.header
$bedtools_dir/bedtools subtract -a $parent2_based_prefix.mummer2vcf.BOTH.annotate.vcf -b $parent2_hardmask_bed > $parent2_based_prefix.mummer2vcf.BOTH.annotate_hardmask.vcf.content
cat $parent2_based_prefix.mummer2vcf.BOTH.annotate.vcf.header $parent2_based_prefix.mummer2vcf.BOTH.annotate_hardmask.vcf.content > $parent2_based_prefix.mummer2vcf.BOTH.annotate_hardmask.vcf
rm $parent2_based_prefix.mummer2vcf.BOTH.annotate.vcf.header
rm $parent2_based_prefix.mummer2vcf.BOTH.annotate_hardmask.vcf.content

perl $RECOMBINEX_HOME/scripts/filter_vcf_by_window.pl -i $parent1_based_prefix.mummer2vcf.BOTH.annotate_hardmask.vcf -w $cluster_window_size -o $parent1_based_prefix.mummer2vcf.BOTH.annotate_hardmask.thin.vcf
perl $RECOMBINEX_HOME/scripts/filter_vcf_by_window.pl -i $parent2_based_prefix.mummer2vcf.BOTH.annotate_hardmask.vcf -w $cluster_window_size -o $parent2_based_prefix.mummer2vcf.BOTH.annotate_hardmask.thin.vcf

$vt_dir/vt view $parent1_based_prefix.mummer2vcf.BOTH.annotate_hardmask.thin.vcf  -f "VTYPE==SNP" -o $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.raw.vcf
$vt_dir/vt view $parent1_based_prefix.mummer2vcf.BOTH.annotate_hardmask.thin.vcf  -f "VTYPE==INDEL" -o $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.raw.vcf

$vt_dir/vt view $parent2_based_prefix.mummer2vcf.BOTH.annotate_hardmask.thin.vcf  -f "VTYPE==SNP" -o $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.raw.vcf
$vt_dir/vt view $parent2_based_prefix.mummer2vcf.BOTH.annotate_hardmask.thin.vcf  -f "VTYPE==INDEL" -o $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.raw.vcf

# vt filter is needed here because the original ref-query coordinate corresponding map for indels will be disrupted during vt decomposition and normalization, so for now, it is safe to eliminate these affected markers to reduce false positive rate, especially for INDEL markers. 
cat $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.raw.vcf | egrep -v "OLD_VARIANT=|OLD_CLUMPED=" > $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.VT_filter.vcf
cat $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.raw.vcf | egrep -v "OLD_VARIANT=|OLD_CLUMPED=" > $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.VT_filter.vcf
cat $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.raw.vcf | egrep -v "OLD_VARIANT=|OLD_CLUMPED=" > $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.VT_filter.vcf
cat $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.raw.vcf | egrep -v "OLD_VARIANT=|OLD_CLUMPED=" > $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.VT_filter.vcf

gzip $parent1_based_prefix.mummer2vcf.BOTH.*.vcf
mv $parent1_based_prefix.mummer2vcf.BOTH.*.vcf.gz $parent1_based_output_dir
gzip $parent2_based_prefix.mummer2vcf.BOTH.*.vcf
mv $parent2_based_prefix.mummer2vcf.BOTH.*.vcf.gz $parent2_based_output_dir

if [[ $genome_aln_mode == "by_chromosome" ]] && [[ $use_centromere_annotation == "yes" ]]
then
    # filter marker pairs locate on the opposite arms of the same chromosomes
    perl $RECOMBINEX_HOME/scripts/filter_markers_from_opposite_chromosome_arms.pl \
	 -i $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.VT_filter.vcf \
	 -rt $parent1_tag -qt $parent2_tag \
	 -rcg $parent1_centromere_gff -qcg $parent2_centromere_gff \
	 -o $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.ChrArm_filter.vcf
    perl $RECOMBINEX_HOME/scripts/filter_markers_from_opposite_chromosome_arms.pl \
	 -i $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.VT_filter.vcf \
	 -rt $parent1_tag -qt $parent2_tag \
	 -rcg $parent1_centromere_gff -qcg $parent2_centromere_gff \
	 -o $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.ChrArm_filter.vcf

    perl $RECOMBINEX_HOME/scripts/filter_markers_from_opposite_chromosome_arms.pl \
	 -i $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.VT_filter.vcf \
	 -rt $parent2_tag -qt $parent1_tag \
	 -rcg $parent2_centromere_gff -qcg $parent1_centromere_gff \
	 -o $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.ChrArm_filter.vcf
    perl $RECOMBINEX_HOME/scripts/filter_markers_from_opposite_chromosome_arms.pl \
	 -i $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.VT_filter.vcf \
	 -rt $parent2_tag -qt $parent1_tag \
	 -rcg $parent2_centromere_gff -qcg $parent1_centromere_gff \
	 -o $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.ChrArm_filter.vcf
    
    # apply the reciprocal filter
    perl $RECOMBINEX_HOME/scripts/filter_nonreciprocal_markers.pl \
	-g1_based_vcf $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.ChrArm_filter.vcf \
	-g2_based_vcf $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.ChrArm_filter.vcf \
	-g1 ${parent1_tag}.genome.fa \
	-g2 ${parent2_tag}.genome.fa \
	-g1_based_prefix  $parent1_based_prefix.final.SNP \
	-g2_based_prefix  $parent2_based_prefix.final.SNP 
    # perl $RECOMBINEX_HOME/scripts/filter_nonreciprocal_markers.pl \
    # 	-g1_based_vcf $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.ChrArm_filter.vcf \
    # 	-g2_based_vcf $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.ChrArm_filter.vcf \
    # 	-g1 ${parent1_tag}.genome.fa \
    # 	-g2 ${parent2_tag}.genome.fa \
    # 	-g1_based_prefix $parent1_based_prefix.final.INDEL \
    # 	-g2_based_prefix $parent2_based_prefix.final.INDEL
    
    gzip $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.ChrArm_filter.vcf
    gzip $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.ChrArm_filter.vcf
    gzip $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.ChrArm_filter.vcf
    gzip $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.ChrArm_filter.vcf
    mv $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.ChrArm_filter.vcf.gz $parent1_based_output_dir
    mv $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.ChrArm_filter.vcf.gz $parent1_based_output_dir
    mv $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.ChrArm_filter.vcf.gz $parent2_based_output_dir
    mv $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.ChrArm_filter.vcf.gz $parent2_based_output_dir

else
    # apply the reciprocal filter
    perl $RECOMBINEX_HOME/scripts/filter_nonreciprocal_markers.pl \
	-g1_based_vcf $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.VT_filter.vcf \
	-g2_based_vcf $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.VT_filter.vcf \
	-g1 ${parent1_tag}.genome.fa \
	-g2 ${parent2_tag}.genome.fa \
	-g1_based_prefix $parent1_based_prefix.final.SNP \
	-g2_based_prefix $parent2_based_prefix.final.SNP

    # perl $RECOMBINEX_HOME/scripts/filter_nonreciprocal_markers.pl \
    # 	-g1_based_vcf $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.VT_filter.vcf \
    # 	-g2_based_vcf $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.VT_filter.vcf \
    # 	-g1 ${parent1_tag}.genome.fa \
    # 	-g2 ${parent2_tag}.genome.fa \
    # 	-g1_based_prefix $parent1_based_prefix.final.INDEL \
    # 	-g2_based_prefix $parent2_based_prefix.final.INDEL

fi

gzip *.vcf
mv $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.raw.vcf.gz $parent1_based_output_dir
mv $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.raw.vcf.gz $parent1_based_output_dir
mv $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.raw.vcf.gz $parent2_based_output_dir
mv $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.raw.vcf.gz $parent2_based_output_dir

mv $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.VT_filter.vcf.gz $parent1_based_output_dir
mv $parent1_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.VT_filter.vcf.gz $parent1_based_output_dir
mv $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.SNP.VT_filter.vcf.gz $parent2_based_output_dir
mv $parent2_based_prefix.mummer2vcf.annotate_hardmask.thin.INDEL.VT_filter.vcf.gz $parent2_based_output_dir

gzip *.markers.txt
gzip $parent1_based_output_dir/*.vcf
gzip $parent2_based_output_dir/*.vcf

# calculate intermarker distance
perl $RECOMBINEX_HOME/scripts/intermarker_distance.pl -i $parent1_based_prefix.final.SNP.markers.txt.gz \
     -p $parent1_based_prefix.final.SNP.markers
# perl $RECOMBINEX_HOME/scripts/intermarker_distance.pl -i $parent1_based_prefix.final.INDEL.markers.txt.gz \
#      -p $parent1_based_prefix.final.INDEL.markers
# gunzip -c $parent1_based_prefix.final.INDEL.markers.txt.gz | tail -n +2 | gzip > $parent1_based_prefix.final.INDEL.markers.content.txt.gz 
# cat $parent1_based_prefix.final.SNP.markers.txt.gz $parent1_based_prefix.final.INDEL.markers.content.txt.gz > $parent1_based_prefix.final.BOTH.markers.txt.gz
# rm $parent1_based_prefix.final.INDEL.markers.content.txt.gz
# perl $RECOMBINEX_HOME/scripts/intermarker_distance.pl -i $parent1_based_prefix.final.BOTH.markers.txt.gz \
#     -p $parent1_based_prefix.final.BOTH.markers

perl $RECOMBINEX_HOME/scripts/intermarker_distance.pl -i $parent2_based_prefix.final.SNP.markers.txt.gz \
     -p $parent2_based_prefix.final.SNP.markers
# perl $RECOMBINEX_HOME/scripts/intermarker_distance.pl -i $parent2_based_prefix.final.INDEL.markers.txt.gz \
#      -p $parent2_based_prefix.final.INDEL.markers
# gunzip -c $parent2_based_prefix.final.INDEL.markers.txt.gz | tail -n +2 | gzip > $parent2_based_prefix.final.INDEL.markers.content.txt.gz 
# cat $parent2_based_prefix.final.SNP.markers.txt.gz $parent2_based_prefix.final.INDEL.markers.content.txt.gz > $parent2_based_prefix.final.BOTH.markers.txt.gz
# rm $parent2_based_prefix.final.INDEL.markers.content.txt.gz
# perl $RECOMBINEX_HOME/scripts/intermarker_distance.pl -i $parent2_based_prefix.final.BOTH.markers.txt.gz \
#     -p $parent2_based_prefix.final.BOTH.markers

cp $parent1_based_prefix.final.*.markers.vcf.gz ./$parent1_based_output_dir/
cp $parent2_based_prefix.final.*.markers.vcf.gz ./$parent2_based_output_dir/

# plot markers

if [[ $use_centromere_annotation == "yes" ]]
then
    Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent1_fai \
	--marker_table $parent1_based_prefix.final.SNP.markers.txt.gz \
	--centromere_gff $parent1_centromere_gff --output_prefix $parent1_based_prefix.final.SNP.markers
    rm Rplots.pdf
    # Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent1_fai \
    # 	--marker_table $parent1_based_prefix.final.INDEL.markers.txt.gz \
    # 	--centromere_gff $parent1_centromere_gff --output_prefix $parent1_based_prefix.final.INDEL.markers
    # rm Rplots.pdf
    # Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent1_fai \
    # 	--marker_table $parent1_based_prefix.final.BOTH.markers.txt.gz \
    # 	--centromere_gff $parent1_centromere_gff --output_prefix $parent1_based_prefix.final.BOTH.markers
    # rm Rplots.pdf

    Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent2_fai \
	--marker_table $parent2_based_prefix.final.SNP.markers.txt.gz \
	--centromere_gff $parent2_centromere_gff --output_prefix $parent2_based_prefix.final.SNP.markers
    rm Rplots.pdf
    # Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent2_fai \
    # 	--marker_table $parent2_based_prefix.final.INDEL.markers.txt.gz \
    # 	--centromere_gff $parent2_centromere_gff --output_prefix $parent2_based_prefix.final.INDEL.markers
    # rm Rplots.pdf
    # Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent2_fai \
    # 	--marker_table $parent2_based_prefix.final.BOTH.markers.txt.gz \
    # 	--centromere_gff $parent2_centromere_gff --output_prefix $parent2_based_prefix.final.BOTH.markers
    # rm Rplots.pdf

else
    Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent1_fai \
	--marker_table $parent1_based_prefix.final.SNP.markers.txt.gz \
	--output_prefix $parent1_based_prefix.final.SNP.markers
    rm Rplots.pdf
    # Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent1_fai \
    # 	--marker_table $parent1_based_prefix.final.INDEL.markers.txt.gz \
    # 	--output_prefix $parent1_based_prefix.final.INDEL.markers
    # rm Rplots.pdf
    # Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R --genome_fai $parent1_fai \
    # 	--marker_table $parent1_based_prefix.final.BOTH.markers.txt.gz \
    # 	--output_prefix $parent1_based_prefix.final.BOTH.markers
    # rm Rplots.pdf

    Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent2_fai \
	--marker_table $parent2_based_prefix.final.SNP.markers.txt.gz \
	--output_prefix $parent2_based_prefix.final.SNP.markers
    rm Rplots.pdf
    # Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent2_fai \
    # 	--marker_table $parent2_based_prefix.final.INDEL.markers.txt.gz \
    # 	--output_prefix $parent2_based_prefix.final.INDEL.markers
    # rm Rplots.pdf
    # Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_markers.R  --genome_fai $parent2_fai \
    # 	--marker_table $parent2_based_prefix.final.BOTH.markers.txt.gz \
    # 	--output_prefix $parent2_based_prefix.final.BOTH.markers
    # rm Rplots.pdf
fi    


# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm *.vcf_header.txt
    rm $parent1_tag.genome.fa
    rm $parent2_tag.genome.fa
    rm $parent1_tag.genome.fa.fai
    rm $parent2_tag.genome.fa.fai
fi



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


