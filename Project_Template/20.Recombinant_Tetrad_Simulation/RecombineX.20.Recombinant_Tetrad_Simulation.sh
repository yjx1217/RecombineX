#!/bin/bash
set -e -o pipefail

##########################################
# load environment variables for RecombineX
source ./../../env.sh

###########################################
# set project-specific variables

coordinate_genome="$RECOMBINEX_HOME/data/ref.genome.raw.relabel.fa" # The input genome file for set up the coordinate system (FASTA format). Default = "$RECOMBINEX_HOME/data/ref.genome.raw.relabel.fa".
parent1_tag="P1"; # The genome tag for simulated genome1. Default = "P1".
parent2_tag="P2"; # The genome tag for simulated genome2. Default = "P2".
marker_table_file="$RECOMBINEX_HOME/data/P1-P2.ref.final.SNP.markers.txt"; # The SNP markers based on the input genome. Default = "$RECOMBINEX_HOME/data/P1-P2.ref.final.SNP.markers.txt".
co_num=90; # The number of crossover (CO) per tetrad. Default = 90. # (mean = 90.5 for S. cerevisiae based on Mancera et al. 2008 Nature paper).
nco_num=65; # The number of noncrossover (NCO) per tetrad. Default = 65. # (mean = 66.1 for S. cerevisiae based on Mancera et al. 2008 Nature paper).
min_interevent_distance=10000; # The minimal genomic distance in basepairs (bp) allowed between two independent recombination events. Default = 10000. (i.e. 10 kb). 
linked_region_bed="";  # genomic regions for complete linkeage (i.e. no recombination) in 3-column BED format (without header).
gc_by_co_ratio=1; # The ratio of the number of crossover with associated gene conversions (GC) among all crossover events. 
co_gc_size_mean=2500; # The mean of crossover-associated gene conversion (CO-GC) size in basepairs (bp). Default = 2500. ("2461.836" for S. cerevisiae based on Mancera et al. 2008 Nature paper).
co_gc_size_stdev=2000; # The standard deviation of crossover-associated gene conversion (CO-GC) size in basepair (bp). Default = 2000. ("2057.673" for S. cerevisiae based on Mancera et al. 2008 Nature paper).
nco_gc_size_mean=2250; # The mean of noncrossover-associated gene conversion (NCO-GC) size in basepairs (bp). Default = 2250. ("2247.716" for NCO-GC for S. cerevisiae based on Mancera et al. 2008 Nature paper).
nco_gc_size_stdev=2200; # The standard deviation of noncrossover-associated gene conversion (NCO-GC) size in basepair (bp). Default = 2200. ("2173.46" for NCO-GC for S. cerevisiae based on Mancera et al. 2008 Nature paper).
min_gc_size=100; # Safe lower bound gene conversion (GC) size in basepair (bp). Default = 100. 
max_gc_size=5000; # Safe upper bound gene conversion (GC) size in basepair (bp). Default = 5000. 
output_prefix="P1-P2"; # The prefix for outputs. Default = "P1-P2".
random_seed="20210210" # The random seed integer for tetrad and read simulation. Default = "20210210". 

## parameters for genotype plotting
color_scheme="$RECOMBINEX_HOME/data/Saccharomyces_cerevisiae.color_scheme.txt" # The color scheme to use for plotting genotypes. Default = "$RECOMBINEX_HOME/data/Saccharomyces_cerevisiae.color_scheme.txt".
plot_centromere="yes" # Whether to plot centromere in the generated genotyping plots. Please note that enable this option requires that you also define the centromere annotation GFF3 file for the specified coordinate genome. Default = "yes".
centromere_gff_for_coordinate_genome="$RECOMBINEX_HOME/data/ref.centromere.relabel.gff" # The centromere annotation GFF3 file for the specified coordinate genome. Only needed when plot_centromere="yes". Default = "$RECOMBINEX_HOME/data/ref.centromere.relabel.gff".

## parameters for reads simulating
Illumina_platform_and_read_length="HiSeq2500L150" # The Illumina sequencing platform and the associated read length. Available options include: "HiSeq2kL100" for HiSeq2000 platform with 100-bp read length; "HiSeq2500L150" for HiSeq2500 platform with 150-bp read length; "HiSeq2500L150" for HiSeq2500 platform with 125-bp read length. "HiSeqXtruSeqL150" for HiSeqX TrueSeq platform with 150-bp read length; "MiSeqv3L250" for MiSeq platform with 250-bp read length; and "NextSeq500v2L75" for NextSeq platform with 75-bp read length. Default = "HiSeq2500L150".

read_coverage="30" # The coverage of simulated parent and gamete reads. Default = "30".

debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
###########################################

# process the pipeline

echo ""
test_file_existence () {
    filename=$1
    if [[ ! -f $filename ]]
    then
        echo "The file $filename does not exists! Process terminated!"
        exit
    fi
}

random_seed_phrase_for_tetrad_simulation="$random_seed"; # The random seed phrase for tetrad simulation. Default = "simulated_tetrad".
random_seed_integer_for_read_simulation="$random_seed"; # The random seed number for read simulation. Default = "20190518".

perl $RECOMBINEX_HOME/scripts/simulate_recombinant_tetrad.pl \
    -cg $coordinate_genome \
    -p1 $parent1_tag \
    -p2 $parent2_tag \
    -m $marker_table_file \
    -co $co_num \
    -nco $nco_num \
    -gc_by_co_ratio $gc_by_co_ratio \
    -co_gc_size_mean $co_gc_size_mean \
    -co_gc_size_stdev $co_gc_size_stdev \
    -nco_gc_size_mean $nco_gc_size_mean \
    -nco_gc_size_stdev $nco_gc_size_stdev \
    -min_gc_size $min_gc_size \
    -max_gc_size $max_gc_size \
    -d $min_interevent_distance \
    -l $linked_region_bed \
    -p $output_prefix \
    -s $random_seed_phrase_for_tetrad_simulation

cp $coordinate_genome coordinate_genome.fa
$samtools_dir/samtools faidx coordinate_genome.fa

Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_tetrad_genotype_for_simulation.R  \
    --input $output_prefix.simulated_tetrad.genotype.for_genotype_plotting.txt \
    --color_scheme $color_scheme \
    --genome1_tag $parent1_tag \
    --genome2_tag $parent2_tag \
    --coordinate_genome_fai coordinate_genome.fa.fai \
    --plot_centromere $plot_centromere \
    --centromere_gff_for_coordinate_genome $centromere_gff_for_coordinate_genome \
    --output $output_prefix.simulated_tetrad.genotype.plot.pdf 

read1_profile="$art_dir/Illumina_profiles/${Illumina_platform_and_read_length}R1.txt"
read2_profile="$art_dir/Illumina_profiles/${Illumina_platform_and_read_length}R2.txt"
test_file_existence $read1_profile
test_file_existence $read2_profile
echo "Simulating Illumina reads with the following read profile:"
echo "$read1_profile"
echo "$read2_profile"
read_length=${Illumina_platform_and_read_length##*L}
mean_fragment_size=500
echo "read_length=$read_length, mean_fragment_size=$mean_fragment_size" 

for parent_tag in $parent1_tag $parent2_tag
do
    parent_genome="$parent_tag.genome.fa"
    test_file_existence $parent_genome
    echo "Simulating Illumina reads for the genome $parent_genome .." 
    echo ""
    $art_dir/art_illumina \
        --qprof1 $read1_profile \
        --qprof2 $read2_profile \
        -f $read_coverage \
        -l $read_length \
        -i $parent_genome \
        -p \
        -na \
        -rs $random_seed_integer_for_read_simulation \
        -m $mean_fragment_size \
        -s 10 \
        -o $parent_tag.simulated_reads.${read_coverage}X.R

    echo ""
    echo "Compress simulated reads for the genome $parent_genome .."
    gzip $parent_tag.simulated_reads.${read_coverage}X.R1.fq
    gzip $parent_tag.simulated_reads.${read_coverage}X.R2.fq
done
echo ""

for spore in a b c d
do
    spore_genome="$output_prefix.simulated_tetrad.$spore.genome.fa" 
    test_file_existence $spore_genome
    echo "Simulating Illumina reads for the genome $spore_genome .." 
    echo ""
    $art_dir/art_illumina \
        --qprof1 $read1_profile \
        --qprof2 $read2_profile \
        -f $read_coverage \
        -l $read_length \
        -i $spore_genome \
        -p \
        -na \
        -rs $random_seed_integer_for_read_simulation \
        -m $mean_fragment_size \
        -s 10 \
        -o $output_prefix.simulated_tetrad.$spore.simulated_reads.${read_coverage}X.R
    echo ""
    echo "Compress simulated reads for the genome $spore_genome .."
    gzip $output_prefix.simulated_tetrad.$spore.simulated_reads.${read_coverage}X.R1.fq
    gzip $output_prefix.simulated_tetrad.$spore.simulated_reads.${read_coverage}X.R2.fq
done

echo ""
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
