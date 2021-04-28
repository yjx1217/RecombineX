#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Env;
use Cwd;

##############################################################
#  script: batch_read_mapping_to_reference_genome.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2021.06.17
#  description: run batch read mapping to the reference genome
#  example: perl batch_read_mapping_to_reference_genome.pl -i Master_Sample_Table.txt -t 4 -o output_dir -reference_genome_assembly_dir ./../01.Reference_Genome_Preprocessing -gamete_reads_dir ./../00.Gamete_Reads -min_mappability 0.85 -window_size 250 -ploidy 1
##############################################################

my $RECOMBINEX_HOME = $ENV{RECOMBINEX_HOME};
my $java_dir = $ENV{java_dir};
my $trimmomatic_dir = $ENV{trimmomatic_dir};
my $bwa_dir = $ENV{bwa_dir};
my $samtools_dir = $ENV{samtools_dir};
my $picard_dir = $ENV{picard_dir};
my $gatk3_dir = $ENV{gatk3_dir};
my $bedtools_dir = $ENV{bedtools_dir};
my $gemtools_dir = $ENV{gemtools_dir};
my $freec_dir = $ENV{freec_dir};
my $sample_table = "Master_Sample_Table.txt";
my $threads = 1;
my $mapping_quality_cutoff_for_mpileup = 30;
my $reference_genome_assembly_dir;
my $gamete_reads_dir;
my $output_dir = "OUTPUT_DIR";
# for aneuploidy detection
my $min_mappability = 0.85;
my $window_size = 250;
my $step_size = 250;
my $ploidy = 1;
my $excluded_chr_list_for_cnv_profiling;
my $debug = "no";

GetOptions('sample_table|i:s' => \$sample_table,
	   'threads|t:i' => \$threads,
	   'mapping_quality_cutoff_for_mpileup|q:i' => \$mapping_quality_cutoff_for_mpileup,
	   'reference_genome_assembly_dir|rga_dir:s' => \$reference_genome_assembly_dir,
	   'gamete_reads_dir|or_dir:s' => \$gamete_reads_dir,
	   'output_dir|o_dir:s' => \$output_dir,
	   'min_mappability|min_map:f' => \$min_mappability,
	   'window|w:i' => \$window_size,
	   'step|s:i' => \$step_size,
	   'ploidy|p:i' => \$ploidy,
           'excluded_chr_list_for_cnv_profiling|e:s' => \$excluded_chr_list_for_cnv_profiling,
           'debug|d:s' => \$debug); 

my $sample_table_fh = read_file($sample_table);
my %sample_table = ();
my @sample_table = ();
parse_sample_table_by_spore($sample_table_fh, \%sample_table, \@sample_table);
my $base_dir = cwd();
system("mkdir $output_dir");

my $all_samples_cnv = "$base_dir/$output_dir/all_samples.CNV_summary.txt";
my $all_samples_cnv_fh = write_file($all_samples_cnv);
print $all_samples_cnv_fh "sample\tref_genome\tchr\tstart\tend\tcopy_number\tstatus\tMWU_test_p_value\tKS_test_p_value\n";

my $sample_count = 0;
my $mapping_count = 0;
my $adapter = "$trimmomatic_dir/adapters/TruSeq3-PE-2.fa";

foreach my $sample_id (@sample_table) {
    my $local_time = localtime();
    print "[$local_time] processing sample $sample_id with short-read mapping\n";
    $sample_count++;
    my $R1_read_file = "$base_dir/$gamete_reads_dir/$sample_table{$sample_id}{'R1_read_file'}";
    my $R2_read_file = "$base_dir/$gamete_reads_dir/$sample_table{$sample_id}{'R2_read_file'}";
    print "processing sample $sample_id\n";
    print "PE_reads_file = $R1_read_file,$R2_read_file\n";
    print "Check the specified short read file:\n";
    if (-e $R1_read_file) {
        print "Successfully located the specified short read file 1: $R1_read_file.\n";
    } else {
        print "Cannot find the specified short read file 1: $R1_read_file!\n";
        print "Exit!\n";
        exit;
    }
    if (-e $R2_read_file) {
        print "Successfully located the specified short read file 2: $R2_read_file.\n";
    } else {
        print "Cannot find the specified short read file 2: $R2_read_file!\n";
        print "Exit!\n";
        exit;
    }

    my $excluded_chr_list_for_cnv_profiling_file = "$base_dir/$excluded_chr_list_for_cnv_profiling";
    print "Check the specified excluded_chr_list_for_cnv_profiling file:\n";
    if (-e $excluded_chr_list_for_cnv_profiling_file) {
        print "Successfully located the specified excluded_chr_list_for_cnv_profiling_file: $excluded_chr_list_for_cnv_profiling_file.\n";
    } else {
        print "Cannot find the specified excluded_chr_list_for_cnv_profiling_file: $excluded_chr_list_for_cnv_profiling_file!\n";
        print "Exit!\n";
        exit;
    }

    # my $R1_read_fh = read_file($R1_read_file);
    # my $raw_read_length = get_read_length($R1_read_fh);
    my $raw_read_length = 100;
    print "raw read length = $raw_read_length\n";
    my $cross_pair = $sample_table{$sample_id}{'cross_pair'};
    my $tetrad_id = $sample_table{$sample_id}{'tetrad_id'};
    my $spore_id = $sample_table{$sample_id}{'spore_id'};
    my $sample_tag = "$cross_pair.$tetrad_id.$spore_id";
    my ($genome1_tag, $genome2_tag) = split /-/, $cross_pair;
    print "genome1: $genome1_tag, genome2: $genome2_tag\n";
    # my @parent_genome_tags = ($genome1_tag, $genome2_tag);
    my @reference_genome_tag = ("ref");
    my $sample_output_dir;
    foreach my $reference_genome_tag (@reference_genome_tag) {
	print "mapping to $reference_genome_tag\n";
	$mapping_count++;
	my $reference_genome_assembly = "$base_dir/$reference_genome_assembly_dir/${reference_genome_tag}.genome.raw.relabel.fa";
	if (-e $reference_genome_assembly) {
            print "Successfully located the expected reference_genome_assembly file: $reference_genome_assembly.\n";
        } else {
            print "Cannot find the expected reference genome assembly file: $reference_genome_assembly!\n";
            print "Exit!\n";
            exit;
        }

	$sample_output_dir = "$base_dir/$output_dir/${sample_tag}.${reference_genome_tag}";
	system("mkdir -p  $sample_output_dir");
	chdir("$sample_output_dir") or die "cannot change directory to: $!\n";
	system("mkdir tmp");
	print "trim the reads by trimmomatic\n";
	system("ln -s $adapter adapter.fa");
	system("ln -s $R1_read_file $sample_tag.R1.raw.fq.gz");
	system("ln -s $R2_read_file $sample_tag.R2.raw.fq.gz");
	system("$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $trimmomatic_dir/trimmomatic.jar PE -threads $threads -phred33 $sample_tag.R1.raw.fq.gz  $sample_tag.R2.raw.fq.gz $sample_tag.R1.trimmed.PE.fq.gz $sample_tag.R1.trimmed.SE.fq.gz  $sample_tag.R2.trimmed.PE.fq.gz $sample_tag.R2.trimmed.SE.fq.gz  ILLUMINACLIP:adapter.fa:2:30:10  SLIDINGWINDOW:5:20 MINLEN:36");
	if ($debug eq "no") {
	    system("rm $sample_tag.R1.raw.fq.gz");
	    system("rm $sample_tag.R2.raw.fq.gz");
	    system("rm $sample_tag.R1.trimmed.SE.fq.gz");
	    system("rm $sample_tag.R2.trimmed.SE.fq.gz");
            system("rm adapter.fa");
	}
	print("mapping the reads by bwa\n");
	system("ln -s $reference_genome_assembly $reference_genome_tag.genome.fa");
	system("$bwa_dir/bwa index $reference_genome_tag.genome.fa");
	system("$bwa_dir/bwa mem -t $threads -M $reference_genome_tag.genome.fa $sample_tag.R1.trimmed.PE.fq.gz $sample_tag.R2.trimmed.PE.fq.gz | $samtools_dir/samtools view -bS -q $mapping_quality_cutoff_for_mpileup -F 3340 -f 2 - >${sample_tag}.${reference_genome_tag}.bam");
	if ($debug eq "no") {
            system("rm $reference_genome_tag.genome.fa.bwt");
            system("rm $reference_genome_tag.genome.fa.pac");
            system("rm $reference_genome_tag.genome.fa.ann");
            system("rm $reference_genome_tag.genome.fa.amb");
            system("rm $reference_genome_tag.genome.fa.sa");
            system("rm $sample_tag.R1.trimmed.PE.fq.gz");
            system("rm $sample_tag.R2.trimmed.PE.fq.gz");
        }
	## index reference sequence
	system("$samtools_dir/samtools faidx $reference_genome_tag.genome.fa");
	system("$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar CreateSequenceDictionary -REFERENCE $reference_genome_tag.genome.fa -OUTPUT $reference_genome_tag.genome.dict");
	## sort bam file by picard-tools: SortSam
	system("$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar SortSam -INPUT ${sample_tag}.${reference_genome_tag}.bam -OUTPUT ${sample_tag}.${reference_genome_tag}.sort.bam -SORT_ORDER coordinate -MAX_RECORDS_IN_RAM 1000000");
	if ($debug eq "no") {
            system("rm ${sample_tag}.${reference_genome_tag}.bam");
        }
	## fixmate
	system("$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar FixMateInformation -INPUT ${sample_tag}.${reference_genome_tag}.sort.bam -OUTPUT ${sample_tag}.${reference_genome_tag}.fixmate.bam");
	if ($debug eq "no") {
            system("rm ${sample_tag}.${reference_genome_tag}.sort.bam");
        }
	## add or replace read groups and sort
	system("$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar AddOrReplaceReadGroups -INPUT=${sample_tag}.${reference_genome_tag}.fixmate.bam -OUTPUT ${sample_tag}.${reference_genome_tag}.rdgrp.bam -SORT_ORDER coordinate -MAX_RECORDS_IN_RAM 1000000 -RGID ${sample_tag}.${reference_genome_tag} -RGLB ${sample_tag}.${reference_genome_tag} -RGPL 'Illumina' -RGPU ${sample_tag}.${reference_genome_tag} -RGSM ${sample_tag} -RGCN 'RGCN'");
	if ($debug eq "no") {
            system("rm ${sample_tag}.${reference_genome_tag}.fixmate.bam");
        }
	# Picard tools remove duplicates
	system("$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar MarkDuplicates -INPUT ${sample_tag}.${reference_genome_tag}.rdgrp.bam -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 -REMOVE_DUPLICATES true -METRICS_FILE ${sample_tag}.${reference_genome_tag}.dedup.matrics -OUTPUT ${sample_tag}.${reference_genome_tag}.dedup.bam"); 
	# index the dedup.bam file
	system("$samtools_dir/samtools index ${sample_tag}.${reference_genome_tag}.dedup.bam");
	if ($debug eq "no") {
            system("rm ${sample_tag}.${reference_genome_tag}.rdgrp.bam");
        }
	# GATK local realign
	# find realigner targets
	system("$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $gatk3_dir/gatk3.jar -nt $threads -R $reference_genome_tag.genome.fa -T RealignerTargetCreator -I ${sample_tag}.${reference_genome_tag}.dedup.bam  -o ${sample_tag}.${reference_genome_tag}.realn.intervals");
	# run realigner
	system("$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $gatk3_dir/gatk3.jar -R $reference_genome_tag.genome.fa -T IndelRealigner -I ${sample_tag}.${reference_genome_tag}.dedup.bam -targetIntervals ${sample_tag}.${reference_genome_tag}.realn.intervals -o ${sample_tag}.${reference_genome_tag}.realn.bam");
	# index final bam file
	system("$samtools_dir/samtools index ${sample_tag}.${reference_genome_tag}.realn.bam");
	if ($debug eq "no") {
            system("rm ${sample_tag}.${reference_genome_tag}.dedup.bam");
            system("rm ${sample_tag}.${reference_genome_tag}.dedup.bam.bai");
            system("rm ${sample_tag}.${reference_genome_tag}.dedup.matrics");
            system("rm ${sample_tag}.${reference_genome_tag}.realn.intervals");
        }

	# generate samtools mpileup 
	system("$samtools_dir/samtools mpileup -Q 0 -C 0 -q $mapping_quality_cutoff_for_mpileup -f $reference_genome_tag.genome.fa  ${sample_tag}.${reference_genome_tag}.realn.bam |gzip -c >${sample_tag}.${reference_genome_tag}.mpileup.gz");	
	# compute basic alignment statistics by samtools
	system("$samtools_dir/samtools flagstat ${sample_tag}.${reference_genome_tag}.realn.bam >${sample_tag}.${reference_genome_tag}.samstat");
	# calculate per-base depth
	system("$samtools_dir/samtools depth -aa ${sample_tag}.${reference_genome_tag}.realn.bam |gzip -c >${sample_tag}.${reference_genome_tag}.depth.txt.gz");
	# compute insert size statistics
	# system("$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar CollectInsertSizeMetrics -I ${sample_tag}.${reference_genome_tag}.realn.bam -O ${sample_tag}.${reference_genome_tag}.insert_size_metrics.txt -H ${sample_tag}.${reference_genome_tag}.insert_size_histogram.pdf -M 0.5");
	# calculate read mapping coverage statistics
	system("perl $RECOMBINEX_HOME/scripts/summarize_mapping_coverage.pl -r $reference_genome_tag.genome.fa -s ${sample_tag}.${reference_genome_tag}.samstat -d ${sample_tag}.${reference_genome_tag}.depth.txt.gz -c 5 -t $sample_id -o ${sample_tag}.${reference_genome_tag}.coverage_summary.txt");
	if ($mapping_count == 1) {
	    system("cp ${sample_tag}.${reference_genome_tag}.coverage_summary.txt ./../all_samples.coverage_summary.txt");
	} else {
	    system("tail -1 ${sample_tag}.${reference_genome_tag}.coverage_summary.txt  >> ./../all_samples.coverage_summary.txt");
	}
	if ($debug eq "no") {
            system("rm ${sample_tag}.${reference_genome_tag}.depth.txt.gz");
        }
	# scan for aneuploidy with FREEC
	my $gc_range_file = "$base_dir/$reference_genome_assembly_dir/ref.FREEC.GC_range.txt";
	my $gc_range_fh = read_file($gc_range_file);
	my ($min_expected_gc, $max_expected_gc) = extract_gc_range($gc_range_fh);
	system("perl $RECOMBINEX_HOME/scripts/run_FREEC_wrapper_lite_for_RecombineX.pl -r $reference_genome_tag -bam ${sample_tag}.${reference_genome_tag}.realn.bam -prefix ${sample_tag}.${reference_genome_tag} -ploidy $ploidy -bedtools $bedtools_dir/bedtools -samtools $samtools_dir/samtools  -freec $freec_dir/freec -window $window_size -step $step_size -read_length_for_mappability $raw_read_length -min_mappability $min_mappability -min_expected_gc $min_expected_gc -max_expected_gc $max_expected_gc -mates_orientation 0 -threads $threads -refseq_genome_preprocessing_dir ./../../$reference_genome_assembly_dir -excluded_chr_list ./../../$excluded_chr_list_for_cnv_profiling");

	system("Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/CNV_segmentation_by_DNAcopy.R --input ${sample_tag}.${reference_genome_tag}.FREEC.bam_ratio.txt --prefix ${sample_tag}.${reference_genome_tag} --window $window_size --ploidy $ploidy --genome_fai $reference_genome_tag.genome.fa.fai");

	system("perl $RECOMBINEX_HOME/scripts/adjust_FREEC_copynumber_by_DNAcopy_copynumber.pl -i ${sample_tag}.${reference_genome_tag}.FREEC.bam_ratio.sorted.txt -a ${sample_tag}.${reference_genome_tag}.FREEC.bam_ratio.sorted.resegmented.lite.txt -o ${sample_tag}.${reference_genome_tag}.FREEC.bam_ratio.sorted.adjusted.txt");

	if ( -s "${sample_tag}.${reference_genome_tag}.FREEC.bam_ratio.sorted.adjusted.txt") {
	    system("Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_CNV_for_FREEC.R --ploidy $ploidy --genome_fai $reference_genome_tag.genome.fa.fai --input ${sample_tag}.${reference_genome_tag}.FREEC.bam_ratio.sorted.adjusted.txt --output ${sample_tag}.${reference_genome_tag}.CNV_plot.pdf");
	    system("rm Rplots.pdf");
	    if ( -s "${sample_tag}.${reference_genome_tag}.FREEC.bam_ratio.sorted.resegmented.CNVs.txt") {
		system("Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/assess_CNV_significance_for_FREEC.R --cnv ${sample_tag}.${reference_genome_tag}.FREEC.bam_ratio.sorted.resegmented.CNVs.txt --ratio ${sample_tag}.${reference_genome_tag}.FREEC.bam_ratio.sorted.adjusted.txt --genome_fai $reference_genome_tag.genome.fa.fai --output ${sample_tag}.${reference_genome_tag}.CNV_significance_test.txt");
		my $sample_cnv_fh = read_file("${sample_tag}.${reference_genome_tag}.CNV_significance_test.txt");
		parse_sample_cnv_file($sample_cnv_fh, $sample_id, "$reference_genome_tag.genome.fa", $all_samples_cnv_fh);
		close $sample_cnv_fh;
	    } else {
		system("echo -e \"chr\tstart\tend\tcopy_number\tstatus\tMWU_test_p_value\tKS_test_p_value\" > ${sample_tag}.${reference_genome_tag}.CNV_significance_test.txt");
	    }
	} else {
	    system(" echo \"Exception encountered for FREEC! Exit! ...\" > ${sample_tag}.${reference_genome_tag}.realn.bam.no_FREEC.txt");
	}
	if ($debug eq "no") {
	    # remove old files
            system("rm for_CNV.bam");
	    system("rm FREEC.bam");
	    system("rm FREEC.bam.header.old.sam");
	    system("rm FREEC.bam.header.new.sam");
	    system("rm GC_profile.${window_size}bp.cnp");
	    system("rm -r tmp");
	}
	$local_time = localtime();
	print "[$local_time] finishing short-read mapping for sample $sample_id \n";
	chdir("./../") or die "cannot change directory to: $!\n";
    }
}

close $all_samples_cnv_fh;

print "A total of $sample_count samples were processed!\n";


sub read_file {
    my $file = shift @_;
    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "gunzip -c $file |") or die "can't open pipe to $file";
    }
    else {
        open($fh, $file) or die "can't open $file";
    }
    return $fh;
}

sub write_file {
    my $file = shift @_;
    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "| gzip -c >$file") or die "can't open $file\n";
    } else {
        open($fh, ">$file") or die "can't open $file\n";
    }
    return $fh;
}  

sub parse_sample_table_by_spore {
    my ($fh, $sample_table_hashref, $sample_table_arrayref) = @_;
    while (<$fh>) {
        chomp;
        /^#/ and next;
	/^\s*$/ and next;
        my ($sample_id, $tetrad_id, $spore_id, $PE_reads_files, $cross_pair, $note) = split /\s+/, $_;
        push @$sample_table_arrayref, $sample_id;
        my ($R1_read_file, $R2_read_file) = split /,/, $PE_reads_files;
        $$sample_table_hashref{$sample_id}{'tetrad_id'} = $tetrad_id;
        $$sample_table_hashref{$sample_id}{'spore_id'} = $spore_id;
        $$sample_table_hashref{$sample_id}{'R1_read_file'} = $R1_read_file;
        $$sample_table_hashref{$sample_id}{'R2_read_file'} = $R2_read_file;
        $$sample_table_hashref{$sample_id}{'cross_pair'} = $cross_pair;
        $$sample_table_hashref{$sample_id}{'note'} = $note;
    }
}

sub parse_sample_cnv_file {
    my ($sample_cnv_fh, $sample, $refseq, $all_samples_cnv_fh) = @_;
    while (<$sample_cnv_fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	/chr\tstart/ and next;
	my ($chr, $start, $end, $copy_number, $status, $MWU_test_p_value, $KS_test_p_value) = split /\t/, $_;
	# if ($MWU_test_p_value < 0.05) {
        print $all_samples_cnv_fh "$sample\t$refseq\t$_\n";
        #}
    }
}


sub get_read_length {
    my $fh = shift @_;
    my $read_length = 100;
    while (<$fh>) {
        chomp;
        if ($. % 4 == 2) {
            $read_length = length $_;
            last;
        }
    }
    return $read_length;
}


sub extract_gc_range {
    my $fh = shift @_;
    my $min_expected_gc;
    my $max_expected_gc;
    while (<$fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	($min_expected_gc, $max_expected_gc) = split /\t/, $_;
    }

    return ($min_expected_gc, $max_expected_gc);
}
