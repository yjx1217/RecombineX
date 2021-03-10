#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Env;

##############################################################
#  script: batch_recombination_plotting_by_reference_genome.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.09.19
#  description: run batch recombination plotting to visualize different crossover and gene conversion events
#  example: perl batch_recombination_plotting_by_reference_genome.pl -s master_sample_table --genome_dir genome_dir --genotype_dir genotype_dir --batch_id batch_id --recombination_profiling_dir recombination_profiling_dir -q 30 -o output_dir -c color_scheme --plot_centromere yes
##############################################################

my $RECOMBINEX_HOME = $ENV{RECOMBINEX_HOME};
my $sample_table = "Master_Sample_Table.txt"; # master sample table
my $genome_dir;
my $genotype_dir;
my $batch_id;
my $recombination_profiling_dir;
my $color_scheme;
my $plot_centromere = "no";
my $flanking = 5000; # recombination event flanking region (bp) for plotting
my $basecall_qual_diff_cutoff = 30; # net quality difference cutoff used for genotyping

GetOptions('sample_table|s:s' => \$sample_table, # master sample table
	   'genome_dir|genome_dir:s' => \$genome_dir, # the path to the 01.Parental_Genome_Preprocessing directory
	   'genotype_dir|genotype_dir:s' => \$genotype_dir, # the path to the genotyping results in the 06.Tetrad_Genotyping directory
	   'batch_id|b:s' => \$batch_id, # the batch name used for genotyping and recombination profiling.
	   'recombination_profiling_dir|recombination_profiling_dir:s' => \$recombination_profiling_dir, # the path to the recombination profiling results in the 07.Recombination_Profiling directory
	   'qual_diff|q:i' => \$basecall_qual_diff_cutoff, # net quality difference used in tetrad genotyping, default = 30
	   'color_scheme|c:s' => \$color_scheme, # color_scheme for genotype plotting
	   'plot_centromere|plot_centromere:s' => \$plot_centromere, # whether to plot centromere
	   'flanking|f:i' => \$flanking); # recombination event flanking region (bp) for plotting

my $sample_table_fh = read_file($sample_table);
my %tetrads = parse_sample_table_by_tetrad($sample_table_fh);
# filter_tetrad_with_inviable_spores(\%tetrads);

foreach my $tetrad_id (sort keys %tetrads) {
    my $cross_pair = $tetrads{$tetrad_id}{'cross_pair'};
    my ($parent1, $parent2) = split "-", $cross_pair;
    print "tetrad: $tetrad_id, cross_pair: $cross_pair, parent1: $parent1, parent2: $parent2\n";
    my @ref = ("ref");
    foreach my $ref (@ref) {
	my @mode = ("raw", "inferred");
	foreach my $mode (@mode) {
	    my $genotype_for_plotting = "$genotype_dir/$batch_id/$cross_pair.$tetrad_id.$ref.q${basecall_qual_diff_cutoff}.genotype.lite.$mode.for_genotype_plotting.txt.gz";
	    my $recombination_events = "$recombination_profiling_dir/$batch_id/$tetrad_id/$cross_pair.$tetrad_id.$ref.q${basecall_qual_diff_cutoff}.genotype.lite.$mode.recombination_profile.events.txt";
	    if (not -e $genotype_for_plotting) {
		die "cannot find the genotype file $genotype_for_plotting!\n";
	    }
	    if (not -e $recombination_events) {
		die "cannot find the recombination events file $recombination_events!\n";
	    }
	    system("mkdir -p $recombination_profiling_dir/$batch_id/$tetrad_id/event_plots");
	    my $recombination_events_fh = read_file($recombination_events);
	    while (<$recombination_events_fh>) {
		chomp;
		/^#/ and next;
		/^\s*$/ and next;
		/^tetrad_id\tevent_id/ and next;
		my ($tetrad_id, $event_id, $event_type, $event_subtype, $chr, $marker_pos_start, $marker_pos_end, $adjusted_pos_start, $adjusted_pos_end, $adjusted_size, $marker_raw_index_start, $marker_raw_index_end, $marker_effective_index_start, $marker_effective_index_end, $affected_spores) = split /\t/, $_;
		if ($plot_centromere eq "yes") {
		    system("Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_tetrad_genotype.R --input $genotype_for_plotting --genome1_tag $parent1 --genome2_tag $parent2 --coordinate_genome_fai $genome_dir/$ref.genome.raw.relabel.fa.fai --coordinate_genome_centromere_gff $genome_dir/$ref.centromere.relabel.gff --color_scheme $color_scheme --output $recombination_profiling_dir/$batch_id/$tetrad_id/event_plots/$cross_pair.$tetrad_id.$ref.q${basecall_qual_diff_cutoff}.genotype.lite.$mode.recombination_event_plot.event_${event_id}.pdf --query_chr $chr --query_start $marker_pos_start --query_end $marker_pos_end --query_flanking $flanking");
		} else {
		    system("Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_tetrad_genotype.R --input $genotype_for_plotting --genome1_tag $parent1 --genome2_tag $parent2 --coordinate_genome_fai $genome_dir/$ref.genome.raw.relabel.fa.fai --color_scheme $color_scheme --output $recombination_profiling_dir/$batch_id/$tetrad_id/event_plots/$cross_pair.$tetrad_id.$ref.q${basecall_qual_diff_cutoff}.genotype.lite.$mode.recombination_event_plot.event_${event_id}.pdf --query_chr $chr --query_start $marker_pos_start --query_end $marker_pos_end --query_flanking $flanking");
		}
	    }
	}
    }
}

sub read_file {
    my $file = shift @_;
    my $fh;
    if ($file =~ /\.gz$/) {
	open($fh, "gunzip -c $file |") or die "can't open pipe to $file";
    } else {
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


sub parse_sample_table_by_tetrad {
    my $fh = shift @_;
    my %sample_table = ();
    while (<$fh>) {
        chomp;
        /^#/ and next;
	    /^\s*$/ and next;
        my ($sample_id, $tetrad_id, $spore_id, $PE_reads_files, $cross_pair, $note) = split /\s+/, $_;
        if (exists $sample_table{$tetrad_id}{'spore_index'}) {
            push @{$sample_table{$tetrad_id}{'spore_index'}}, $spore_id;
        } else {
            @{$sample_table{$tetrad_id}{'spore_index'}} = ($spore_id);
        }
        $sample_table{$tetrad_id}{'cross_pair'} = $cross_pair;
        $sample_table{$tetrad_id}{'note'} = $note;
    }
    return %sample_table;
}

	
sub filter_tetrad_with_inviable_spores {
    my $tetrads_hashref = shift @_;
    foreach my $tetrad_id (sort keys %$tetrads_hashref) {
	my $viable_spore_count = scalar @{$$tetrads_hashref{$tetrad_id}{'spore_index'}};
	if ($viable_spore_count < 4) {
	    delete $$tetrads_hashref{$tetrad_id};
	}
    }
}
