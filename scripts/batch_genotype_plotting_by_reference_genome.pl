#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Env;

##############################################################
#  script: batch_genotype_plotting_by_reference_genome.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.09.19
#  description: plot genotype for the tetrad
#  example: perl batch_genotype_plotting_by_reference_genome.pl -s Master_Sample_Table.txt -m raw -genome_dir genome_dir -q 30 -genotype_dir genotype_dir -batch_id $batch_id -output_dir output_dir -color_scheme color_scheme -plot_centromere yes
##############################################################

my $RECOMBINEX_HOME = $ENV{RECOMBINEX_HOME};
my $sample_table = "Master_Sample_Table.txt"; # master sample table
my $mode = "raw"; # "raw" or "inferred"
my $basecall_qual_diff_cutoff = 30; # net quality difference cutoff used for genotyping
my $genome_dir;
my $genotype_dir;
my $batch_id;
my $color_scheme;
my $plot_centromere = "no";
my $output_dir;

GetOptions('sample_table|s:s' => \$sample_table,
	   'mode|m:s' => \$mode,
	   'qual|q:f' => \$basecall_qual_diff_cutoff,
	   'genome_dir|genome_dir:s' => \$genome_dir,
	   'genotype_dir|genotype_dir:s' => \$genotype_dir,
	   'batch_id|b:s' => \$batch_id,
	   'output_dir|output_dir:s' => \$output_dir,
	   'color_scheme|c:s' => \$color_scheme,
	   'plot_centromere|plot_centromere:s' => \$plot_centromere, # whether to plot centromere
    );

my $sample_table_fh = read_file($sample_table);
my %tetrads = parse_sample_table_by_tetrad($sample_table_fh);
close $sample_table_fh;

foreach my $tetrad_id (sort keys %tetrads) {
    my $cross_pair = $tetrads{$tetrad_id}{'cross_pair'};
    my ($genome1_tag, $genome2_tag) = split "-", $cross_pair;
    print "tetrad: $tetrad_id, cross_pair: $cross_pair, parent1: $genome1_tag, parent2: $genome2_tag\n";
    # my @refseq_tags = ($genome1_tag, $genome2_tag);
    my @refseq_tags = ("ref");
    foreach my $refseq_tag (@refseq_tags) {
	my $prefix = "$cross_pair.$tetrad_id.$refseq_tag.q${basecall_qual_diff_cutoff}";
	my $genotype_input = "$genotype_dir/$batch_id/$prefix.genotype.lite.$mode.txt.gz";
	my $genotype_input_fh = read_file($genotype_input);
	my %genotypes = parse_genotype_file($genotype_input_fh);
	close $genotype_input_fh;
	my $output = "$output_dir/$prefix.genotype.lite.$mode.for_genotype_plotting.txt.gz";
	my $output_fh = write_file($output);
	print $output_fh "chr\tmarker_id\traw_marker_start\traw_marker_end\tadjusted_marker_start\tadjusted_marker_end\tspore_genotype\tspore_id\n";
	my @spores = qw(a b c d);
	foreach my $chr (sort keys %genotypes) {
	    my @pos = sort {$a <=> $b} keys %{$genotypes{$chr}};
	    my $first_marker_pos = $pos[0];
	    my $last_marker_pos = $pos[-1];
	    foreach my $spore_id (@spores) {
		my $marker_id = 1;
		for (my $i = 0; $i < scalar @pos; $i++) {
		    my $raw_marker_start = $pos[$i];
		    my $raw_marker_end = $pos[$i];
		    my $adjusted_marker_start;
		    my $adjusted_marker_end;
		    my $marker_genotype;
		    if ($pos[$i] == $first_marker_pos) {
		    	$adjusted_marker_start = $pos[$i];
		    	$adjusted_marker_end = ($pos[$i] + $pos[$i+1])/2;
		    } elsif ($pos[$i] == $last_marker_pos) {
		    	$adjusted_marker_start = ($pos[$i-1] + $pos[$i])/2;
		    	$adjusted_marker_end = $pos[$i];
		    } else {
		    	$adjusted_marker_start = ($pos[$i-1] + $pos[$i])/2;
		    	$adjusted_marker_end = ($pos[$i] + $pos[$i+1])/2;
		    }
		    $marker_genotype = $genotypes{$chr}{$pos[$i]}{$spore_id};
		    print $output_fh "$chr\t$marker_id\t$raw_marker_start\t$raw_marker_end\t$adjusted_marker_start\t$adjusted_marker_end\t$marker_genotype\t$spore_id\n";
		    $marker_id++;
		}
	    }
	}
	close $output_fh;
	my $genotype_for_plotting="$output_dir/$prefix.genotype.lite.$mode.for_genotype_plotting.txt.gz";
	if ($plot_centromere eq "yes") {
	    system("Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_tetrad_genotype.R --input $genotype_for_plotting --genome1_tag $genome1_tag --genome2_tag $genome2_tag --coordinate_genome_fai $genome_dir/$refseq_tag.genome.raw.relabel.fa.fai --coordinate_genome_centromere_gff  $genome_dir/$refseq_tag.centromere.relabel.gff --color_scheme $color_scheme --output $output_dir/$cross_pair.$tetrad_id.$refseq_tag.q${basecall_qual_diff_cutoff}.genotype.lite.$mode.genotype_plot.pdf");
	} else {
	    system("Rscript --vanilla --slave $RECOMBINEX_HOME/scripts/plot_tetrad_genotype.R --input $genotype_for_plotting --genome1_tag $genome1_tag --genome2_tag $genome2_tag --coordinate_genome_fai $genome_dir/$refseq_tag.genome.raw.relabel.fa.fai --color_scheme $color_scheme --output $output_dir/$cross_pair.$tetrad_id.$refseq_tag.q${basecall_qual_diff_cutoff}.genotype.lite.$mode.genotype_plot.pdf");
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

sub parse_genotype_file {
    my $fh = shift @_;
    my %genotypes = ();
    while (<$fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	my ($chr, $pos, $refseq_tag, $genotype_a, $genotype_b, $genotype_c, $genotype_d) = split /\t/, $_;
	my @genotypes = ($genotype_a, $genotype_b, $genotype_c, $genotype_d);
	# print "chr=$chr, pos=$pos, refseq_tag=$reftag, genotypes = @genotypes\n";
	my @spores = qw(a b c d);
	foreach my $spore_id (@spores) {
	    my $gt = shift @genotypes;
	    $genotypes{$chr}{$pos}{$spore_id} = $gt;
	}
    }
    return %genotypes;
}
