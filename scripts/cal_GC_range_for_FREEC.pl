#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Statistics::Descriptive;

my $input = "FREEC.refseq.GC_content.txt";
my $output = "FREEC.refseq.GC_range.txt";
my $window_size = 250;
my $lower_quantile = 15;
my $upper_quantile = 100 - $lower_quantile;
my $min_mappability = 0.85;
GetOptions('input|i:s' => \$input,
	   'output|o:s' => \$output,
	   'window_size|w:s' => \$window_size,
	   'lower_quantile|lower_quantile:s' => \$lower_quantile,
	   'upper_quantile|upper_quantile:s' => \$upper_quantile,
	   'min_mappability|min_mappability:s' => \$min_mappability
    );

my $GC_fh = read_file($input);
my %GC = parse_GC_file($GC_fh, $window_size);
my $output_fh = write_file($output);

my @GC_genome = ();    
foreach my $chr (sort keys %GC) {
    my @GC_chr = ();
    foreach my $i (sort {$a <=> $b} keys %{$GC{$chr}}) {
	if ($GC{$chr}{$i}{'effective_ratio'} > $min_mappability) {
	    push @GC_chr, $GC{$chr}{$i}{'GC'};
	}
    }
    if ((scalar @GC_chr) >= 10) {
	my $gc_stat_chr = Statistics::Descriptive::Full->new();
	$gc_stat_chr->add_data(@GC_chr);
	$gc_stat_chr->sort_data();
	$gc_stat_chr->presorted(1);
	# my $gc_mean_chr = sprintf("%.3f", $gc_stat_chr->mean());
	# my $gc_median_chr = sprintf("%.3f", $gc_stat_chr->median());
	my $gc_min_chr = sprintf("%.3f", $gc_stat_chr->min());
	my $gc_max_chr = sprintf("%.3f", $gc_stat_chr->max());
	# my $gc_stdev_chr = sprintf("%.3f", $gc_stat_chr->standard_deviation());
	my $gc_quantile_lower_chr = $gc_stat_chr->percentile($lower_quantile);
	$gc_quantile_lower_chr = sprintf("%.3f", $gc_quantile_lower_chr);
	my $gc_quantile_upper_chr = $gc_stat_chr->percentile($upper_quantile);
	$gc_quantile_upper_chr = sprintf("%.3f", $gc_quantile_upper_chr);
	push @GC_genome, $gc_quantile_upper_chr;
	push @GC_genome, $gc_quantile_lower_chr;
	# if ($gc_quantile_lower_chr < $min_expected_gc) {
	#     $min_expected_gc = $gc_quantile_lower_chr;
	# }
	# if ($gc_quantile_upper_chr > $max_expected_gc) {
	#     $max_expected_gc = $gc_quantile_upper_chr;
	# }
	print "chr=$chr\n";
	print "gc_quantile_upper_chr=$gc_quantile_upper_chr, gc_quantile_lower_chr=$gc_quantile_lower_chr\n";
	print "\n";
    }
}

my $gc_stat_genome = Statistics::Descriptive::Full->new();
$gc_stat_genome->add_data(@GC_genome);
$gc_stat_genome->sort_data();
$gc_stat_genome->presorted(1);
my $gc_min_genome = sprintf("%.3f", $gc_stat_genome->min());
my $gc_max_genome = sprintf("%.3f", $gc_stat_genome->max());

print "min_expected_gc=$gc_min_genome, max_expected_gc=$gc_max_genome\n";
print $output_fh "#min_expected_gc\tmax_expected_gc\n";
print $output_fh "$gc_min_genome\t$gc_max_genome\n";


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

sub parse_GC_file {
    my ($fh, $window_size) = @_;
    my %GC= ();
    my $i = 0;
    while (<$fh>) {
        chomp;
        /^#/ and next;
	    /^\s*$/ and next;
        my ($chr, $window_start, $window_end, $AT_pct, $GC_pct, $A_count, $C_count, $G_count, $T_count, $N_count, $other_count, $window_seq_length) = split /\s+/, $_;
	if (not exists $GC{$chr}) {
	    $i = 0;
	} else {
	    $i++;
	}

	my $AT_count = $A_count + $T_count;
	my $GC_count = $G_count + $C_count;
	$GC{$chr}{$i}{'start'} = $window_start;
	$GC{$chr}{$i}{'end'} = $window_end;
	$GC{$chr}{$i}{'effective_window_size'} = $AT_count + $GC_count;
	$GC{$chr}{$i}{'effective_ratio'} = sprintf("%.3f", $GC{$chr}{$i}{'effective_window_size'}/$window_size);
	if ($GC{$chr}{$i}{'effective_window_size'} > 0) {
	    $GC{$chr}{$i}{'GC'} = $GC_count/$GC{$chr}{$i}{'effective_window_size'};
	} else {
	    $GC{$chr}{$i}{'GC'} = "NA";
	}
    }
    return %GC;
}
