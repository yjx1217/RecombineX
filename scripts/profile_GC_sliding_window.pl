#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Statistics::Descriptive;

##############################################################
#  script: profile_GC_sliding_window.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.10.10
#  description: generate GC-cotent profile for the input multi-fasta file.
#  example: perl profile_GC_sliding_window.pl  -i input.fa(.gz) -p prefix -w window_size -s step_size -m minimal_effective_ratio -q upper_quantile
#  example: perl profile_GC_sliding_window.pl  -i input.fa(.gz) -p prefix -w 250 -s 250 -m 0.5 -q 85
##############################################################

my ($input, $prefix, $window_size, $step_size, $upper_quantile, $min_effective_ratio);

$window_size = 250;
$step_size = $window_size;
$upper_quantile = 85;
$min_effective_ratio = 0.5;

GetOptions('input|i:s' => \$input,
           'prefix|p:s' => \$prefix,
	   'window|w:i' => \$window_size,
	   'step|s:i' => \$step_size,
	   'quantile|q:i' => \$upper_quantile, # e.g. 90
	   'min_effective_ratio|m:f' => \$min_effective_ratio); 

my $lower_quantile = 100 - $upper_quantile;


my $input_fh = read_file($input);
my %input = ();
my @input = ();

parse_fasta_file($input_fh, \%input, \@input);

my $output_detail = "$prefix.GC_profile.detail.txt";
my $output_detail_fh = write_file($output_detail);
print $output_detail_fh "#window=$window_size, step=$step_size\n";
print $output_detail_fh "#chromosome\twindow_index\tGC_content\teffective_window_size\teffective_window\teffective_ratio\n";

my $output_summary = "$prefix.GC_profile.summary.txt";
my $output_summary_fh = write_file($output_summary);
print $output_summary_fh "#window=$window_size, step=$step_size, min_effective_ratio=$min_effective_ratio\n";
print $output_summary_fh "#chromosome\tmean_GC\tmedian_GC\tstdev_GC\t${lower_quantile}%_quantile_GC\t${upper_quantile}%_quantile_GC\tmin_GC\tmax_GC\n";

my $output_lite = "$prefix.GC_profile.lite.W${window_size}S${step_size}.txt";
my $output_lite_fh = write_file($output_lite);

# define scanning window
# for a given seq with a length of l, the total number of windows with a size of w: 
# = [l - w]/s + 1

my %GC = ();
foreach my $chr (@input) {
    my $remainder = $input{$chr};
    my $i = 0;
    my $GC_stat_chr = Statistics::Descriptive::Full->new();
    my @GC_chr = ();
    while ((length $remainder) >= $window_size) {
	$i++;
	my $window_start = $step_size * ($i - 1) + 1;
	my $window_end = $window_start + $window_size - 1;
	my $subseq = substr $remainder, 0, $window_size;
	$subseq = uc $subseq;
	$GC{$chr}{$i}{'start'} = $window_start;
	$GC{$chr}{$i}{'end'} = $window_end;
	my $GC_count  = () = $subseq =~ /(G|C)/g;
	my $AT_count  = () = $subseq =~ /(A|T)/g;
	$GC{$chr}{$i}{'effective_window_size'} = $AT_count + $GC_count;
	$GC{$chr}{$i}{'effective_ratio'} = sprintf("%.3f", $GC{$chr}{$i}{'effective_window_size'}/$window_size);
	my $GC_window;
	if ($GC{$chr}{$i}{'effective_window_size'} > 0) {
	    $GC{$chr}{$i}{'GC'} = $GC_count/$GC{$chr}{$i}{'effective_window_size'};
	    $GC_window = sprintf("%.3f", $GC{$chr}{$i}{'GC'});
	} else {
	    $GC{$chr}{$i}{'GC'} = "NA";
	    $GC_window = "NA";
	}
	print $output_detail_fh "$chr\t$i\t$GC_window\t$GC{$chr}{$i}{'effective_window_size'}\t$GC{$chr}{$i}{'effective_ratio'}\n";
	my $window_midpoint = ($window_start + $window_end)/2;
	print $output_lite_fh "$chr\t$window_midpoint\t$GC_window\n";
	if (($AT_count + $GC_count)/$window_size >= $min_effective_ratio) {
	    push @GC_chr, $GC{$chr}{$i}{'GC'};
	}
	$remainder = substr $remainder, $step_size;
    }
    $GC_stat_chr->add_data(@GC_chr);
    $GC_stat_chr->sort_data();
    $GC_stat_chr->presorted(1);
    my $GC_mean_chr = sprintf("%.3f", $GC_stat_chr->mean());
    my $GC_median_chr = sprintf("%.3f", $GC_stat_chr->median());
    my $GC_min_chr = sprintf("%.3f", $GC_stat_chr->min());
    my $GC_max_chr = sprintf("%.3f", $GC_stat_chr->max());
    my $GC_stdev_chr = sprintf("%.3f", $GC_stat_chr->standard_deviation());
    my $GC_quantile_lower_chr = $GC_stat_chr->percentile($lower_quantile);
    $GC_quantile_lower_chr = sprintf("%.3f", $GC_quantile_lower_chr);
    my $GC_quantile_upper_chr = $GC_stat_chr->percentile($upper_quantile);
    $GC_quantile_upper_chr = sprintf("%.3f", $GC_quantile_upper_chr);

    print $output_summary_fh "$chr\t$GC_mean_chr\t$GC_median_chr\t$GC_stdev_chr\t$GC_quantile_lower_chr\t$GC_quantile_upper_chr\t$GC_min_chr\t$GC_max_chr\n";
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
    if ($file =~ /.gz$/) {
        open($fh, "| gzip -c >$file") or die "can't open $file\n";
    } else {
        open($fh, ">$file") or die "can't open $file\n";
    }
    return $fh;
}  

sub parse_fasta_file {
    my ($fh, $input_hashref, $input_arrayref) = @_;
    my $seq_name = "";
    while (<$fh>) {
        chomp;
        if (/^\s*$/) {
            next;
        } elsif (/^\s*#/) {
            next;
        } elsif (/^>(.*)/) {
            $seq_name = $1;
            push @$input_arrayref, $seq_name;
            $$input_hashref{$seq_name} = "";
        } else {
            $$input_hashref{$seq_name} .= $_;
        }
    }
}


