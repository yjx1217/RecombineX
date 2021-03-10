#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Statistics::Descriptive;

##############################################################
#  script: profile_masking_sliding_window.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.10.10
#  description: generate masked-content profile for the input multi-fasta file.
#  example: perl profile_MaskedRegion_sliding_window.pl  -i input.fa(.gz) -p prefix -w window_size -s step_size 
#  example: perl profile_MaskedRegion_sliding_window.pl  -i input.fa(.gz) -p prefix -w 250 -s 250 
##############################################################

my ($input, $prefix, $window_size, $step_size);

$window_size = 250;
$step_size = $window_size;

GetOptions('input|i:s' => \$input,
           'prefix|p:s' => \$prefix,
	   'window|w:i' => \$window_size,
	   'step|s:i' => \$step_size);

my $input_fh = read_file($input);
my %input = ();
my @input = ();
parse_fasta_file($input_fh, \%input, \@input);

my $output_lite = "$prefix.masking_profile.lite.W${window_size}S${step_size}.txt";
my $output_lite_fh = write_file($output_lite);

# define scanning window
# for a given seq with a length of l, the total number of windows with a size of w: 
# = [l - w]/s + 1

my %masking = ();
foreach my $chr (@input) {
    my $remainder = $input{$chr};
    my $i = 0;
    while ((length $remainder) >= $window_size) {
	$i++;
	my $window_start = $step_size * ($i - 1) + 1;
	my $window_end = $window_start + $window_size - 1;
	my $subseq = substr $remainder, 0, $window_size;
	$subseq = uc $subseq;
	$masking{$chr}{$i}{'start'} = $window_start;
	$masking{$chr}{$i}{'end'} = $window_end;
	my $GC_count  = () = $subseq =~ /(G|C)/g;
	my $AT_count  = () = $subseq =~ /(A|T)/g;
	$masking{$chr}{$i}{'masking_size'} = $window_size - ($AT_count + $GC_count);
	$masking{$chr}{$i}{'masking_ratio'} = sprintf("%.3f", $masking{$chr}{$i}{'masking_size'}/$window_size);
	my $masking_window = $masking{$chr}{$i}{'masking_ratio'};
	my $window_midpoint = ($window_start + $window_end)/2;
	print $output_lite_fh "$chr\t$window_midpoint\t$masking_window\n";
	$remainder = substr $remainder, $step_size;
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


