#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
#use Statistics::Descriptive;

##############################################################
#  script: profile_marker_density_sliding_window.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.10.10
#  description: generate marker density profile for the input multi-fasta file.
#  example: perl profile_marker_density_sliding_window.pl -r $ref.genome.raw.relabel.fa -i ./../14.Polymorphic_Markers_by_Consensus/$cross_pair.$ref.final.SNP.markers.txt.gz -p $prefix -w $window_size -s $step_size
##############################################################

my ($input, $refseq, $prefix, $window_size, $step_size);

$window_size = 250;
$step_size = $window_size;

GetOptions('input|i:s' => \$input,
           'refseq|r:s' => \$refseq,
           'prefix|p:s' => \$prefix,
	   'window|w:i' => \$window_size,
	   'step|s:i' => \$step_size);

my $input_fh = read_file($input);

my $refseq_fh = read_file($refseq);
my %refseq = ();
my @refseq = ();
parse_fasta_file($refseq_fh, \%refseq, \@refseq);

my $output_lite = "$prefix.marker_density_profile.lite.W${window_size}S${step_size}.txt";
my $output_lite_fh = write_file($output_lite);

# define scanning window
# for a given seq with a length of l, the total number of windows with a size of w: 
# = [l - w]/s + 1

my %window = ();
foreach my $chr (@refseq) {
    my $remainder = $refseq{$chr};
    my $i = 0;
    while ((length $remainder) >= $window_size) {
	$i++;
	my $window_start = $step_size * ($i - 1) + 1;
	my $window_end = $window_start + $window_size - 1;
	$window{$chr}{$i}{'window_start'} = $window_start;
	$window{$chr}{$i}{'window_end'} = $window_end;
	$window{$chr}{$i}{'window_midpoint'} = ($window_start + $window_end)/2;
	$window{$chr}{$i}{'marker_count'} = 0;
	$remainder = substr $remainder, $step_size;
    }
}

my %markers = ();
while (<$input_fh>) {
    chomp;
    /^\s*$/ and next;
    /^#/ and next;
    /^parent\w_chr\tparent\w_start/ and next;
    my @line = split /\t/, $_;
    my $chr = $line[0];
    my $start = $line[1];
    my $end = $line[2];
    foreach my $i (sort {$a <=> $b} keys %{$window{$chr}}) {
	my $window_start = $window{$chr}{$i}{'window_start'};
	my $window_end = $window{$chr}{$i}{'window_end'};
	if (($start >= $window_start) and ($end <= $window_end)) {
	    $window{$chr}{$i}{'marker_count'}++;
	    last;
	}
    }
}

foreach my $chr (@refseq) {
    foreach  my $i (sort {$a <=> $b} keys %{$window{$chr}}) {
	my $window_midpoint = $window{$chr}{$i}{'window_midpoint'};
	my $window_marker_density = $window{$chr}{$i}{'marker_count'}/$window_size;
	print $output_lite_fh "$chr\t$window_midpoint\t$window_marker_density\n";
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


