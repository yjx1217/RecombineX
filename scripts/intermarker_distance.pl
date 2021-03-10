#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: intermarker_distance.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.02.28
#  description: calculate distances between consecutive markers
#  example: perl intermarker_distance.pl -i input.SNP.markers.txt(.gz) -p $prefix
##############################################################


my ($input, $prefix);

GetOptions('i|input:s' => \$input,
	   'p|prefix:s' => \$prefix);

my $input_fh = read_file($input);
my %markers = ();
my $marker_count = 0;
$marker_count = parse_markers($input_fh, \%markers);

my $offset;
my $output = "$prefix.intermarker_distance.txt";
my $output_fh = write_file($output);

my @distance = ();
foreach my $chr (sort keys %{$markers{'ref'}}) {
    $offset = "";
    foreach my $midpoint (sort {$a <=> $b} keys %{$markers{'ref'}{$chr}}) {
	if ($offset eq "") {
	    $offset = $midpoint;
	} else {
	    my $d = $midpoint - $offset;
	    push @distance, $d;
	    $offset = $midpoint;
	}
    }
}

my $mean_distance = cal_mean(@distance);
my $median_distance = cal_median(@distance);
my $stdev_distance = cal_stdev(@distance);
$mean_distance = sprintf("%.2f", $mean_distance);
$median_distance = sprintf("%.2f", $median_distance);
$stdev_distance = sprintf("%.2f", $stdev_distance);

print $output_fh "prefix\tmarker_count\tmean_distance\tmedian_distance\tstdev_distance\n";
print $output_fh "$prefix\t$marker_count\t$mean_distance\t$median_distance\t$stdev_distance\n";


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


sub parse_markers {
    my ($fh, $marker_hashref) = @_;
    my $marker_count = 0;
    while (<$fh>) {
	chomp;
	/^chr\tstart\tend/ and next;
	/^ref_chr\tref_start\tref_end/ and next;
	/^parent1_chr\tparent1_start\tparent1_end/ and next;
	/^parent2_chr\tparent2_start\tparent2_end/ and next;
	my ($ref_chr, $ref_start, $ref_end, $ref_allele, $query_allele, $query_chr, $query_start, $query_end, $query_orientation) = split /\t/, $_;
	my $ref_midpoint = ($ref_start + $ref_end)/2;
	my $query_midpoint = ($query_start + $query_end)/2;
	$$marker_hashref{'ref'}{$ref_chr}{$ref_midpoint} = $_;
	$$marker_hashref{'query'}{$query_chr}{$query_midpoint} = $_;
	$marker_count++;
    }
    return $marker_count;
}

sub cal_mean {
    my @data = @_;
    my $data_length = scalar @data;
    my $mean;
    if ($data_length == 0) {
	$mean = "NA";
    } else {
	my $sum = 0;
	foreach my $element (@data) {
	    $sum += $element;
	}
	$mean = $sum/$data_length;
    }
    return $mean;
}

sub cal_median {
    my @data = @_;
    my $data_length = scalar @data;
    my $median;
    if ($data_length == 0) {
	$median = "NA";
    } else {
	my @data_sorted = sort {$a <=> $b} @data;
	if ($data_length % 2 == 1) {
	    $median = $data_sorted[int($data_length/2)];
	} else {
	    my ($upper, $lower);
	    $lower = $data_sorted[int($data_length/2)];
	    $upper = $data_sorted[int($data_length/2)-1];
	    $median = ($lower + $upper)/2;
	}
    }
    return $median;
}

sub cal_stdev {
    my @data = @_;
    my $data_length = scalar @data;
    my $stdev;
    if ($data_length == 0) {
	$stdev = "NA";
    } else {
	my $mean = cal_mean(@data);
	my $sqdev_sum = 0;
	foreach my $element (@data) {
	    $sqdev_sum += ($element - $mean)**2;
	}
	$stdev = sqrt($sqdev_sum/($data_length-1));
    }
    return $stdev;
}

