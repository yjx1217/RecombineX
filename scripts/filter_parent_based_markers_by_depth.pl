#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: filter_parent_based_markers_by_depth.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2021.03.10
#  description: filter out markers with unexpected depth (<0.5 genome-wide median depth or >1.5 genome-wide median depth)
#  example: perl filter_parent_based_markers_by_depth.pl -depth_detail depth.txt.gz -depth_summary depth.summary.txt -i markers.vcf.gz -o filtered.markers.vcf.gz
##############################################################

my ($input, $output, $depth_summary, $depth_detail);

GetOptions('i|input:s' => \$input,
	   'depth_summary|depth_summary:s' => \$depth_summary,
	   'depth_detail|depth_detail:s' => \$depth_detail,
	   'o|output:s' => \$output);

my $depth_summary_fh = read_file($depth_summary);
my %depth_summary = parse_depth_summary_file($depth_summary_fh);

my $depth_detail_fh = read_file($depth_detail);
my %depth_detail = ();
parse_depth_detail_file($depth_detail_fh, \%depth_detail);

my $output_fh = write_file($output);

my $filtered_count = 0;

my $input_fh = read_file($input);
while (<$input_fh>) {
    chomp;
    if (/^#/) {
	print $output_fh "$_\n";
	next;
    } elsif (/^\s*$/) {
	print $output_fh "$_\n";
	next;
    } elsif (/^chr\tstart/) {
	print $output_fh "$_\n";
	next;
    } elsif (/^ref_chr\tref_start/) {
	print $output_fh "$_\n";
	next;
    } elsif (/^parent1_chr\tparent1_start/) {
	print $output_fh "$_\n";
	next;
    } elsif (/^parent2_chr\tparent2_start/) {
	print $output_fh "$_\n";
	next;
    }
    my @line = split /\t/, $_;
    my $ref_chr = $line[0];
    my $ref_start = $line[1];
    my ($ref_end) = ($_ =~ /ref_end=(\d+)/);
    my $flag = 1;
    if (($depth_detail{$ref_chr}{$ref_start} > $depth_summary{$ref_chr} * 0.5) and ($depth_detail{$ref_chr}{$ref_start} < $depth_summary{$ref_chr} * 1.5)) {
	if (($depth_detail{$ref_chr}{$ref_end} > $depth_summary{$ref_chr} * 0.5) and ($depth_detail{$ref_chr}{$ref_end} < $depth_summary{$ref_chr} * 1.5)) {
	    $flag = 0;
	}
    }
    if ($flag == 0) {
	print $output_fh "$_\n";
    } else {
	$filtered_count++;
    }
}

print "filtered $filtered_count markers due to substantial depth bias!\n";

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

sub parse_depth_summary_file {
    my $fh = shift @_;
    my @header = ();
    my @content = ();
    my %median_depth_by_chr = ();
    while (<$fh>) {
        chomp;
        /^#/ and next;
	/^\s*$/ and next;
	if (/sample\tref_genome\ttotal_reads/) {
	    @header = split /\t/, $_;
	} else {
	    @content = split /\t/, $_;
	}
    }
    my $column_num = scalar @header;
    for (my $i = 10; $i < $column_num; $i += 2) {
	my $chr = $header[$i];
	$chr =~ s/_median_depth//gi;
	$median_depth_by_chr{$chr} = $content[$i];
	# print "chr=$chr, median_depth=$content[$i]\n";
    }
    return %median_depth_by_chr;
}


sub parse_depth_detail_file {
    my ($fh, $depth_hashref) = @_;
    while (<$fh>) {
        chomp;
        /^#/ and next;
	/^\s*$/ and next;
        my ($chr, $pos, $depth) = split /\t/, $_;
	$$depth_hashref{$chr}{$pos} = $depth;
    }
}
