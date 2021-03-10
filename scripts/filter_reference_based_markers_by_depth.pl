#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: filter_markers_by_depth.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2021.03.10
#  description: filter out markers with unexpected depth (<0.5 genome-wide median depth or >1.5 genome-wide median depth)
#  example: perl filter_markers_by_depth.pl -ref_depth_detail ref.depth.txt.gz -ref_depth_summary ref.depth.summary.txt -query_depth_detail query.depth.txt.gz -query_depth_summary query.depth.summary.txt_ -i marker_table -o filtered_marker_table
##############################################################

my ($input, $output, $ref_depth_summary, $ref_depth_detail, $query_depth_summary, $query_depth_detail);

GetOptions('i|input:s' => \$input,
	   'ref_depth_summary|ref_depth_summary:s' => \$ref_depth_summary,
	   'query_depth_summary|query_depth_summary:s' => \$query_depth_summary,
	   'ref_depth_detail|ref_depth_detail:s' => \$ref_depth_detail,
	   'query_depth_detail|query_depth_detail:s' => \$query_depth_detail,
	   'o|output:s' => \$output);

my $ref_depth_summary_fh = read_file($ref_depth_summary);
my %ref_depth_summary = parse_depth_summary_file($ref_depth_summary_fh);

my $query_depth_summary_fh = read_file($query_depth_summary);
my %query_depth_summary = parse_depth_summary_file($query_depth_summary_fh); 

my $ref_depth_detail_fh = read_file($ref_depth_detail);
my %ref_depth_detail = ();
parse_depth_detail_file($ref_depth_detail_fh, \%ref_depth_detail);

my $query_depth_detail_fh = read_file($query_depth_detail);
my %query_depth_detail = ();
parse_depth_detail_file($query_depth_detail_fh, \%query_depth_detail);

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
    my ($ref_chr, $ref_start, $ref_end, $ref_allele, $query_allele, $query_chr, $query_start, $query_end, $match_orientation) = split /\t/, $_;
    my $flag_by_ref = 1;
    my $flag_by_query = 1;
    if (($ref_depth_detail{$ref_chr}{$ref_start} > $ref_depth_summary{$ref_chr} * 0.5) and ($ref_depth_detail{$ref_chr}{$ref_start} < $ref_depth_summary{$ref_chr} * 1.5)) {
	if (($ref_depth_detail{$ref_chr}{$ref_end} > $ref_depth_summary{$ref_chr} * 0.5) and ($ref_depth_detail{$ref_chr}{$ref_end} < $ref_depth_summary{$ref_chr} * 1.5)) {
	    $flag_by_ref = 0;
	}
    }
    if (($query_depth_detail{$query_chr}{$query_start} > $query_depth_summary{$query_chr} * 0.5) and ($query_depth_detail{$query_chr}{$query_start} < $query_depth_summary{$query_chr} * 1.5)) {
	if (($query_depth_detail{$query_chr}{$query_end} > $query_depth_summary{$query_chr} * 0.5) and ($query_depth_detail{$query_chr}{$query_end} < $query_depth_summary{$query_chr} * 1.5)) {
	    $flag_by_query = 0;
	}
    }
    if (($flag_by_ref == 0) and ($flag_by_query == 0)) {
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
