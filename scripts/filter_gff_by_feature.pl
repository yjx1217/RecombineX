#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: filter_gff_by_feature.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2019.01.28
#  description: filter gff by genomic feature type(s) (e.g. "CDS" or "TY,tRNA"), "-m keep" to retain specified features whereas "-m remove" to remove specified features.  
#  example: perl filter_gff_by_feature.pl -i input.gff(.gz) -o output.gff(.gz) -f centromere,tRNA -m keep 
##############################################################

my ($input, $output, $feature, $mode);
$feature = "all"; #
$mode = "keep"; # or "remove"

GetOptions('input|i:s' => \$input,
	   'output|o:s' => \$output,
	   'feature|f:s' => \$feature, # a space-separated list of features. e.g. "TY,tRNA"
	   'mode|m:s' => \$mode); # mode: keep, remove

my $input_fh = read_file($input);
# print "feature=$feature\n";
my @feature;
my $matching_pattern;
if ($feature eq "all") {
    $matching_pattern = ".*";
} else {
    @feature = split ",", $feature; 
    $matching_pattern = join "|", @feature;
    $matching_pattern = "($matching_pattern)";
}


# print "matching_pattern=$matching_pattern\n"; 

my $output_fh = write_file($output);

while (<$input_fh>) {
    chomp;
    /^\##FASTA/ and last;
    /^\#/ and next;
    /^\s*$/ and next;
    my @line = split /\s+/, $_;
    my $type = $line[2];
    if ($type =~ /\b$matching_pattern\b/) {
	if ($mode eq "keep") {
	    print $output_fh "$_\n";
	}
    } else {
	if ($mode eq "remove") {
	    print $output_fh "$_\n";
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

