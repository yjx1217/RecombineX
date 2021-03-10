#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: relabel_chr_in_gff.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2017.09.11
#  description: attach a prefix for the chr/contig name in the gff file
#  example: perl relabel_chr_in_gff.pl -i input.gff -t prefix_tag -o output.relabeled.gff
##############################################################

my ($input, $tag, $output);

GetOptions('input|i:s' => \$input,
	   'tag|t:s' => \$tag,
	   'output|o:s' => \$output);

my $input_fh = read_file($input);
my $output_fh = write_file($output);

while (<$input_fh>) {
    chomp;
    /^\##FASTA/ and last;
    /^\#/ and next;
    /^\s*$/ and next;
    my @line = split /\s+/, $_;
    $line[0] = $tag . "_". $line[0];
    my $line = join "\t", @line;
    print $output_fh "$line\n";
}
close $input_fh;

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


