#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: tidy_SGDref_genome.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.02.09
#  description: simplify the chromosome tags in the SGD genome release multi-fasta file 
#  example: perl tidy_SGDref_genome.pl  -i S288C_reference_sequence_R64-2-1_20150113.fsa(.gz) -o SGDref.fa(.gz) 
##############################################################

my ($input, $output);

GetOptions('input|i:s' => \$input,
	   'output|o:s' => \$output);

my $input_fh = read_file($input);
my @input = ();
my %input = ();
parse_fasta_file($input_fh, \%input, \@input);

my $output_fh = write_file($output);
    
foreach my $id (@input) {
    my $new_id;
    if ($id =~ /\[(?:chromosome|location)=([^\]]+)\]/) {
	$new_id = "chr". $1;
	print $output_fh ">$new_id\n$input{$id}\n";
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


