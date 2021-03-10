#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;


my ($input, $prefix, $excluded_chr_list);
$excluded_chr_list = "";
GetOptions('input|i:s' => \$input,
	   'prefix|p:s' => \$prefix,
	   'excluded_chr_list|e:s' => \$excluded_chr_list);


my $input_fh = read_file($input);
my %input = ();
my @input = ();
parse_fasta_file($input_fh, \%input, \@input);
    
my $output = "$prefix.FREEC.fa";
my $output_fh = write_file($output);

my @excluded_chr_list_regexp = "";
if ($excluded_chr_list ne "") {
    print "excluded_chr_list = $excluded_chr_list\n";
    my $excluded_chr_list_fh = read_file($excluded_chr_list);
    my %excluded_chr = parse_list_file($excluded_chr_list_fh);
    foreach my $chr (sort keys %excluded_chr) {
        push @excluded_chr_list_regexp, "/$chr/{d;b}";
    }

    my @input_filtered = ();
    foreach my $chr (@input) {
        if (not exists $excluded_chr{$chr}) {
            push @input_filtered, $chr;
        } else {
            delete $input{$chr};
        }
    }
    @input = @input_filtered;
}

my %FREEC_chr_dict = ();
my @FREEC_chr = ();
my $FREEC_chr_dict_output = "$prefix.FREEC_chr.dict";
my $FREEC_chr_dict_output_fh = write_file($FREEC_chr_dict_output);

# setup chr, chrLength   
system("mkdir ${prefix}_FREEC_chr");
my $chr_index = 0;
foreach my $chr (@input) {
    $chr_index++;
    my $FREEC_chr = "chr". $chr_index;
    $FREEC_chr_dict{'FREEC_to_original'}{$chr_index} = $chr;
    $FREEC_chr_dict{'original_to_FREEC'}{$chr} = $FREEC_chr;
    push @FREEC_chr, $FREEC_chr;
    print $output_fh ">$FREEC_chr\n$input{$chr}\n";
    my $FREEC_chr_fa = "./${prefix}_FREEC_chr/${FREEC_chr}.fa";
    my $FREEC_chr_fa_fh = write_file($FREEC_chr_fa);
    print $FREEC_chr_fa_fh ">$FREEC_chr\n$input{$chr}\n";
    close $FREEC_chr_fa_fh;
    print $FREEC_chr_dict_output_fh "original_to_FREEC\t$chr\t$FREEC_chr\n";
    print $FREEC_chr_dict_output_fh "FREEC_to_original\t$chr_index\t$chr\n";
}
close $output_fh;
close $FREEC_chr_dict_output_fh;

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

sub parse_list_file {
    my $fh = shift @_;
    my %list = ();
    while (<$fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	if (exists $list{$_}) {
	    $list{$_}++;
	} else {
	    $list{$_} = 1;
	}
    }
    return %list
}
