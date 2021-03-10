#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: filter_dubious_subtelomeric_makers.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.02.28
#  description: filter out markers involved in subtelomeric duplication blocks.
#  example: perl filter_dubious_subtelomeric_makers.pl -i raw.vcf -rt ref -qt query  -rcg ref_centromere_gff -qcg query_centromere_gff -o output.vcf
##############################################################

my ($input, $output, $ref_tag, $query_tag, $ref_centromere_gff, $query_centromere_gff);

GetOptions('i|input:s' => \$input,
	   'rt|ref_tag:s' => \$ref_tag,
	   'qt|query_tag:s' => \$query_tag,
	   'rcg|ref_centromere_gff:s' => \$ref_centromere_gff,
	   'qcg|query_centromere_gff:s' => \$query_centromere_gff,
	   'o|output:s' => \$output);

my $input_fh = read_file($input);

my $ref_centromere_gff_fh = read_file($ref_centromere_gff);
my %ref_centromere_gff = parse_gff_file($ref_centromere_gff_fh);

my $query_centromere_gff_fh = read_file($query_centromere_gff);
my %query_centromere_gff = parse_gff_file($query_centromere_gff_fh);

my $output_fh = write_file($output);
my $filtered_count = 0;

while (<$input_fh>) {
    chomp;
    /^\s*$/ and next;
    if (/^##source=/) {
	print $output_fh "##source=filter_dubious_subtelomeric_makers.pl\n";
    } elsif (/^#/) {
	print $output_fh "$_\n";
    } else {
	my ($chr, $ref_pos, $var_id, $ref_base, $alt_base, $qual, $filter, $info) = split /\t/, $_;
	# print "info=$info\n";
	my ($type, $ref_chr, $ref_start, $ref_end, $query_chr, $query_start, $query_end, $query_orientation) = ($info =~ /type=([^;]+);ref_chr=([^;]+);ref_start=([^;]+);ref_end=([^;]+);query_chr=([^;]+);query_start=([^;]+);query_end=([^;]+);query_orientation=([^;]+)/);
	# print "ref_id=$ref_id, query_id=$query_id, type=$type\n";
	my $adjusted_ref_chr;
	my $adjusted_query_chr;

	if ($ref_chr =~ /^${ref_tag}_(\S+)/) {
	    $adjusted_ref_chr = $1; 
	}
	if ($query_chr =~ /^${query_tag}_(\S+)/) {
	    $adjusted_query_chr = $1; 
	}

	# check if the paired ref and query markers locate on the same chromosome
	if ($adjusted_ref_chr eq $adjusted_query_chr) {
	    # then check if the paired ref and query markers locate on the same side of this chromosome arm
	    # if no, this is a dubious marker pair and will be excluded
	    my $flag = 0;
	    if (($ref_end < $ref_centromere_gff{$ref_chr}{'start'}) and ($query_start > $query_centromere_gff{$query_chr}{'end'})) {
		$flag = 1;
	    } elsif (($ref_start > $ref_centromere_gff{$ref_chr}{'end'}) and ($query_end < $query_centromere_gff{$query_chr}{'start'})) {
		$flag = 1;
	    }
	    
	    if ($flag == 0) {
		print $output_fh "$_\n";
	    } else {
		# print "filter out dubious marker: $_\n";
		$filtered_count++;
	    }
	} 
    }
}

print "\nfiltered out $filtered_count markers that fall on opposite sides of chromosome arms!\n";


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

sub parse_gff_file {
    my $fh = shift @_;
    my %gff = ();
    while (<$fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	my ($chr, $source, $type, $start, $end, $score, $strand, $frame, $info) = split /\t/, $_;
	my ($id) = ($info =~ /ID=([^;]+)/);
	$gff{$chr}{'chr'} = $chr;
	$gff{$chr}{'start'} = $start;
	$gff{$chr}{'end'} = $end;
	$gff{$chr}{'strand'} = $strand;
    }
    return %gff;
}



    



