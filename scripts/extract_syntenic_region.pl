#!/usr/bin/perl
use warnings FATAL => 'all';
use strict;

my @chr = qw(chrI chrII chrIII chrIV chrV chrVI chrVII chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI);

my @cross_pair = qw(DBVPG6044-DBVPG6765 DBVPG6044-Y12 DBVPG6044-YPS128 DBVPG6765-Y12 DBVPG6765-YPS128 Y12-YPS128);

foreach my $cross_pair (@cross_pair) {
    my @ref = split /-/, $cross_pair;
    foreach my $ref (@ref) {
	my $prefix = "$cross_pair.$ref";
	my $output = "$prefix.colinear_region.gff3";
	my $output_fh = write_file($output);
	foreach my $chr (@chr) {
	    my $input = "./../12.Polymorphic_Markers_by_Cross_Parent_Genome_Alignment/${cross_pair}_${ref}_based/$cross_pair.$ref.$chr.1coords";
	    my $input_fh = read_file($input);
	    while (<$input_fh>) {
		next if $. <= 4; # Skip first 4 lines
		chomp;
		my ($ref_start, $ref_end, $query_start, $query_end, $ref_aln_length, $query_aln_length, $aln_identity, $ref_length, $query_length, $ref_coverage, $query_coverage, $ref_orientation, $query_orientation, $ref_chr, $query_chr) = split /\s+/, $_;
		if (($ref_orientation eq "1") and ($query_orientation eq "1")) {
		    print $output_fh "$ref_chr\tparent_genome_aln\tcollinear_region\t$ref_start\t$ref_end\t.\t+\t.\tquery_chr=$query_chr;query_start=$query_start;query_end=$query_end;aln_identity=$aln_identity;\n";
		}
	    }
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
    if ($file =~ /.gz$/) {
        open($fh, "| gzip -c >$file") or die "can't open $file\n";
    } else {
        open($fh, ">$file") or die "can't open $file\n";
    }
    return $fh;
}  
