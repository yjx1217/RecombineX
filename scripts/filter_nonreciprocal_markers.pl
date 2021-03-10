#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: filter_nonreciprocal_makers.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.08.17
#  description: filter out nonreciprocal markers in the ref-query and query-ref comparison.
#  example: perl filter_nonreciprocal_makers.pl -g1_based_vcf genome1_based.vcf(.gz) -g2_based_vcf genome2_based.vcf(.gz) -g1 genome1.fa(.gz) -g2 genome2.fa(.gz)  -g1_based_prefix g1_based -g2_based_prefix g2_based
##############################################################

my ($genome1_based_vcf, $genome2_based_vcf, $genome1_fasta, $genome2_fasta, $g1_based_prefix, $g2_based_prefix);

GetOptions('g1_based_vcf|genome1_based_vcf:s' => \$genome1_based_vcf,
	   'g2_based_vcf|genome2_based_vcf:s' => \$genome2_based_vcf,
	   'g1|genome1_fasta:s' => \$genome1_fasta,
	   'g2|genome2_fasta:s' => \$genome2_fasta,
	   'g1_based_prefix|g1_based_prefix:s' => \$g1_based_prefix,
	   'g2_based_prefix|g2_based_prefix:s' => \$g2_based_prefix,
    );

my $genome1_based_vcf_fh = read_file($genome1_based_vcf);
my $genome2_based_vcf_fh = read_file($genome2_based_vcf);
my %genome1_based_vcf = parse_vcf_file($genome1_based_vcf_fh);
close $genome1_based_vcf_fh;
my %genome2_based_vcf = parse_vcf_file($genome2_based_vcf_fh);
close $genome2_based_vcf_fh;

my $genome1_fasta_fh = read_file($genome1_fasta);
my %genome1_fasta = ();
my @genome1_fasta = ();
parse_fasta_file($genome1_fasta_fh, \%genome1_fasta, \@genome1_fasta);

my $genome2_fasta_fh = read_file($genome2_fasta);
my %genome2_fasta = ();
my @genome2_fasta = ();
parse_fasta_file($genome2_fasta_fh, \%genome2_fasta, \@genome2_fasta);


my $output_genome1_based_txt = "$g1_based_prefix.markers.txt";
my $output_genome1_based_txt_fh = write_file($output_genome1_based_txt);
print $output_genome1_based_txt_fh "parent1_chr\tparent1_start\tparent1_end\tparent1_allele\tparent2_allele\tparent2_chr\tparent2_start\tparent2_end\trelative_orientation\tmarker_type\n";

my $output_genome1_based_vcf = "$g1_based_prefix.markers.vcf";
my $output_genome1_based_vcf_fh = write_file($output_genome1_based_vcf);

$genome1_based_vcf_fh = read_file($genome1_based_vcf);
while (<$genome1_based_vcf_fh>) {
    chomp;
    /^\s*$/ and next;
     if (/^##source=/) {
        print $output_genome1_based_vcf_fh "##source=filter_nonreciprocal_markers.pl\n";
    } elsif (/^#/) {
	print $output_genome1_based_vcf_fh "$_\n";
    } else {
	my ($chr, $pos, $id, $ref_allele, $alt_allele, $qual, $filter, $info) = split /\t/, $_;
	my ($type, $ref_chr, $ref_start, $ref_end, $query_chr, $query_start, $query_end, $query_orientation) = ($info =~ /type=([^;]+);ref_chr=([^;]+);ref_start=([^;]+);ref_end=([^;]+);query_chr=([^;]+);query_start=([^;]+);query_end=([^;]+);query_orientation=([^;]+)/);
	my $match_key;
	my $marker_type;
	if ($type eq "SNP") {
	    $marker_type = "SNP";
	} else {
	    $marker_type = "INDEL";
	}
	if ($query_orientation eq "1") {
	    $match_key = "${query_chr}:${query_start}-${query_end}";
	    if (exists $genome2_based_vcf{$match_key}) {
		if ($ref_chr eq $genome2_based_vcf{$match_key}{'query_chr'}) {
		    if (($ref_start eq $genome2_based_vcf{$match_key}{'query_start'}) and ($ref_end eq $genome2_based_vcf{$match_key}{'query_end'})) {
			if ($query_orientation eq $genome2_based_vcf{$match_key}{'query_orientation'}) {
			    if (($ref_start eq $genome2_based_vcf{$match_key}{'query_start'}) and ($ref_end eq $genome2_based_vcf{$match_key}{'query_end'})) {
				if (($ref_allele eq $genome2_based_vcf{$match_key}{'alt_allele'}) and ($alt_allele eq $genome2_based_vcf{$match_key}{'ref_allele'})) {
				    print $output_genome1_based_txt_fh "$ref_chr\t$ref_start\t$ref_end\t$ref_allele\t$alt_allele\t$query_chr\t$query_start\t$query_end\t$query_orientation\t$marker_type\n";
				    print $output_genome1_based_vcf_fh "$_\n";
				}
			    }
			}
		    }
		}
	    }
	} else {
	    # adjustment for inversions
	    if ($type eq "SNP") {
		$match_key = "${query_chr}:${query_start}-${query_end}";
		if (exists $genome2_based_vcf{$match_key}) {
		    if ($ref_chr eq $genome2_based_vcf{$match_key}{'query_chr'}) {
			if (($ref_start eq $genome2_based_vcf{$match_key}{'query_start'}) and ($ref_end eq $genome2_based_vcf{$match_key}{'query_end'})) {
			    if ($query_orientation eq $genome2_based_vcf{$match_key}{'query_orientation'}) {
				if ((revcom($ref_allele) eq $genome2_based_vcf{$match_key}{'alt_allele'}) and (revcom($alt_allele) eq $genome2_based_vcf{$match_key}{'ref_allele'})) {
				    print $output_genome1_based_txt_fh "$ref_chr\t$ref_start\t$ref_end\t$ref_allele\t$alt_allele\t$query_chr\t$query_start\t$query_end\t$query_orientation\t$marker_type\n";
				    print $output_genome1_based_vcf_fh "$_\n";
				}
			    }
			}
		    }
		}
	    } elsif ($type =~ /(INSERTION|DELETION)/) {
		# further adjustment for indels in on inverted regions
		my $adjusted_ref_start = $ref_start + 1;
		my $adjusted_ref_end = $ref_end + 1;
		my $adjusted_query_start = $query_start - 1;
		my $adjusted_query_end = $query_end - 1;
		my $adjusted_ref_allele = substr $genome1_fasta{$ref_chr}, $adjusted_ref_start-1, $adjusted_ref_end-$adjusted_ref_start+1;
		my $adjusted_alt_allele = substr $genome2_fasta{$query_chr}, $adjusted_query_start-1, $adjusted_query_end-$adjusted_query_start+1;
		$adjusted_alt_allele = revcom($adjusted_alt_allele);
		$match_key = "${query_chr}:${adjusted_query_start}-${adjusted_query_end}";
		if (exists $genome2_based_vcf{$match_key}) {
		    if ($ref_chr eq $genome2_based_vcf{$match_key}{'query_chr'}) {
			if (($adjusted_ref_start eq $genome2_based_vcf{$match_key}{'query_start'}) and ($adjusted_ref_end eq $genome2_based_vcf{$match_key}{'query_end'})) {
			    if ($query_orientation eq $genome2_based_vcf{$match_key}{'query_orientation'}) {
				if ((revcom($adjusted_ref_allele) eq $genome2_based_vcf{$match_key}{'alt_allele'}) and (revcom($adjusted_alt_allele) eq $genome2_based_vcf{$match_key}{'ref_allele'})) {
				    print $output_genome1_based_txt_fh "$ref_chr\t$ref_start\t$ref_end\t$ref_allele\t$alt_allele\t$query_chr\t$query_start\t$query_end\t$query_orientation\t$marker_type\n";
				    print $output_genome1_based_vcf_fh "$_\n";
				}
			    }
			}		    
		    }
		}
	    }
	}
    }
}


#########


my $output_genome2_based_txt = "$g2_based_prefix.markers.txt";
my $output_genome2_based_txt_fh = write_file($output_genome2_based_txt);
print $output_genome2_based_txt_fh "parent2_chr\tparent2_start\tparent2_end\tparent2_allele\tparent1_allele\tparent1_chr\tparent1_start\tparent1_end\trelative_orientation\tmarker_type\n";

my $output_genome2_based_vcf = "$g2_based_prefix.markers.vcf";
my $output_genome2_based_vcf_fh = write_file($output_genome2_based_vcf);

$genome2_based_vcf_fh = read_file($genome2_based_vcf);
while (<$genome2_based_vcf_fh>) {
    chomp;
    /^\s*$/ and next;
    if (/^##source=/) {
        print $output_genome2_based_vcf_fh "##source=filter_nonreciprocal_markers.pl\n";
    } elsif (/^#/) {
	print $output_genome2_based_vcf_fh "$_\n";
    } else {
	my ($chr, $pos, $id, $ref_allele, $alt_allele, $qual, $filter, $info) = split /\t/, $_;
	my ($type, $ref_chr, $ref_start, $ref_end, $query_chr, $query_start, $query_end, $query_orientation) = ($info =~ /type=([^;]+);ref_chr=([^;]+);ref_start=([^;]+);ref_end=([^;]+);query_chr=([^;]+);query_start=([^;]+);query_end=([^;]+);query_orientation=([^;]+)/);
	my $match_key;
	my $marker_type;
	if ($type eq "SNP") {
	    $marker_type = "SNP";
	} else {
	    $marker_type = "INDEL";
	}
	if ($query_orientation eq "1") {
	    $match_key = "${query_chr}:${query_start}-${query_end}";
	    if (exists $genome1_based_vcf{$match_key}) {
		if ($ref_chr eq $genome1_based_vcf{$match_key}{'query_chr'}) {
		    if (($ref_start eq $genome1_based_vcf{$match_key}{'query_start'}) and ($ref_end eq $genome1_based_vcf{$match_key}{'query_end'})) {
			
			if ($query_orientation eq $genome1_based_vcf{$match_key}{'query_orientation'}) {
			    if (($ref_start eq $genome1_based_vcf{$match_key}{'query_start'}) and ($ref_end eq $genome1_based_vcf{$match_key}{'query_end'})) {
				if (($ref_allele eq $genome1_based_vcf{$match_key}{'alt_allele'}) and ($alt_allele eq $genome1_based_vcf{$match_key}{'ref_allele'})) {
				    print $output_genome2_based_txt_fh "$ref_chr\t$ref_start\t$ref_end\t$ref_allele\t$alt_allele\t$query_chr\t$query_start\t$query_end\t$query_orientation\t$marker_type\n";
				    print $output_genome2_based_vcf_fh "$_\n";
				}
			    }
			}
		    }
		}
	    }
	} else {
	    # adjustment for inversions
	    if ($type eq "SNP") {
		$match_key = "${query_chr}:${query_start}-${query_end}";
		if (exists $genome1_based_vcf{$match_key}) {
		    if ($ref_chr eq $genome1_based_vcf{$match_key}{'query_chr'}) {
			if (($ref_start eq $genome1_based_vcf{$match_key}{'query_start'}) and ($ref_end eq $genome1_based_vcf{$match_key}{'query_end'})) {
			    
			    if ($query_orientation eq $genome1_based_vcf{$match_key}{'query_orientation'}) {
				if ((revcom($ref_allele) eq $genome1_based_vcf{$match_key}{'alt_allele'}) and (revcom($alt_allele) eq $genome1_based_vcf{$match_key}{'ref_allele'})) {
				    print $output_genome2_based_txt_fh "$ref_chr\t$ref_start\t$ref_end\t$ref_allele\t$alt_allele\t$query_chr\t$query_start\t$query_end\t$query_orientation\t$marker_type\n";
				    print $output_genome2_based_vcf_fh "$_\n";
				}
			    }
			}
		    }
		}
	    } elsif ($type =~ /(INSERTION|DELETION)/) {
		# further adjustment for indels in on inverted regions
		my $adjusted_ref_start = $ref_start + 1;
		my $adjusted_ref_end = $ref_end + 1;
		my $adjusted_query_start = $query_start - 1;
		my $adjusted_query_end = $query_end - 1;
		my $adjusted_ref_allele = substr $genome2_fasta{$ref_chr}, $adjusted_ref_start-1, $adjusted_ref_end-$adjusted_ref_start+1;
		my $adjusted_alt_allele = substr $genome1_fasta{$query_chr}, $adjusted_query_start-1, $adjusted_query_end-$adjusted_query_start+1;
		$adjusted_alt_allele = revcom($adjusted_alt_allele);
		$match_key = "${query_chr}:${adjusted_query_start}-${adjusted_query_end}";
		if (exists $genome1_based_vcf{$match_key}) {
		    if ($ref_chr eq $genome1_based_vcf{$match_key}{'query_chr'}) {
			if (($adjusted_ref_start eq $genome1_based_vcf{$match_key}{'query_start'}) and ($adjusted_ref_end eq $genome1_based_vcf{$match_key}{'query_end'})) {
			    
			    if ($query_orientation eq $genome1_based_vcf{$match_key}{'query_orientation'}) {
				if ((revcom($adjusted_ref_allele) eq $genome1_based_vcf{$match_key}{'alt_allele'}) and (revcom($adjusted_alt_allele) eq $genome1_based_vcf{$match_key}{'ref_allele'})) {
				    print $output_genome2_based_txt_fh "$ref_chr\t$ref_start\t$ref_end\t$ref_allele\t$alt_allele\t$query_chr\t$query_start\t$query_end\t$query_orientation\t$marker_type\n";
				    print $output_genome2_based_vcf_fh "$_\n";
				}
			    }
			}
		    }			    
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
    if ($file =~ /\.gz$/) {
        open($fh, "| gzip -c >$file") or die "can't open $file\n";
    } else {
        open($fh, ">$file") or die "can't open $file\n";
    }
    return $fh;
}

sub parse_vcf_file {
    my $fh = shift @_;
    my %vcf = ();
    while (<$fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	my ($chr, $pos, $id, $ref_allele, $alt_allele, $qual, $filter, $info) = split /\t/, $_;
	my ($type, $ref_chr, $ref_start, $ref_end, $query_chr, $query_start, $query_end, $query_orientation) = ($info =~ /type=([^;]+);ref_chr=([^;]+);ref_start=([^;]+);ref_end=([^;]+);query_chr=([^;]+);query_start=([^;]+);query_end=([^;]+);query_orientation=([^;]+)/);
	my $match_key = "${ref_chr}:${ref_start}-${ref_end}";
	$vcf{$match_key}{'ref_allele'} = $ref_allele;
	$vcf{$match_key}{'alt_allele'} = $alt_allele;
	$vcf{$match_key}{'qual'} = $qual;
	$vcf{$match_key}{'ref_chr'} = $ref_chr;
	$vcf{$match_key}{'ref_start'} = $ref_start;
	$vcf{$match_key}{'ref_end'} = $ref_end;
	$vcf{$match_key}{'query_chr'} = $query_chr;
	$vcf{$match_key}{'query_start'} = $query_start;
	$vcf{$match_key}{'query_end'} = $query_end;
	$vcf{$match_key}{'query_orientation'} = $query_orientation;
    }
    return %vcf;
}


sub revcom {
    my $seq = shift @_;
    my $seq_revcom = reverse $seq;
    $seq_revcom =~ tr/ATGCatgc/TACGtacg/;
    return $seq_revcom;
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



