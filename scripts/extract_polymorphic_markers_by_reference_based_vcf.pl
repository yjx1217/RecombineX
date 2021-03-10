#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: extract_polymorphic_markers_by_reference_based_vcf.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2019.09.06
#  description: extract polymorphic markers based on two VCF files against the same reference genome.
#  example: perl extract_polymorphic_markers_by_reference_based_vcf.pl -r refseq.fa(.gz) -v1 parent1.vcf(.gz) -v2 parent2.vcf(.gz) -b0 ref.repeatmask.bed(.gz) -b1 parent1.cnv.bed -b2 parent2.cnv.bed -q 30  -o output.markers.txt(.gz)
##############################################################

my ($refseq, $vcf1, $vcf2, $bed0, $bed1, $bed2, $q, $output, $check_repeat);
$q = 30;
$check_repeat = "no";
GetOptions('r|refseq:s' => \$refseq,
	   'v1|vcf1:s' => \$vcf1,
	   'v2|vcf2:s' => \$vcf2,
	   'b0|bed0:s' => \$bed0,
	   'b1|bed1:s' => \$bed1,
	   'b2|bed2:s' => \$bed2,
	   'check_repeat|check_repeat:s' => \$check_repeat,
	   'q|quality:s' => \$q,
	   'o|output:s' => \$output);

my $refseq_fh = read_file($refseq);
my %refseq = ();
my @refseq = ();
parse_fasta_file($refseq_fh, \%refseq, \@refseq);
    
my $vcf1_fh = read_file($vcf1);
my %vcf1 = parse_vcf_file($vcf1_fh);
my $vcf2_fh = read_file($vcf2);
my %vcf2 = parse_vcf_file($vcf2_fh);

my $bed0_fh = read_file($bed0);
my %bed0 = parse_bed_file($bed0_fh);

my $bed1_fh = read_file($bed1);
my %bed1 = parse_bed_file($bed1_fh);
my $bed2_fh = read_file($bed2);
my %bed2 = parse_bed_file($bed2_fh);

my $output_fh = write_file($output);
print $output_fh "ref_chr\tref_start\tref_end\tparent1_allele\tparent2_allele\tref_chr\tref_start\tref_end\trelative_orientation\tmarker_type\n";


foreach my $chr (@refseq) {
    my @start1 = sort {$a <=> $b} keys %{$vcf1{$chr}};
    my @start2 = sort {$a <=> $b} keys %{$vcf2{$chr}};
    my @start = sort {$a <=> $b} keys %{{map {($_ => 1)} (@start1, @start2)}};
    foreach my $start (@start) {
	if ((exists $vcf1{$chr}{$start}) and (exists $vcf2{$chr}{$start})) {
	    if (($vcf1{$chr}{$start}{'ref_allele'} eq $vcf2{$chr}{$start}{'ref_allele'}) and ($vcf1{$chr}{$start}{'alt_allele'} ne $vcf2{$chr}{$start}{'alt_allele'})) {
		if (($vcf1{$chr}{$start}{'qual'} >= $q) and ($vcf1{$chr}{$start}{'qual'} >= $q)) {
		    if ($vcf1{$chr}{$start}{'marker_type'} eq $vcf1{$chr}{$start}{'marker_type'}) {
			my $end = $start + length($vcf1{$chr}{$start}{'ref_allele'}) - 1;
			my $repeat_flag = 0;
			if ($check_repeat eq "yes") {
			    $repeat_flag = check_overlap_with_bed(\%bed0, $chr, $start, $end);
			}
			if ($repeat_flag == 0) {
			    my $cnv_flag1 = check_overlap_with_bed(\%bed1, $chr, $start, $end);
			    my $cnv_flag2 = check_overlap_with_bed(\%bed2, $chr, $start, $end);
			    if (($cnv_flag1 == 0) and ($cnv_flag2 == 0)) {
				print $output_fh "$chr\t$start\t$end\t$vcf1{$chr}{$start}{'alt_allele'}\t$vcf2{$chr}{$start}{'alt_allele'}\t$chr\t$start\t$end\t1\t$vcf1{$chr}{$start}{'marker_type'}\n";
			    }
			}
		    }
		}
	    }
	} elsif ((exists $vcf1{$chr}{$start}) and (not exists $vcf2{$chr}{$start})) {
	    if ($vcf1{$chr}{$start}{'qual'} >= $q) {
		my $end = $start + length($vcf1{$chr}{$start}{'ref_allele'}) - 1;
		my $repeat_flag = 0;
		if ($check_repeat eq "yes") {
		    $repeat_flag = check_overlap_with_bed(\%bed0, $chr, $start, $end);
		}
		if ($repeat_flag == 0) {
		    my $cnv_flag1 = check_overlap_with_bed(\%bed1, $chr, $start, $end);
		    my $cnv_flag2 = check_overlap_with_bed(\%bed2, $chr, $start, $end);
		    if (($cnv_flag1 == 0) and ($cnv_flag2 == 0)) {
			print $output_fh "$chr\t$start\t$end\t$vcf1{$chr}{$start}{'alt_allele'}\t$vcf1{$chr}{$start}{'ref_allele'}\t$chr\t$start\t$end\t1\t$vcf1{$chr}{$start}{'marker_type'}\n";
		    }
		}
	    }
	} else {
	    if ($vcf2{$chr}{$start}{'qual'} >= $q) {
		my $end = $start + length($vcf2{$chr}{$start}{'ref_allele'}) - 1;
		my $repeat_flag = 0;
		if ($check_repeat eq "yes") {
		    $repeat_flag = check_overlap_with_bed(\%bed0, $chr, $start, $end);
		}
		if ($repeat_flag == 0) {
		    my $cnv_flag1 = check_overlap_with_bed(\%bed1, $chr, $start, $end);
		    my $cnv_flag2 = check_overlap_with_bed(\%bed2, $chr, $start, $end);
		    if (($cnv_flag1 == 0) and ($cnv_flag2 == 0)) {
			print $output_fh "$chr\t$start\t$end\t$vcf2{$chr}{$start}{'ref_allele'}\t$vcf2{$chr}{$start}{'alt_allele'}\t$chr\t$start\t$end\t1\t$vcf2{$chr}{$start}{'marker_type'}\n";
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


sub parse_vcf_file {
    my $fh = shift @_;
    my %vcf = ();
    while (<$fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	my ($chr, $pos, $id, $ref_allele, $alt_allele, $qual, $filter, $info) = split /\t/, $_;
	$vcf{$chr}{$pos}{'ref_allele'} = $ref_allele;
	$vcf{$chr}{$pos}{'alt_allele'} = $alt_allele;
	$vcf{$chr}{$pos}{'qual'} = $qual;
	$vcf{$chr}{$pos}{'info'} = $info;
	$vcf{$chr}{$pos}{'marker_type'} = "NA";
	if ($info =~ /VT=(SNP|INDEL)/) {
	    $vcf{$chr}{$pos}{'marker_type'} = $1;
	}
    }
    return %vcf;
}

sub parse_bed_file {
    my $fh = shift @_;
    my %bed = ();
    while (<$fh>) {
        chomp;
        /^#/ and next;
	/^\s*$/ and next;
	/^chr\tstart\tend/ and next;
	my ($chr, $start, $end) = split /\t/, $_;
	$bed{$chr}{$start + 1}{$end} = 1;
    }
    return %bed;
}

sub check_overlap_with_bed {
    my ($bed_hashref, $chr, $start, $end) = @_;
    my $cnv_flag = 0;
    if (exists $$bed_hashref{$chr}) {
	foreach my $s (sort {$a <=> $b} keys %{$$bed_hashref{$chr}}) {
	    if (($cnv_flag == 0) and ($s <= $end)) {
		foreach my $e (sort {$a <=> $b} keys %{$$bed_hashref{$chr}{$s}}) {
		    if ($e >= $start) {
			$cnv_flag = 1;
			last;
		    } 
		}
	    }
	}
    }
    return $cnv_flag;
}
