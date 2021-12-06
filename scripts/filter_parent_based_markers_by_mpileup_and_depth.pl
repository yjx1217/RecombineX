#!/usr/bin/perl
use warnings FATAL => 'all';;
use strict;
use Getopt::Long;

##############################################################
#  script: filter_parent_based_markers_by_mpileup_and_depth.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2021.06.17
#  description: filter out markers with unexpected depth (<0.5 genome-wide median depth or >1.5 genome-wide median depth)
#  example: perl filter_parent_based_markers_by_mpileup_and_depth.pl -depth_detail mpileup.gz -depth_summary depth.summary.txt -i markers.vcf.gz -o filtered.markers.vcf.gz
##############################################################

my ($input, $output, $depth_summary, $mpileup);

GetOptions('i|input:s' => \$input,
	   'depth_summary|depth_summary:s' => \$depth_summary,
	   'mpileup|mpileup:s' => \$mpileup,
	   'o|output:s' => \$output);

my $marker_purity_cutoff = 0.9;
my $depth_summary_fh = read_file($depth_summary);
my %depth_summary = parse_depth_summary_file($depth_summary_fh);

my $output_fh = write_file($output);

my $filtered_count = 0;

my $i = 0;
my %marker_list = ();
my %check_list = ();

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

    $i++;
    $marker_list{$i} = $_;
    $check_list{$ref_chr}{$ref_start} = 1;
    $check_list{$ref_chr}{$ref_end} = 1;
}


my $mpileup_fh = read_file($mpileup);
my %mpileup = ();
parse_mpileup_file($mpileup_fh, $marker_purity_cutoff, \%check_list, \%mpileup);

for (my $j = 1; $j <= $i; $j++) {
    my @line = split /\t/, $marker_list{$j};
    my $ref_chr = $line[0];
    my $ref_start = $line[1];
    my $ref_allele = $line[3];
    my $query_allele = $line[4];
    my ($query_chr) = ($marker_list{$j} =~ /query_chr=([^;]+)/);
    my ($query_start) = ($marker_list{$j} =~ /query_start=([^;]+)/);
    my ($query_end) = ($marker_list{$j} =~ /query_end=([^;]+)/);
    my ($ref_end) = ($marker_list{$j} =~ /ref_end=([^;]+)/);
    # print "ref_chr=$ref_chr, ref_start=$ref_start, ref_allele=$ref_allele, query_allele=$query_allele\n";
    if (exists $mpileup{$ref_chr}{$ref_start}{'allele'}) {
	if ($query_allele eq $mpileup{$ref_chr}{$ref_start}{'allele'}) {
	    # print "PASS: ref_chr=$ref_chr, ref_start=$ref_start, ref_allele=$ref_allele, query_allele=$query_allele, mpileup_allele=$mpileup{$ref_chr}{$ref_start}{'allele'}\n";
	    my $flag = 1;
	    if ((exists $mpileup{$ref_chr}{$ref_start}{'depth'}) and ($mpileup{$ref_chr}{$ref_end}{'depth'})) {
		if (($mpileup{$ref_chr}{$ref_start}{'depth'} > $depth_summary{$ref_chr} * 0.5) and ($mpileup{$ref_chr}{$ref_start}{'depth'} < $depth_summary{$ref_chr} * 1.5)) {
		    if (($mpileup{$ref_chr}{$ref_end}{'depth'} > $depth_summary{$ref_chr} * 0.5) and ($mpileup{$ref_chr}{$ref_end}{'depth'} < $depth_summary{$ref_chr} * 1.5)) {
			$flag = 0;
		    }
		}
	    }
	    if ($flag == 0) {
		print $output_fh "$marker_list{$j}\n";
	    } else {
		$filtered_count++;
	    }
	} else {
	    # print "Filter out this marker: ref_chr=$ref_chr, ref_start=$ref_start, ref_allele=$ref_allele, query_allele=$query_allele, mpileup_allele=$mpileup{$ref_chr}{$ref_start}{'allele'}\n";
	    $filtered_count++;
	}
    }
}

print "filtered $filtered_count markers due to false negative in variant calling against referenc and/or substantial depth bias!\n";

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


sub parse_mpileup_file {
    my ($fh, $marker_purity_cutoff, $markers_hashref, $mpileup_hashref) = @_;
    while (<$fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	my ($chr, $pos, $ref_allele, $depth, $reads_bases, $reads_quals) = split /\t/, $_;
	$$mpileup_hashref{$chr}{$pos}{'depth'} = $depth;
	if (exists $$markers_hashref{$chr}{$pos}) {
	    if ($depth > 0) {
		# print "raw mpileup: $_\n";                                                                                                                                        # remove the start and end marks of reads mapping
		$reads_bases =~ s/\^\S//gi;
		$reads_bases =~ s/\$//gi;
		# print "reads_bases: $reads_bases\n";
		my %indels = ();
		my $j = 0; # indel index         
		while ($reads_bases =~ /((?:\S)(?:\+|\-)[0-9]+[ACGTNacgtn]+)/g) {
		    $j++;
		    my $indel_raw = $&;
		    # my $match_start = $-[0] + 1;
		    # my $match_end = $+[0];
		    $indels{$j}{'indel_raw'} = $&;
		    ($indels{$j}{'lead_base'}, $indels{$j}{'type'}, $indels{$j}{'length'}, $indels{$j}{'indel_adjusted'}) = ($indels{$j}{'indel_raw'} =~ /(\S)(\+|\-)([0-9]+)([ACGTNacgtn]+)/);
		    $indels{$j}{'indel_adjusted'} = substr($indels{$j}{'indel_adjusted'}, 0, $indels{$j}{'length'});
		    if ($indels{$j}{'lead_base'} =~ /(\.|\,)/) {
			$indels{$j}{'lead_base'} = $ref_allele;
		    } else {
			# print "corner case! lead_base = $indels{$j}{'lead_base'}\n";
			$indels{$j}{'lead_base'} = uc $indels{$j}{'lead_base'};
		    }
		    if ($indels{$j}{'type'} eq '+') {
			$indels{$j}{'indel_adjusted'} = $indels{$j}{'lead_base'}. $indels{$j}{'indel_adjusted'};
			$reads_bases =~ s/(?:\S)(?:\+|\-)[0-9]+[ACGTNacgtn]{$indels{$j}{'length'}}/I/;
		    } else {
			$indels{$j}{'indel_adjusted'} = $indels{$j}{'lead_base'};
			$reads_bases =~ s/(?:\S)(?:\+|\-)[0-9]+[ACGTNacgtn]{$indels{$j}{'length'}}/D/;
		    }
		    # print "reads_bases after replacement: $reads_bases\n";
		}
		my %basecall = ();
                $j = 0;
                for (my $i = 1; $i <= $depth; $i++) {
                    my $base = substr $reads_bases, $i - 1, 1;
                    my $qual = substr $reads_quals, $i - 1, 1;
                    $qual = ord($qual) - 33;
		    # set safe bound for $qual
		    if ($qual > 40) {
                        $qual =	40;
                    } elsif ($qual < 0)	{
                        $qual =	0;
                    }
                    if (($base eq '.') or ($base eq ',')) {
                        $base = $ref_allele; # 20210226
		    } elsif ($base =~ /(I|D)/) {
                        # if indel
			$j++;
                        $base = uc $indels{$j}{'indel_adjusted'};
                    } elsif ($base =~ /[ATGCatgc]/) {
                        $base = uc $base;
                    } else {
                        # print "corner case: $base\n";
                        $base = 'NA';
                    }
		    
                    if (not exists $basecall{$base}) {
                        @{$basecall{$base}{'index'}} = ($i);
                        @{$basecall{$base}{'qual'}} = ($qual);
                        $basecall{$base}{'count'} = 1;
                    } else {
                        push @{$basecall{$base}{'index'}}, $i;
                        push @{$basecall{$base}{'qual'}}, $qual;
                        $basecall{$base}{'count'}++;
                    }
                    # print "i=$i: base=$base, basecall_index = @{$basecall{$base}{'index'}}, basecall_qual=@{$basecall{$base}{'qual'}}\n";
                }
                # print "\n";
                foreach my $base (sort keys %basecall) {
                    $basecall{$base}{'qual_sum'} = cal_sum(@{$basecall{$base}{'qual'}});
                }
		
                my @basecall_sorted = sort {$basecall{$b}{'qual_sum'} <=> $basecall{$a}{'qual_sum'} or $basecall{$b}{'count'} <=> $basecall{$a}{'count'}} keys %basecall;
                # print "basecall_sorted = @basecall_sorted\n";
                my $basecall_major = shift @basecall_sorted;
                # print "basecall_major = $basecall_major\n";
		my $basecall_major_count = $basecall{$basecall_major}{'count'};
		my $basecall_major_purity = $basecall_major_count/$depth;
		if ($basecall_major_purity > $marker_purity_cutoff) {
		    $$mpileup_hashref{$chr}{$pos}{'allele'} = $basecall_major;
		} 
	    }
        }
    }
}




sub cal_sum {
    my @input = @_;
    my $sum = 0;
    foreach my $n (@input) {
        $sum += $n;
    }
    return $sum;
}
