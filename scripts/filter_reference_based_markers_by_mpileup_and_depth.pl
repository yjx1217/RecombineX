#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: filter_reference_based_markers_by_mpileup_and_depth.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.05.03
#  description: filter out markers with false negative variant call against reference and/or with unexpected depth (<0.5 genome-wide median depth or >1.5 genome-wide median depth)
#  example: perl filter_reference_based_markers_by_mpileup_and_depth.pl -ref_mpileup ref.mpileup.gz -ref_depth_summary ref.depth.summary.txt -query_mpileup query.mpileup.gz -query_depth_summary query.depth.summary.txt_ -i marker_table -o filtered_marker_table -filter_by_depth_variation yes
##############################################################


my ($input, $output, $filter_by_depth_variation, $ref_depth_summary, $ref_mpileup, $query_depth_summary, $query_mpileup);

GetOptions('i|input:s' => \$input,
           'filter_by_depth_variation|filter_by_depth_variation:s' => \$filter_by_depth_variation,
	   'ref_depth_summary|ref_depth_summary:s' => \$ref_depth_summary,
	   'query_depth_summary|query_depth_summary:s' => \$query_depth_summary,
	   'ref_mpileup|ref_mpileup:s' => \$ref_mpileup,
	   'query_mpileup|query_mpileup:s' => \$query_mpileup,
	   'o|output:s' => \$output);


my $marker_purity_cutoff = 0.9;
my $ref_depth_summary_fh;
my %ref_depth_summary;
my $query_depth_summary_fh;
my %query_depth_summary;

if ($filter_by_depth_variation eq "yes") {
    $ref_depth_summary_fh = read_file($ref_depth_summary);
    %ref_depth_summary = parse_depth_summary_file($ref_depth_summary_fh);

    $query_depth_summary_fh = read_file($query_depth_summary);
    %query_depth_summary = parse_depth_summary_file($query_depth_summary_fh); 
}

my $output_fh = write_file($output);

my $filtered_count = 0;

my $input_fh = read_file($input);
my $i = 0;
my %marker_list = ();
my %ref_check_list = ();
my %query_check_list= ();

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
    $i++;
    $marker_list{$i} = $_;
    $ref_check_list{$ref_chr}{$ref_start} = 1;
    $ref_check_list{$ref_chr}{$ref_end} = 1;
    $query_check_list{$ref_chr}{$query_start} = 1;
    $query_check_list{$ref_chr}{$query_end} = 1;
}


my $ref_mpileup_fh = read_file($ref_mpileup);
my %ref_mpileup = ();
parse_mpileup_file($ref_mpileup_fh, $marker_purity_cutoff, \%ref_check_list, \%ref_mpileup);

my $query_mpileup_fh = read_file($query_mpileup);
my %query_mpileup = ();
parse_mpileup_file($query_mpileup_fh, $marker_purity_cutoff, \%query_check_list, \%query_mpileup);

for (my $j = 1; $j <= $i; $j++) {
    my ($ref_chr, $ref_start, $ref_end, $ref_allele, $query_allele, $query_chr, $query_start, $query_end, $match_orientation) = split /\t/, $marker_list{$j};
    if ((exists $ref_mpileup{$ref_chr}{$ref_start}{'allele'}) and (exists $query_mpileup{$query_chr}{$query_start}{'allele'})) {
	if (($ref_allele eq $ref_mpileup{$ref_chr}{$ref_start}{'allele'}) and ($query_allele eq $query_mpileup{$query_chr}{$query_start}{'allele'})) {
            if ($filter_by_depth_variation eq "yes") {
		my $flag_by_ref_depth = 1;
		my $flag_by_query_depth = 1;
		if ((exists $ref_mpileup{$ref_chr}{$ref_start}{'depth'}) and (exists $ref_mpileup{$ref_chr}{$ref_end}{'depth'})) {
		    if (($ref_mpileup{$ref_chr}{$ref_start}{'depth'} > $ref_depth_summary{$ref_chr} * 0.5) and ($ref_mpileup{$ref_chr}{$ref_start}{'depth'} < $ref_depth_summary{$ref_chr} * 1.5)) {
			if (($ref_mpileup{$ref_chr}{$ref_end}{'depth'} > $ref_depth_summary{$ref_chr} * 0.5) and ($ref_mpileup{$ref_chr}{$ref_end}{'depth'} < $ref_depth_summary{$ref_chr} * 1.5)) {
			    $flag_by_ref_depth = 0;
			}
		    }
		}
		if ((exists $query_mpileup{$query_chr}{$query_start}{'depth'}) and (exists $query_mpileup{$query_chr}{$query_end}{'depth'})) {
		    if (($query_mpileup{$query_chr}{$query_start}{'depth'} > $query_depth_summary{$query_chr} * 0.5) and ($query_mpileup{$query_chr}{$query_start}{'depth'} < $query_depth_summary{$query_chr} * 1.5)) {
			if (($query_mpileup{$query_chr}{$query_end}{'depth'} > $query_depth_summary{$query_chr} * 0.5) and ($query_mpileup{$query_chr}{$query_end}{'depth'} < $query_depth_summary{$query_chr} * 1.5)) {
			    $flag_by_query_depth = 0;
			}
		    }
		}
		if (($flag_by_ref_depth == 0) and ($flag_by_query_depth == 0)) {
		    print $output_fh "$marker_list{$j}\n";
		} else {
		    $filtered_count++;
		}
	    } else {
                print $output_fh "$marker_list{$j}\n";
            }
	} else {
            # print "Filter out this marker: ref_chr=$ref_chr, ref_start=$ref_start, ref_allele=$ref_allele, query_allele=$query_allele, mpileup_allele=$mpileup{$ref_chr}{$ref_start}{'allele'}\n";
	    $filtered_count++;
	}
    }
}

print "filtered $filtered_count markers due to false negative in variant calling against reference and/or due to substantial depth bias!\n";

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
                # print "raw mpileup: $_\n";
                # remove the start and end marks of reads mapping
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
			$qual = 40;
		    } elsif ($qual < 0) {
			$qual = 0;
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
