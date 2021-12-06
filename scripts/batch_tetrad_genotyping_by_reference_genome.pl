#!/usr/bin/perl
use warnings FATAL => 'all';
use strict;
use Getopt::Long;
use Cwd;
# use Data::Dumper;   

##############################################################
#  script: batch_tetrad_genotying_by_reference_genome.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2019.09.09
#  description: run batch tetrad-based genotyping based on the mpileup files of spores
#  example: perl batch_tetrad_genotyping_by_reference_genome.pl -s Master_Sample_Table.txt -m SNP -c ./../../data/chr_list.yeast.txt -q 30 -marker_dir ./../02.Polymorphic_Markers_by_Reference_based_Read_Mapping -b batch_id -gamete_read_mapping_dir ./../03.Gamete_Read_Mapping_to_Reference_Genome -output_dir output_dir
##############################################################

my ($sample_table, $marker_type, $chr_list, $basecall_qual_diff_cutoff, $basecall_purity_cutoff, $marker_dir, $batch_id, $gamete_read_mapping_dir, $output_dir, $apply_cnv_filter, $allow_heteroduplex);
$marker_dir = "./../02.Polymorphic_Markers_by_Reference_based_Read_Mapping";
$gamete_read_mapping_dir = "./../03.Gamete_Read_Mapping_to_Reference_Genome";

$batch_id = "Batch_TEST";
$marker_type = "SNP";
$basecall_qual_diff_cutoff = 30;
$basecall_purity_cutoff = 0.9;
$output_dir = "$batch_id";
$apply_cnv_filter = "no";
$allow_heteroduplex = "no";

GetOptions('sample_table|s:s' => \$sample_table, # master sample list
           'marker_type|m:s' => \$marker_type, # SNP or INDEL or BOTH
	   'qual|q:i' => \$basecall_qual_diff_cutoff, # 30
	   'purity|p:f' => \$basecall_purity_cutoff, # 0.9
	   'chr_list|c:s' => \$chr_list,
	   'marker_dir|m_dir:s' => \$marker_dir,
	   'batch_id|b:s' => \$batch_id, 
	   'apply_cnv_filter:s' => \$apply_cnv_filter,
	   'allow_heteroduplex:s' => \$allow_heteroduplex,
	   'gamete_read_mapping_dir|orm_dir:s' => \$gamete_read_mapping_dir,
	   'output_dir|o:s' => \$output_dir);

my $sample_table_fh = read_file($sample_table);
my $chr_list_fh = read_file($chr_list);
my @chr = parse_list_file($chr_list_fh);
my %tetrads = parse_sample_table_by_tetrad($sample_table_fh);
close $sample_table_fh;

my $base_dir = cwd();
system("mkdir $output_dir");

foreach my $tetrad_id (sort keys %tetrads) {
    my $cross_pair = $tetrads{$tetrad_id}{'cross_pair'};
    my ($genome1_tag, $genome2_tag) = split "-", $cross_pair;
    print "tetrad: $tetrad_id, cross_pair: $cross_pair, parent1: $genome1_tag, parent2: $genome2_tag\n";
    $marker_type = uc $marker_type;
    my %markers = ();
    print "get markers ...\n";
    get_markers($genome1_tag, $genome2_tag, $marker_dir, $marker_type, \%markers);
    my %genotypes = ();
    my %mutations = ();
    my %heteroduplex = ();
    # my @spore_indices = qw(a b c d);
    my @spore_indices = @{$tetrads{$tetrad_id}{'spore_index'}};

    # raw genotype calling
    print "raw genotype calling ...\n";
    foreach my $spore_index (@spore_indices) {
	my $spore_id = "$tetrad_id.$spore_index";
	my $ref_based_spore_dir = "${gamete_read_mapping_dir}/${batch_id}/${cross_pair}.${spore_id}.ref";
	my $ref_based_spore_mpileup = "${ref_based_spore_dir}/${cross_pair}.${spore_id}.ref.mpileup.gz";
	if (-f $ref_based_spore_mpileup) {
	    my $ref_based_spore_mpileup_fh = read_file($ref_based_spore_mpileup);
	    parse_mpileup_file($spore_id, $spore_index, $cross_pair, $genome1_tag, $genome2_tag, $ref_based_spore_mpileup_fh, $allow_heteroduplex, \%markers, \%genotypes, \%mutations, \%heteroduplex);
	    if ($apply_cnv_filter eq "yes") {
		my $ref_based_spore_cnv = "${ref_based_spore_dir}/${cross_pair}.${spore_id}.ref.CNV_significance_test.txt";
		my $ref_based_spore_cnv_fh = read_file($ref_based_spore_cnv);
		filter_genotype_by_cnv($spore_id, $spore_index, $cross_pair, $genome1_tag, $genome2_tag, $ref_based_spore_cnv_fh, \%genotypes, \%mutations, \%heteroduplex);
	    }
	}
    }
    my $prefix;
    print "output genotypes based parent1: $genome1_tag\n";
    $prefix = "$base_dir/$output_dir/$cross_pair.$tetrad_id.ref.q${basecall_qual_diff_cutoff}";
    output_genotypes($prefix, $tetrad_id, $cross_pair, $genome1_tag, $genome2_tag, \@chr, \%genotypes, \%mutations, \%heteroduplex);
    
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


sub parse_sample_table_by_tetrad {
    my $fh = shift @_;
    my %sample_table = ();
    while (<$fh>) {
        chomp;
        /^#/ and next;
	/^\s*$/ and next;
        my ($sample_id, $tetrad_id, $spore_id, $PE_reads_files, $cross_pair, $note) = split /\s+/, $_;
        if (exists $sample_table{$tetrad_id}{'spore_index'}) {
            push @{$sample_table{$tetrad_id}{'spore_index'}}, $spore_id;
        } else {
            @{$sample_table{$tetrad_id}{'spore_index'}} = ($spore_id);
        }
        $sample_table{$tetrad_id}{'cross_pair'} = $cross_pair;
        $sample_table{$tetrad_id}{'note'} = $note;
    }
    return %sample_table;
}

sub parse_list_file {
    my $fh = shift @_;
    my @chr = ();
    while (<$fh>) {
        chomp;
        /^#/ and next;
	/^\s*$/ and next;
	push @chr, $_;
    }
    return @chr;
}

sub get_markers {
    my ($genome1_tag, $genome2_tag, $marker_dir, $marker_type, $markers_hashref) = @_;
    my @marker_type = ();
    if ($marker_type eq "SNP") {
	push @marker_type, "SNP";
    } elsif ($marker_type eq "INDEL") {
	push @marker_type, "INDEL";
    } elsif ($marker_type eq "BOTH") {
	push @marker_type, ("SNP", "INDEL");
    } else {
	die "Unrecognized marker type: $marker_type! Please input 'SNP' or 'INDEL' or'BOTH'. Exit!\n";
    }
    foreach my $t (@marker_type) {
	# print "extract segregating $t markers ...\n";
	my $ref_based_markers = "$marker_dir/${genome1_tag}-${genome2_tag}.ref.final.$t.markers.txt.gz";
	print "ref_based_markers: $ref_based_markers\n";
	unless (-e $ref_based_markers) {
	    die "Cannot find the marker files: $ref_based_markers! Exit!\n";
	}
	my $ref_based_markers_fh = read_file($ref_based_markers);
	parse_markers_table_file($ref_based_markers_fh, $genome1_tag, $genome2_tag, $markers_hashref);
	
	close $ref_based_markers_fh;
    }
}


sub parse_markers_table_file {
    my ($fh, $ref, $query, $markers_hashref) = @_;
    my $marker_count = 0;
    while (<$fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	/^chr\tstart/ and next;
	/^ref_chr\tref_start/ and next;
	/^parent1_chr\tparent1_start/ and next;
	/^parent2_chr\tparent2_start/ and next;
	my ($ref_chr, $ref_start, $ref_end, $ref_allele, $query_allele, $query_chr, $query_start, $query_end, $match_orientation) = split /\t/, $_;
	$$markers_hashref{'ref'}{$ref_chr}{$ref_start}{'ref_chr'} = $ref_chr;
	$$markers_hashref{'ref'}{$ref_chr}{$ref_start}{'ref_start'} = $ref_start;
	$$markers_hashref{'ref'}{$ref_chr}{$ref_start}{'ref_end'} = $ref_end;
	$$markers_hashref{'ref'}{$ref_chr}{$ref_start}{'ref_allele'} = $ref_allele;
	$$markers_hashref{'ref'}{$ref_chr}{$ref_start}{'query_chr'} = $query_chr;
	$$markers_hashref{'ref'}{$ref_chr}{$ref_start}{'query_start'} = $query_start;
	$$markers_hashref{'ref'}{$ref_chr}{$ref_start}{'query_end'} = $query_end;
	$$markers_hashref{'ref'}{$ref_chr}{$ref_start}{'query_allele'} = $query_allele;
	$$markers_hashref{'ref'}{$ref_chr}{$ref_start}{'match_orientation'} = $match_orientation;
	$marker_count++
    }
    print "total_marker_count: $marker_count\n";
}


sub parse_mpileup_file {
    my ($spore_id, $spore_index, $cross_pair, $ref, $query, $fh, $allow_heteroduplex, $markers_hashref, $genotypes_hashref, $mutations_hashref, $heteroduplex_hashref) = @_;
    while (<$fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	my ($chr, $pos, $ref_allele, $cov, $reads_bases, $reads_quals) = split /\t/, $_;
	if (exists $$markers_hashref{'ref'}{$chr}{$pos}) {

	    my $marker_ref_allele = $$markers_hashref{'ref'}{$chr}{$pos}{'ref_allele'};
	    my $marker_query_allele = $$markers_hashref{'ref'}{$chr}{$pos}{'query_allele'};
	    # print "chr=$chr, pos=$pos, ref_allele=$ref_allele, marker_parent1_allele=$marker_ref_allele, marker_parent2_allele=$marker_query_allele\n";
	    my $marker_ref_allele_length = length $marker_ref_allele;
	    my $marker_query_allele_length = length $marker_query_allele;
	    my $marker_type;
	    if (($marker_ref_allele_length == 1) and ($marker_query_allele_length == 1)) {
		$marker_type = "SNP";
	    } elsif ($marker_ref_allele_length > $marker_query_allele_length) {
		$marker_type = "DEL";
	    } else {
		$marker_type = "INS";
	    }

	    if ($cov > 0) {
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
		for (my $i = 1; $i <= $cov; $i++) {
		    my $base = substr $reads_bases, $i - 1, 1;
		    my $qual = substr $reads_quals, $i - 1, 1;
		    $qual = ord($qual) - 33;
		    # set safe bound for $qual
		    if ($qual > 40) {
			$qual = 40;
		    } elsif ($qual < 0) {
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
		my $basecall_major = $basecall_sorted[0];
		# print "basecall_major = $basecall_major\n";
		my $basecall_major_count = $basecall{$basecall_major}{'count'};

		my $basecall_secondary_major;
		my $basecall_secondary_major_count = 0;
		if ((scalar @basecall_sorted) > 1) {
		    $basecall_secondary_major = $basecall_sorted[1];
		    $basecall_secondary_major_count = $basecall{$basecall_secondary_major}{'count'};
		}
		my $basecall_major_purity = $basecall_major_count/$cov;
		my $basecall_secondary_major_purity = $basecall_secondary_major_count/$cov;

		# check for potential heteroduplex
		my $heteroduplex_flag = 0;
		if ($basecall_secondary_major_count > 0) {
		    if (($basecall_major eq $marker_ref_allele) and ($basecall_secondary_major eq $marker_query_allele)) {
			$heteroduplex_flag = 1;
		    } elsif (($basecall_major eq $marker_query_allele) and ($basecall_secondary_major eq $marker_ref_allele)) {
			$heteroduplex_flag = 1;
		    } else {
			$heteroduplex_flag = 0;
		    }
		}

		my $basecall_leveraged;
		if (($heteroduplex_flag == 1) and ($allow_heteroduplex eq "yes")) {
		    my $basecall_major_qual_sum = $basecall{$basecall_major}{'qual_sum'};
		    my $basecall_secondary_major_qual_sum = $basecall{$basecall_secondary_major}{'qual_sum'};
		    my $basecall_alt_qual_sum = 0;

		    if (scalar @basecall_sorted > 2 ) {
			for(my $b = 2; $b < (scalar @basecall_sorted); $b++) {
			    my $base = $basecall_sorted[$b];
			    $basecall_alt_qual_sum += $basecall{$base}{'qual_sum'};
			}
		    }
		    
		    my $basecall_qual_diff = $basecall_major_qual_sum + $basecall_secondary_major_qual_sum - $basecall_alt_qual_sum;
		    if (($basecall_qual_diff >= $basecall_qual_diff_cutoff) and ($basecall_major_purity + $basecall_secondary_major_purity >= $basecall_purity_cutoff)) {
			$basecall_leveraged = "heteroduplex";
		    } else {
			$basecall_leveraged = "NA";
		    }
		} else {
		    my $basecall_major_qual_sum = $basecall{$basecall_major}{'qual_sum'};
		    my $basecall_alt_qual_sum = 0;
		    if (scalar @basecall_sorted > 1 ) {
			for(my $b = 1; $b < (scalar @basecall_sorted); $b++) {
			    my $base = $basecall_sorted[$b];
			    $basecall_alt_qual_sum += $basecall{$base}{'qual_sum'};
			}
		    }
		    my $basecall_qual_diff = $basecall_major_qual_sum - $basecall_alt_qual_sum;
		    if (($basecall_qual_diff >= $basecall_qual_diff_cutoff) and ($basecall_major_purity >= $basecall_purity_cutoff)) {
			$basecall_leveraged = $basecall_major;
		    } else {
			$basecall_leveraged = "NA";
		    }

		}


		# print "basecall_major=$basecall_major, basecall_leveraged=$basecall_leveraged, basecall_qual_diff=$basecall_qual_diff, basecall_major_purity=$basecall_major_purity\n\n";
		if ($basecall_leveraged eq $marker_ref_allele) {
		    $$genotypes_hashref{'ref'}{$chr}{$pos}{$spore_index} = $ref;
		} elsif ($basecall_leveraged eq $marker_query_allele) {
                    $$genotypes_hashref{'ref'}{$chr}{$pos}{$spore_index} = $query;
		} elsif ($basecall_leveraged eq 'NA') {
		    $$genotypes_hashref{'ref'}{$chr}{$pos}{$spore_index} = "NA";
		} elsif ($basecall_leveraged eq 'heteroduplex') {
		    $$genotypes_hashref{'ref'}{$chr}{$pos}{$spore_index} = "heteroduplex";
		    $$heteroduplex_hashref{'ref'}{$chr}{$pos}{$spore_index}{'spore_id'} = $spore_id;
		    $$heteroduplex_hashref{'ref'}{$chr}{$pos}{$spore_index}{'cross_pair'} = $cross_pair;
		    $$heteroduplex_hashref{'ref'}{$chr}{$pos}{$spore_index}{'ref_allele'} = $marker_ref_allele;
		    $$heteroduplex_hashref{'ref'}{$chr}{$pos}{$spore_index}{'query_allele'} = $marker_query_allele;
		    $$heteroduplex_hashref{'ref'}{$chr}{$pos}{$spore_index}{'spore_allele'} = "heteroduplex"; #"$basecall_major,$basecall_secondary_major";
		} else {
		    $$genotypes_hashref{'ref'}{$chr}{$pos}{$spore_index} = "NA";
		    $$mutations_hashref{'ref'}{$chr}{$pos}{$spore_index}{'spore_id'} = $spore_id;
		    $$mutations_hashref{'ref'}{$chr}{$pos}{$spore_index}{'cross_pair'} = $cross_pair;
		    $$mutations_hashref{'ref'}{$chr}{$pos}{$spore_index}{'ref_allele'} = $marker_ref_allele;
		    $$mutations_hashref{'ref'}{$chr}{$pos}{$spore_index}{'query_allele'} = $marker_query_allele;
		    $$mutations_hashref{'ref'}{$chr}{$pos}{$spore_index}{'spore_allele'} = $basecall_leveraged;
		}
		# print "genotype = $$genotypes_hashref{'ref'}{$chr}{$pos}{$spore_index}\n\n";
	    } else {
		$$genotypes_hashref{'ref'}{$chr}{$pos}{$spore_index} = "NA";
	    }
	}
    }

    foreach my $chr (sort keys %{$$markers_hashref{'ref'}}) {
	foreach my $pos (sort {$a <=> $b} keys %{$$markers_hashref{'ref'}{$chr}}) {
	    if (not exists $$genotypes_hashref{'ref'}{$chr}{$pos}{$spore_index}) {
		# print "undifined markers in genotype reading!\n";
		# print "chr=$chr, pos=$pos, spore_index=$spore_index\n";
		$$genotypes_hashref{'ref'}{$chr}{$pos}{$spore_index} = "NA";
	    }
	}
    }
}


sub filter_genotype_by_cnv {
    my ($spore_id, $spore_index, $cross_pair, $genome1_tag, $genome2_tag, $fh, $genotypes_hashref, $mutations_hashref, $heteroduplex_hashref) = @_;
    my %cnv = ();
    while (<$fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	/^chr\tstart\tend\tcopy_number\tstatus\tMWU_test_p_value\tKS_test_p_value/ and next;
	my ($chr, $start, $end, $copy_number, $status, $MWU_test_p_value, $KS_test_p_value) = split /\t/, $_;
	$cnv{'ref'}{$chr}{$start}{'end'} = $end;
	$cnv{'ref'}{$chr}{$start}{'copy_number'} = $copy_number;
	$cnv{'ref'}{$chr}{$start}{'status'} = $status;
    }
    foreach my $chr (sort keys %{$cnv{'ref'}}) {
	foreach my $pos (sort {$a <=> $b} keys %{$$genotypes_hashref{'ref'}{$chr}}) {
	    foreach my $s (sort {$a <=> $b} keys %{$cnv{'ref'}{$chr}}) {
		my $e = $cnv{'ref'}{$chr}{$s}{'end'};
		if (($pos >= $s) and ($pos <= $e)) {
		    # print "convert the genotype of $spore_id at $chr:$pos from $$genotypes_hashref{'ref'}{$chr}{$pos}{$spore_index} to 'NA' due to CNV filter: $chr:$s-$e\n"; 
		    $$genotypes_hashref{'ref'}{$chr}{$pos}{$spore_index} = "NA";
		    if (exists $$mutations_hashref{'ref'}{$chr}{$pos}{$spore_index}) {
			delete $$mutations_hashref{'ref'}{$chr}{$pos}{$spore_index};
		    }
		    if (exists $$heteroduplex_hashref{'ref'}{$chr}{$pos}{$spore_index}) {
			delete $$heteroduplex_hashref{'ref'}{$chr}{$pos}{$spore_index};
		    }
		    last;
		} elsif ($pos < $s) {
		    last;
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

sub infer_genotypes {
    my ($genome1_tag, $genome2_tag, $gt_raw_hashref) = @_;
    my %allele_count = (
	'allele_parent1' => 0,
	'allele_parent2' => 0,
	'allele_na' => 0,
	'allele_mut' => 0,
	);
    foreach my $spore_index (sort keys %$gt_raw_hashref) {
	if ($$gt_raw_hashref{$spore_index} eq $genome1_tag) {
	    $allele_count{'allele_parent1'}++;
	} elsif ($$gt_raw_hashref{$spore_index} eq $genome2_tag) {
	    $allele_count{'allele_parent2'}++;
	} elsif ($$gt_raw_hashref{$spore_index} eq "NA") {
	    $allele_count{'allele_na'}++;
	} else {
	    $allele_count{'allele_mut'}++;
	}
    }
    # begin inference
    my %gt_inferred = %$gt_raw_hashref;  
    if (($allele_count{'allele_na'} == 1) and ($allele_count{'allele_mut'} == 0)) {
	foreach my $spore_index (sort keys %$gt_raw_hashref) {
	    if ($$gt_raw_hashref{$spore_index} eq "NA") {
		if (($allele_count{'allele_parent1'} == 1) and ($allele_count{'allele_parent2'} == 2)) {
		    $gt_inferred{$spore_index} = $genome1_tag;
		} elsif (($allele_count{'allele_parent2'} == 1) and ($allele_count{'allele_parent1'} == 2)) {
		    $gt_inferred{$spore_index} = $genome2_tag;
		} else {
		    # pass, no inferrence
		}
	    }
	}
    } elsif (($allele_count{'allele_na'} == 2) and ($allele_count{'allele_mut'} == 0)) {
	foreach my $spore_index (sort keys %$gt_raw_hashref) {
	    if ($$gt_raw_hashref{$spore_index} eq "NA") {
		if ($allele_count{'allele_parent1'} == 2) {
		    $gt_inferred{$spore_index} = $genome2_tag;
		} elsif ($allele_count{'allele_parent2'} == 2) {
		    $gt_inferred{$spore_index} = $genome1_tag;
		} 
	    }
	}
    }
    return %gt_inferred;
}

sub output_genotypes {
    my ($prefix, $tetrad_id, $cross_pair, $ref, $query, $chr_arrayref, $genotypes_hashref, $mutations_hashref, $heteroduplex_hashref) = @_;
    my $ref_based_genotypes_detailed_output = "$prefix.genotype.detailed.txt.gz";
    my $ref_based_genotypes_detailed_output_fh = write_file($ref_based_genotypes_detailed_output);
    print $ref_based_genotypes_detailed_output_fh "tetrad_id\tcross_pair\tref\tchr\tpos\ta_raw_GT\tb_raw_GT\tc_raw_GT\td_raw_GT\ta_inferred_GT\tb_inferred_GT\tc_inferred_GT\td_inferred_GT\n";

    my $ref_based_genotypes_lite_raw_output = "$prefix.genotype.lite.raw.txt.gz";
    my $ref_based_genotypes_lite_raw_output_fh = write_file($ref_based_genotypes_lite_raw_output);
    my $ref_based_genotypes_lite_inferred_output = "$prefix.genotype.lite.inferred.txt.gz";
    my $ref_based_genotypes_lite_inferred_output_fh = write_file($ref_based_genotypes_lite_inferred_output);
    # my $ref_based_mutation_output = "$prefix.mutations.txt.gz";
    # my $ref_based_mutation_output_fh = write_file($ref_based_mutation_output);
    # print $ref_based_mutation_output_fh "spore_id\tcross_pair\tref\tchr\tpos\ttype\tref_allele\tquery_allele\tspore_allele\n";
    my $ref_based_heteroduplex_output = "$prefix.heteroduplex.txt.gz";
    my $ref_based_heteroduplex_output_fh = write_file($ref_based_heteroduplex_output);
    print $ref_based_heteroduplex_output_fh "spore_id\tcross_pair\tref\tchr\tpos\ttype\tref_allele\tquery_allele\tspore_allele\n";
    
    my @spore_indices = qw(a b c d);
    foreach my $c (@$chr_arrayref) {
	my $chr = "ref_${c}";
	foreach my $pos (sort {$a <=> $b} keys %{$$genotypes_hashref{'ref'}{$chr}}) {
	    my %gt_raw = ();
	    foreach my $spore_index (@spore_indices) {
		if (exists $$genotypes_hashref{'ref'}{$chr}{$pos}{$spore_index}) {
		    $gt_raw{$spore_index} = $$genotypes_hashref{'ref'}{$chr}{$pos}{$spore_index};
		} else {
		    $gt_raw{$spore_index} = "NA";
		}
		# if (exists $$mutations_hashref{'ref'}{$chr}{$pos}{$spore_index}) {
		# 	print $ref_based_mutation_output_fh "$$mutations_hashref{'ref'}{$chr}{$pos}{$spore_index}{'spore_id'}\t";
		# 	print $ref_based_mutation_output_fh "$cross_pair\t$ref\t$chr\t$pos\t";
		# 	print $ref_based_mutation_output_fh "$$mutations_hashref{'ref'}{$chr}{$pos}{$spore_index}{'ref_allele'}\t";
		# 	print $ref_based_mutation_output_fh "$$mutations_hashref{'ref'}{$chr}{$pos}{$spore_index}{'query_allele'}\t";
		# 	print $ref_based_mutation_output_fh "$$mutations_hashref{'ref'}{$chr}{$pos}{$spore_index}{'spore_allele'}\t";
		# 	print $ref_based_mutation_output_fh "$$mutations_hashref{'ref'}{$chr}{$pos}{$spore_index}{'confidence_score'}\n";
		# }
		if (exists $$heteroduplex_hashref{'ref'}{$chr}{$pos}{$spore_index}) {
		    print $ref_based_heteroduplex_output_fh "$$heteroduplex_hashref{'ref'}{$chr}{$pos}{$spore_index}{'spore_id'}\t";
		    print $ref_based_heteroduplex_output_fh "$cross_pair\t$ref\t$chr\t$pos\t";
		    print $ref_based_heteroduplex_output_fh "$$heteroduplex_hashref{'ref'}{$chr}{$pos}{$spore_index}{'ref_allele'}\t";
		    print $ref_based_heteroduplex_output_fh "$$heteroduplex_hashref{'ref'}{$chr}{$pos}{$spore_index}{'query_allele'}\t";
		    print $ref_based_heteroduplex_output_fh "$$heteroduplex_hashref{'ref'}{$chr}{$pos}{$spore_index}{'spore_allele'}\n";
		}

	    }
	    
	    my %gt_inferred = infer_genotypes($ref, $query, \%gt_raw);
	    # generate lite output based on %gt_raw
	    output_lite($ref_based_genotypes_lite_raw_output_fh, $ref, $query, $chr, $pos, \%gt_raw);
	    # generate lite output for based on %gt_inferred
	    output_lite($ref_based_genotypes_lite_inferred_output_fh, $ref, $query, $chr, $pos, \%gt_inferred);
	    print $ref_based_genotypes_detailed_output_fh "$tetrad_id\t$cross_pair\tref\t$chr\t$pos\t";
	    
	    foreach my $spore_index (@spore_indices) {
		print $ref_based_genotypes_detailed_output_fh "$gt_raw{$spore_index}\t";
	    }
	    
	    foreach my $spore_index (@spore_indices) {
		print $ref_based_genotypes_detailed_output_fh "$gt_inferred{$spore_index}\t";
	    }
	    print $ref_based_genotypes_detailed_output_fh "\n";
	}
    }
}
    

sub output_lite {
    my ($fh, $ref, $query, $chr, $pos, $gt_hashref) = @_;
    my @spore_indices = qw(a b c d);
    my %gt_lite = ();
    my $parental_allele_count = 0;
    foreach my $spore_index (@spore_indices) {
	if ($$gt_hashref{$spore_index} eq $ref) {
	    $gt_lite{$spore_index} = $ref;
	    $parental_allele_count++;
	} elsif ($$gt_hashref{$spore_index} eq $query) {
            $gt_lite{$spore_index} = $query;
	    $parental_allele_count++;
	} else {
	    $gt_lite{$spore_index} = $$gt_hashref{$spore_index};
	}
    }
    print $fh "$chr\t$pos\tref";
    foreach my $spore_index (@spore_indices) {
	print $fh "\t$gt_lite{$spore_index}";
    }
    print $fh "\n";
}


