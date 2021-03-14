#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Cwd;

##############################################################
#  script: batch_tetrad_genotying_by_parent_genomes.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2019.09.09
#  description: run batch tetrad-based genotyping based on the mpileup files of spores
#  example: perl batch_tetrad_genotyping_by_parent_genomes.pl -s Master_Sample_Table.txt -m SNP -c ./../../data/chr_list.yeast.txt -q 30 -marker_dir ./../14.Polymorphic_Markers_by_Consensus -b batch_id -gamete_read_mapping_dir ./../15.Gamete_Read_Mapping_to_Parent_Genomes -output_dir output_dir
##############################################################

my ($sample_table, $marker_type, $chr_list, $basecall_qual_diff_cutoff, $basecall_purity_cutoff, $marker_dir, $batch_id, $gamete_read_mapping_dir, $output_dir);
$batch_id = "Batch_TEST";
$marker_type = "SNP";
$basecall_qual_diff_cutoff = 30;
$basecall_purity_cutoff = 0.9;
$output_dir = "$batch_id";
$marker_dir = "./../14.Polymorphic_Markers_by_Consensus";
$gamete_read_mapping_dir = "./../15.Gamete_Read_Mapping_to_Parent_Genomes";

GetOptions('sample_table|s:s' => \$sample_table, # master sample list
           'marker_type|m:s' => \$marker_type, # SNP or INDEL or BOTH
	   'qual|q:i' => \$basecall_qual_diff_cutoff, # 30
	   'purity|p:f' => \$basecall_purity_cutoff, # 30
	   'chr_list|c:s' => \$chr_list,
	   'marker_dir|m_dir:s' => \$marker_dir,
	   'batch_id|b:s' => \$batch_id, 
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
    # my @spore_indices = qw(a b c d);
    my @spore_indices = @{$tetrads{$tetrad_id}{'spore_index'}};

    # raw genotype calling
    print "raw genotype calling ...\n";
    foreach my $spore_index (@spore_indices) {
	my $spore_id = "$tetrad_id.$spore_index";
	my $genome1_based_spore_dir = "${gamete_read_mapping_dir}/${batch_id}/${cross_pair}.${spore_id}.${genome1_tag}";
	my $genome1_based_spore_mpileup = "${genome1_based_spore_dir}/${cross_pair}.${spore_id}.${genome1_tag}.mpileup.gz";
	if (-f $genome1_based_spore_mpileup) {
	    my $genome1_based_spore_mpileup_fh = read_file($genome1_based_spore_mpileup);
	    parse_mpileup_file($spore_id, $spore_index, $cross_pair, $genome1_tag, $genome2_tag, $genome1_based_spore_mpileup_fh, \%markers, \%genotypes, \%mutations);
	    my $genome1_based_spore_cnv = "${genome1_based_spore_dir}/${cross_pair}.${spore_id}.${genome1_tag}.CNV_significance_test.txt";
            my $genome1_based_spore_cnv_fh = read_file($genome1_based_spore_cnv);
            filter_genotype_by_cnv($spore_id, $spore_index, $cross_pair, $genome1_tag, $genome2_tag, $genome1_based_spore_cnv_fh, \%genotypes, \%mutations);
	}
	my $genome2_based_spore_dir = "${gamete_read_mapping_dir}/${batch_id}/${cross_pair}.${spore_id}.${genome2_tag}";
	my $genome2_based_spore_mpileup = "${genome2_based_spore_dir}/${cross_pair}.${spore_id}.${genome2_tag}.mpileup.gz";
	if (-f $genome2_based_spore_mpileup) {
	    my $genome2_based_spore_mpileup_fh = read_file($genome2_based_spore_mpileup);
	    parse_mpileup_file($spore_id, $spore_index, $cross_pair, $genome2_tag, $genome1_tag, $genome2_based_spore_mpileup_fh, \%markers, \%genotypes, \%mutations);
	    my $genome2_based_spore_cnv = "${genome2_based_spore_dir}/${cross_pair}.${spore_id}.${genome2_tag}.CNV_significance_test.txt";
            my $genome2_based_spore_cnv_fh = read_file($genome2_based_spore_cnv);
            filter_genotype_by_cnv($spore_id, $spore_index, $cross_pair, $genome2_tag, $genome1_tag, $genome2_based_spore_cnv_fh, \%genotypes, \%mutations);
	}
    }

    # cross check between ${genome1_tag}-${genome2_tag}.${genome1_tag} based signals and ${genome1_tag}-${genome2_tag}.${genome2_tag} based signals
    crosscheck_genotypes($tetrad_id, $genome1_tag, $genome2_tag, \%markers, \%genotypes);
    crosscheck_mutations($tetrad_id, $genome1_tag, $genome2_tag, \%markers, \%mutations);

    my $prefix;
    # output_by_parent1
    print "output genotypes based parent1: $genome1_tag\n";
    $prefix = "$base_dir/$output_dir/$cross_pair.$tetrad_id.$genome1_tag.q${basecall_qual_diff_cutoff}";
    output_genotypes($prefix, $tetrad_id, $cross_pair, $genome1_tag, $genome2_tag, \@chr, \%genotypes, \%mutations);
    
    # output_by_parent2
    print "output genotypes based parent2: $genome2_tag\n";
    $prefix = "$base_dir/$output_dir/$cross_pair.$tetrad_id.$genome2_tag.q${basecall_qual_diff_cutoff}";
    output_genotypes($prefix, $tetrad_id, $cross_pair, $genome2_tag, $genome1_tag, \@chr, \%genotypes, \%mutations);
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
	my $genome1_based_markers = "$marker_dir/${genome1_tag}-${genome2_tag}.$genome1_tag.final.$t.markers.txt.gz";
	my $genome2_based_markers = "$marker_dir/${genome1_tag}-${genome2_tag}.$genome2_tag.final.$t.markers.txt.gz";

	unless ((-e $genome1_based_markers) and (-e $genome2_based_markers)) {
	    die "Cannot find the marker files: $genome1_based_markers and $genome2_based_markers! Exit!\n";
	}
	my $genome1_based_markers_fh = read_file($genome1_based_markers);
	parse_markers_table_file($genome1_based_markers_fh, $genome1_tag, $genome2_tag, $markers_hashref);
	close $genome1_based_markers_fh;
	my $genome2_based_markers_fh = read_file($genome2_based_markers);
	parse_markers_table_file($genome2_based_markers_fh, $genome2_tag, $genome1_tag, $markers_hashref);
	close $genome2_based_markers_fh;
    }
}


sub parse_markers_table_file {
    my ($fh, $ref, $query, $markers_hashref) = @_;
    while (<$fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	/^chr\tstart/ and next;
	/^ref_chr\tref_start/ and next;
	/^parent1_chr\tparent1_start/ and next;
	/^parent2_chr\tparent2_start/ and next;
	my ($ref_chr, $ref_start, $ref_end, $ref_allele, $query_allele, $query_chr, $query_start, $query_end, $match_orientation) = split /\t/, $_;
	$$markers_hashref{$ref}{$query}{$ref_chr}{$ref_start}{'ref_chr'} = $ref_chr;
	$$markers_hashref{$ref}{$query}{$ref_chr}{$ref_start}{'ref_start'} = $ref_start;
	$$markers_hashref{$ref}{$query}{$ref_chr}{$ref_start}{'ref_end'} = $ref_end;
	$$markers_hashref{$ref}{$query}{$ref_chr}{$ref_start}{'ref_allele'} = $ref_allele;
	$$markers_hashref{$ref}{$query}{$ref_chr}{$ref_start}{'query_chr'} = $query_chr;
	$$markers_hashref{$ref}{$query}{$ref_chr}{$ref_start}{'query_start'} = $query_start;
	$$markers_hashref{$ref}{$query}{$ref_chr}{$ref_start}{'query_end'} = $query_end;
	$$markers_hashref{$ref}{$query}{$ref_chr}{$ref_start}{'query_allele'} = $query_allele;
	$$markers_hashref{$ref}{$query}{$ref_chr}{$ref_start}{'match_orientation'} = $match_orientation;
    }
}


sub parse_mpileup_file {
    my ($spore_id, $spore_index, $cross_pair, $ref, $query, $fh, $markers_hashref, $genotypes_hashref, $mutations_hashref) = @_;
    while (<$fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	my ($chr, $pos, $ref_allele, $cov, $reads_bases, $reads_quals) = split /\t/, $_;
	if (exists $$markers_hashref{$ref}{$query}{$chr}{$pos}) {
	    my $marker_ref_allele = $$markers_hashref{$ref}{$query}{$chr}{$pos}{'ref_allele'};
	    my $marker_query_allele = $$markers_hashref{$ref}{$query}{$chr}{$pos}{'query_allele'};
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
		    if (($base eq '.') or ($base eq ',')) {
			$base = $marker_ref_allele;
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
 		my $basecall_major_purity = $basecall_major_count/$cov;

		my $basecall_major_qual_sum = $basecall{$basecall_major}{'qual_sum'};
		my $basecall_alt_qual_sum = 0;
		foreach my $base (@basecall_sorted) {
		    $basecall_alt_qual_sum += $basecall{$base}{'qual_sum'};
		}
		my $basecall_leveraged;
		my $basecall_qual_diff = $basecall_major_qual_sum - $basecall_alt_qual_sum;
		if (($basecall_qual_diff >= $basecall_qual_diff_cutoff) and ($basecall_major_purity >= 0.9)) {
		    $basecall_leveraged = $basecall_major;
		} else {
		    $basecall_leveraged = "NA";
		}	
		# print "basecall_leveraged = $basecall_leveraged, basecall_qual_diff = $basecall_qual_diff\n\n";


		if ($basecall_leveraged eq $marker_ref_allele) {
		    $$genotypes_hashref{$ref}{$query}{$chr}{$pos}{$spore_index} = $ref;
		} elsif ($basecall_leveraged eq $marker_query_allele) {
                    $$genotypes_hashref{$ref}{$query}{$chr}{$pos}{$spore_index} = $query;
		} elsif ($basecall_leveraged eq 'NA') {
		    $$genotypes_hashref{$ref}{$query}{$chr}{$pos}{$spore_index} = "NA";
		} else {
		    $$genotypes_hashref{$ref}{$query}{$chr}{$pos}{$spore_index} = "NA";
		    $$mutations_hashref{$ref}{$query}{$chr}{$pos}{$spore_index}{'spore_id'} = $spore_id;
		    $$mutations_hashref{$ref}{$query}{$chr}{$pos}{$spore_index}{'cross_pair'} = $cross_pair;
		    $$mutations_hashref{$ref}{$query}{$chr}{$pos}{$spore_index}{'ref_allele'} = $marker_ref_allele;
		    $$mutations_hashref{$ref}{$query}{$chr}{$pos}{$spore_index}{'query_allele'} = $marker_query_allele;
		    $$mutations_hashref{$ref}{$query}{$chr}{$pos}{$spore_index}{'spore_allele'} = $basecall_leveraged;
		}
		# print "genotype = $$genotypes_hashref{$ref}{$query}{$chr}{$pos}{$spore_index}\n\n";
		#sleep(2);
	    } else {
		$$genotypes_hashref{$ref}{$query}{$chr}{$pos}{$spore_index} = "NA";
	    }
	}
    }

    foreach my $chr (sort keys %{$$markers_hashref{$ref}{$query}}) {
        foreach my $pos (sort {$a <=> $b} keys %{$$markers_hashref{$ref}{$query}{$chr}}) {
            if (not exists $$genotypes_hashref{$ref}{$query}{$chr}{$pos}{$spore_index}) {
                # print "undifined markers in genotype reading!\n";                                                                                                    
                # print "chr=$chr, pos=$pos, spore_index=$spore_index\n";                                                                                              
                $$genotypes_hashref{$ref}{$query}{$chr}{$pos}{$spore_index} = "NA";
            }
        }
    }

}

sub filter_genotype_by_cnv {
    my ($spore_id, $spore_index, $cross_pair, $ref, $query, $fh, $genotypes_hashref, $mutations_hashref) = @_;
    my %cnv = ();
    while (<$fh>) {
        chomp;
        /^#/ and next;
	/^\s*$/ and next;
        /^chr\tstart\tend\tcopy_number\tstatus\tMWU_test_p_value\tKS_test_p_value/ and next;
        my ($chr, $start, $end, $copy_number, $status, $MWU_test_p_value, $KS_test_p_value) = split /\t/, $_;
        $cnv{$ref}{$chr}{$start}{'end'} = $end;
        $cnv{$ref}{$chr}{$start}{'copy_number'} = $copy_number;
        $cnv{$ref}{$chr}{$start}{'status'} = $status;
    }
    foreach my $chr (sort keys %{$cnv{'ref'}}) {
        foreach my $pos (sort {$a <=> $b} keys %{$$genotypes_hashref{'ref'}{$chr}}) {
            foreach my $s (sort {$a <=> $b} keys %{$cnv{$ref}{$chr}}) {
		my $e = $cnv{$ref}{$chr}{$s}{'end'};
		if (($pos >= $s) and ($pos <= $e)) {
                    # print "convert the genotype of $spore_id at $chr:$pos to 'NA' due to CNV filter: $chr:$s-$e\n";
                    $$genotypes_hashref{$ref}{$chr}{$pos}{$spore_index} = "NA";
                    if (exists $$mutations_hashref{$ref}{$chr}{$pos}{$spore_index}{'spore_allele'}) {
                        $$mutations_hashref{$ref}{$chr}{$pos}{$spore_index}{'spore_allele'}= "NA";
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


sub crosscheck_genotypes {
    my ($tetrad_id, $ref, $query, $markers_hashref, $genotypes_hashref) = @_;
    my @spore_indices = qw(a b c d);
    my %markers_stat = ();
    foreach my $spore_index (@spore_indices) {
	$markers_stat{$spore_index}{'total_count'} = 0;
	$markers_stat{$spore_index}{'total_cat1_count'} = 0; # genotype resolvable for both parents (i.e. no "NA" calls)
	$markers_stat{$spore_index}{'total_cat1_culled'} = 0; # genotype resolvable for both parents but with conflicting calls 
	$markers_stat{$spore_index}{'total_cat2_count'} = 0; # only parent1 give resolvable genotype
	$markers_stat{$spore_index}{'total_cat3_count'} = 0; # only parent2 give resolvable genotype
	$markers_stat{$spore_index}{'total_cat4_count'} = 0; # genotype is not resolvale for both parents
    }

    foreach my $chr_ref (sort keys %{$$markers_hashref{$ref}{$query}}) {
	# print "chr=$chr_ref\n";
	foreach my $pos_ref (sort {$a <=> $b} keys %{$$markers_hashref{$ref}{$query}{$chr_ref}}) {
	    my $chr_query = $$markers_hashref{$ref}{$query}{$chr_ref}{$pos_ref}{'query_chr'};
	    my $pos_query = $$markers_hashref{$ref}{$query}{$chr_ref}{$pos_ref}{'query_start'};
	    foreach my $spore_index (@spore_indices) {
		if (not exists $$genotypes_hashref{$ref}{$query}{$chr_ref}{$pos_ref}{$spore_index}) {
		    $$genotypes_hashref{$ref}{$query}{$chr_ref}{$pos_ref}{$spore_index} = "NA";
		}
		if (not exists $$genotypes_hashref{$query}{$ref}{$chr_query}{$pos_query}{$spore_index}) {
		    $$genotypes_hashref{$query}{$ref}{$chr_query}{$pos_query}{$spore_index} = "NA";
		}
		$markers_stat{$spore_index}{'total_count'}++;
		if (($$genotypes_hashref{$ref}{$query}{$chr_ref}{$pos_ref}{$spore_index} ne "NA") and ($$genotypes_hashref{$query}{$ref}{$chr_query}{$pos_query}{$spore_index} ne "NA")) {
		    $markers_stat{$spore_index}{'total_cat1_count'}++;
		    if ($$genotypes_hashref{$ref}{$query}{$chr_ref}{$pos_ref}{$spore_index} ne $$genotypes_hashref{$query}{$ref}{$chr_query}{$pos_query}{$spore_index}) {
			# print "conflicting genotypes typeI for spore $tetrad_id.$spore_index, chr_ref = $chr_ref, pos_ref = $pos_ref, chr_query = $chr_query, pos_query = $pos_query >\n";
			# print "based on $ref: $$genotypes_hashref{$ref}{$query}{$chr_ref}{$pos_ref}{$spore_index} ";
			# print "based on $query: $$genotypes_hashref{$query}{$ref}{$chr_query}{$pos_query}{$spore_index}\n";
			$$genotypes_hashref{$ref}{$query}{$chr_ref}{$pos_ref}{$spore_index} = "NA";
			$$genotypes_hashref{$query}{$ref}{$chr_query}{$pos_query}{$spore_index} = "NA";
			$markers_stat{$spore_index}{'total_cat1_culled'}++;
		    }
		} elsif (($$genotypes_hashref{$ref}{$query}{$chr_ref}{$pos_ref}{$spore_index} ne "NA") and ($$genotypes_hashref{$query}{$ref}{$chr_query}{$pos_query}{$spore_index} eq "NA")) {
		    $$genotypes_hashref{$ref}{$query}{$chr_ref}{$pos_ref}{$spore_index} = "NA";
		    $$genotypes_hashref{$query}{$ref}{$chr_query}{$pos_query}{$spore_index} = "NA";
		    $markers_stat{$spore_index}{'total_cat2_count'}++;
		} elsif (($$genotypes_hashref{$ref}{$query}{$chr_ref}{$pos_ref}{$spore_index} eq "NA") and ($$genotypes_hashref{$query}{$ref}{$chr_query}{$pos_query}{$spore_index} ne "NA")) {
		    $$genotypes_hashref{$ref}{$query}{$chr_ref}{$pos_ref}{$spore_index} = "NA";
		    $$genotypes_hashref{$query}{$ref}{$chr_query}{$pos_query}{$spore_index} = "NA";
		    $markers_stat{$spore_index}{'total_cat3_count'}++;
		} else {
		    $markers_stat{$spore_index}{'total_cat4_count'}++;
		}
	    }
	}
    }
    foreach my $spore_index (@spore_indices) {
	my $retained_markers_count = $markers_stat{$spore_index}{'total_cat1_count'} - $markers_stat{$spore_index}{'total_cat1_culled'};
	my $retaining_pct = (100 * $retained_markers_count)/$markers_stat{$spore_index}{'total_count'};
	$retaining_pct = sprintf("%.2f", $retaining_pct);
	print "crosscheck_genotypes: spore_id = $tetrad_id.$spore_index, ";
	print "total_markers = $markers_stat{$spore_index}{'total_count'}, total_retained_markers = $retained_markers_count, retaining_percentage = $retaining_pct\%\n";
	print "total_cat1_markers = $markers_stat{$spore_index}{'total_cat1_count'}, total_cat1_culled = $markers_stat{$spore_index}{'total_cat1_culled'}, total_cat2_markers = $markers_stat{$spore_index}{'total_cat2_count'}, total_cat3_markers = $markers_stat{$spore_index}{'total_cat3_count'}, total_cat4_markers = $markers_stat{$spore_index}{'total_cat4_count'}\n\n";
    }
}

sub crosscheck_mutations {
    my ($tetrad_id, $ref, $query, $markers_hashref, $mutations_hashref) = @_;
    my @spore_indices = qw(a b c d);
    foreach my $chr_ref (sort keys %{$$markers_hashref{$ref}{$query}}) {
	foreach my $pos_ref (sort {$a <=> $b} keys %{$$markers_hashref{$ref}{$query}{$chr_ref}}) {
	    my $chr_query = $$markers_hashref{$ref}{$query}{$chr_ref}{$pos_ref}{'query_chr'};
	    my $pos_query = $$markers_hashref{$ref}{$query}{$chr_ref}{$pos_ref}{'query_start'};
	    foreach my $spore_index (@spore_indices) {	    
		if ((exists $$mutations_hashref{$ref}{$query}{$chr_ref}{$pos_ref}{$spore_index}) and (exists $$mutations_hashref{$query}{$ref}{$chr_query}{$pos_query}{$spore_index})) {
		    # if the mutation at this marker position have been infered based both parent strains in the same spore, then the mutation confidence score = 1;
		    $$mutations_hashref{$ref}{$query}{$chr_ref}{$pos_ref}{$spore_index}{'confidence_score'} = 1;
		    $$mutations_hashref{$query}{$ref}{$chr_query}{$pos_query}{$spore_index}{'confidence_score'} = 1;
		} elsif ((exists $$mutations_hashref{$ref}{$query}{$chr_ref}{$pos_ref}{$spore_index}) and (not exists $$mutations_hashref{$query}{$ref}{$chr_query}{$pos_query}{$spore_index})) {
		    $$mutations_hashref{$ref}{$query}{$chr_ref}{$pos_ref}{$spore_index}{'confidence_score'} = 0.5;
		} elsif ((not exists $$mutations_hashref{$ref}{$query}{$chr_ref}{$pos_ref}{$spore_index}) and (exists $$mutations_hashref{$query}{$ref}{$chr_query}{$pos_query}{$spore_index})) {
		    $$mutations_hashref{$query}{$ref}{$chr_query}{$pos_query}{$spore_index}{'confidence_score'} = 0.5;
		}
	    }
	}
    }
}


sub output_genotypes {
    my ($prefix, $tetrad_id, $cross_pair, $ref, $query, $chr_arrayref, $genotypes_hashref, $mutations_hashref) = @_;
    my $output_for_mutation = "no"; # by default, we turn of this output
    my $ref_based_genotypes_detailed_output = "$prefix.genotype.detailed.txt.gz";
    my $ref_based_genotypes_detailed_output_fh = write_file($ref_based_genotypes_detailed_output);
    print $ref_based_genotypes_detailed_output_fh "tetrad_id\tcross_pair\tref\tchr\tpos\ta_raw_GT\tb_raw_GT\tc_raw_GT\td_raw_GT\ta_inferred_GT\tb_inferred_GT\tc_inferred_GT\td_inferred_GT\n";

    my $ref_based_genotypes_lite_raw_output = "$prefix.genotype.lite.raw.txt.gz";
    my $ref_based_genotypes_lite_raw_output_fh = write_file($ref_based_genotypes_lite_raw_output);
    my $ref_based_genotypes_lite_inferred_output = "$prefix.genotype.lite.inferred.txt.gz";
    my $ref_based_genotypes_lite_inferred_output_fh = write_file($ref_based_genotypes_lite_inferred_output);
    # my $ref_based_mutation_output = "$prefix.mutations.txt.gz";
    # my $ref_based_mutation_output_fh = write_file($ref_based_mutation_output);
    # print $ref_based_mutation_output_fh "spore_id\tcross_pair\tref\tchr\tpos\ttype\tref_allele\tquery_allele\tspore_allele\tconfidence_score\n";
    my @spore_indices = qw(a b c d);
    foreach my $c (@$chr_arrayref) {
	my $chr = "${ref}_${c}";
	foreach my $pos (sort {$a <=> $b} keys %{$$genotypes_hashref{$ref}{$query}{$chr}}) {
	    my %gt_raw = ();
	    foreach my $spore_index (@spore_indices) {
		if (exists $$genotypes_hashref{$ref}{$query}{$chr}{$pos}{$spore_index}) {
		    $gt_raw{$spore_index} = $$genotypes_hashref{$ref}{$query}{$chr}{$pos}{$spore_index};
		} else {
		    $gt_raw{$spore_index} = "NA";
		}
		# if (exists $$mutations_hashref{$ref}{$query}{$chr}{$pos}{$spore_index}) {
		#     # set to >= 1 if want to only report mutations supported by mapping to both parents
		#     if ($$mutations_hashref{$ref}{$query}{$chr}{$pos}{$spore_index}{'confidence_score'} >= 0.5) {
		# 	print $ref_based_mutation_output_fh "$$mutations_hashref{$ref}{$query}{$chr}{$pos}{$spore_index}{'spore_id'}\t";
		# 	print $ref_based_mutation_output_fh "$cross_pair\t$ref\t$chr\t$pos\t";
		# 	print $ref_based_mutation_output_fh "$$mutations_hashref{$ref}{$query}{$chr}{$pos}{$spore_index}{'ref_allele'}\t";
		# 	print $ref_based_mutation_output_fh "$$mutations_hashref{$ref}{$query}{$chr}{$pos}{$spore_index}{'query_allele'}\t";
		# 	print $ref_based_mutation_output_fh "$$mutations_hashref{$ref}{$query}{$chr}{$pos}{$spore_index}{'spore_allele'}\t";
		# 	print $ref_based_mutation_output_fh "$$mutations_hashref{$ref}{$query}{$chr}{$pos}{$spore_index}{'confidence_score'}\n";
		#     }
		# }
	    }
	    
	    my %gt_inferred = infer_genotypes($ref, $query, \%gt_raw);
	    # generate lite output based on %gt_raw
	    output_lite($ref_based_genotypes_lite_raw_output_fh, $ref, $query, $chr, $pos, \%gt_raw);
	    # generate lite output for based on %gt_inferred
	    output_lite($ref_based_genotypes_lite_inferred_output_fh, $ref, $query, $chr, $pos, \%gt_inferred);
	    print $ref_based_genotypes_detailed_output_fh "$tetrad_id\t$cross_pair\t$ref\t$chr\t$pos\t";
	    
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
    my $parent_allele_count = 0;
    foreach my $spore_index (@spore_indices) {
	if ($$gt_hashref{$spore_index} eq $ref) {
	    $gt_lite{$spore_index} = $ref;
	    $parent_allele_count++;
	} elsif ($$gt_hashref{$spore_index} eq $query) {
            $gt_lite{$spore_index} = $query;
	    $parent_allele_count++;
	} else {
	    $gt_lite{$spore_index} = "NA";
	}
    }
    print $fh "$chr\t$pos\t$ref";
    foreach my $spore_index (@spore_indices) {
	print $fh "\t$gt_lite{$spore_index}";
    }
    print $fh "\n";
}


