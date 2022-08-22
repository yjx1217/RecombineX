#!/usr/bin/perl
use warnings FATAL => 'all';
#use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: batch_gamete_based_recombination_profiling_by_parent_genomes.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.07.26
#  description: run batch recombination profiling to characterize different crossover and gene conversion events
#  example: perl batch_gamete_based_recombination_profiling_by_parent_genomes.pl -s Master_Sample_Table.txt -n 3 -l 100 -d 5000 -q 30 -b batch_id
##############################################################

my $sample_table = "Master_Sample_Table.txt"; # master sample table
my $batch_id = "Batch_TEST";
my $genotype_dir;
my $n = 1; # minimal number of markers to be considered for linkage blocks
my $l = 10; # minimal marker-bounded block size to be considered for linkage blocks
my $d = 5000; # distance range for merging COs
my $qual_diff = 20; # net quality difference cutoff used for genotyping
my $output_dir = "$batch_id"; # the output directory
GetOptions('sample_table|s:s' => \$sample_table, # master sample table
	   'genotype_dir|g:s' => \$genotype_dir, # input directory for genotyping txt file
	   'batch_id|b:s' => \$batch_id, # batch name used for genotyping 
	   'min_marker_number|n:i' => \$n, # minimal enclosed markers to be considered for defining a linkage block, default = 1
	   'min_block_size|l:i' => \$l, # minimal linkage block size (bp) to be considered, default = 10 bp
	   'qual_diff|q:i' => \$qual_diff, # net quality difference used in tetrad genotyping, default = 20
	   'distance|d:s' => \$d, # event merging distance (bp), default = 5000
	   'output_dir|o:s' => \$output_dir); # the output directory

print "parameter setting table for this run>>\n";
print "sample_table = $sample_table\n";
print "genotype_dir = $genotype_dir\n";
print "batch_id = $batch_id\n";
print "min_marker_number = $n\n";
print "min_block_size = $l\n";
print "qual_diff = $qual_diff\n";
print "merging distance = $d\n";

my $sample_table_fh = read_file($sample_table);
my %tetrads = parse_sample_table_by_tetrad($sample_table_fh);
# filter_tetrad_with_inviable_spores(\%tetrads);
calculate_spore_availibility(\%tetrads);

foreach my $tetrad_id (sort keys %tetrads) {
    system("mkdir -p $output_dir/$tetrad_id");
    my $cross_pair = $tetrads{$tetrad_id}{'cross_pair'};
    my $spore_availibility = $tetrads{$tetrad_id}{'spore_availibility'};
    my ($parent1, $parent2) = split "-", $cross_pair;
    print "tetrad: $tetrad_id, cross_pair: $cross_pair, parent1: $parent1, parent2: $parent2, spore_availibility: $spore_availibility\n";
    my @ref = ($parent1, $parent2);
    foreach my $ref (@ref) {
	my @mode = ("raw", "inferred");
	foreach my $mode (@mode) {
	    my $genotype_input = "$genotype_dir/$batch_id/$cross_pair.$tetrad_id.$ref.q${qual_diff}.genotype.lite.$mode.txt.gz";
	    my $genotype_input_fh = read_file($genotype_input);
	    my $prefix = "$cross_pair.$tetrad_id.$ref.q${qual_diff}.genotype.lite.$mode.recombination_profile";
	    my %tetrad_genotypes = ();
	    my %preliminary_linkage_blocks = ();
	    my %markers = ();
	    parse_tetrad_genotypes_file($genotype_input_fh, $parent1, $parent2, $ref, $spore_availibility, \%tetrad_genotypes, \%preliminary_linkage_blocks, \%markers);
	    my %linkage_blocks = adjust_linkage_blocks($n, $l, \%preliminary_linkage_blocks);
	    my %filtered_linkage_blocks = filter_linkage_blocks(\%linkage_blocks);
	    my @merging_log = ();	
	    my $recombination_event_id = 0;
	    my %recombination_events = ();
	    my %CO_associated_GC = ();

	    foreach my $chr (sort keys %filtered_linkage_blocks) {
		my $block_total_count = scalar keys %{$filtered_linkage_blocks{$chr}};
		if ($block_total_count < 2) {
		    # print "no recombination on chromosome: $chr\n";
		} else {
		    my @block_id = sort {$a <=> $b} keys %{$filtered_linkage_blocks{$chr}};
		    for (my $i = 0; $i < $block_total_count - 1; $i++) {
			# slide through blocks with the same segregation pattern, two neighbouring blocks a time
			my $current_block_id = $block_id[$i];
			my $next_block_id = $block_id[$i + 1];
			my $current_block_genotype_pattern = $filtered_linkage_blocks{$chr}{$current_block_id}{'genotype_pattern'};
			my $next_block_genotype_pattern = $filtered_linkage_blocks{$chr}{$next_block_id}{'genotype_pattern'};
			my $current_block_segregation_pattern = $filtered_linkage_blocks{$chr}{$current_block_id}{'segregation_pattern'};
			my $next_block_segregation_pattern = $filtered_linkage_blocks{$chr}{$next_block_id}{'segregation_pattern'};
			my ($genotype_pattern_total_diff, $genotype_pattern_diff_by_NA, $affected_spores) = compare_genotype_patterns_with_NA($current_block_genotype_pattern, $next_block_genotype_pattern);
			my $genotype_pattern_actual_diff = $genotype_pattern_total_diff - $genotype_pattern_diff_by_NA;
			if ($genotype_pattern_actual_diff <= 2) {
			    # single CO
			    $recombination_event_id++;
			    my $event_type = "CO";
			    my $event_subtype = "Type1_CO"; # single CO without associated GC
			    @{$CO_associated_GC{$recombination_event_id}} = ();
			    my $type1_CO_block_start = $current_block_id;
			    my $type1_CO_block_end = $next_block_id;
			    my $type1_CO_affected_spores = $affected_spores;
			    record_recombination_event($recombination_event_id, $event_type, $event_subtype, $chr, $type1_CO_block_start, $type1_CO_block_end, $type1_CO_affected_spores, \%recombination_events, \%filtered_linkage_blocks);
			} else {
			    # double CO
			    $recombination_event_id++;
			    my $event_type = "CO"; 
			    my $event_subtype = "Type7_CO"; # single CO with associated GC on chromatids involved in the CO and no other GC existed
			    @{$CO_associated_GC{$recombination_event_id}} = ();
			    my $type7_CO_block_start = $current_block_id;
			    my $type7_CO_block_end = $next_block_id;
			    my $type7_CO_affected_spores = $affected_spores;
			    record_recombination_event($recombination_event_id, $event_type, $event_subtype, $chr, $type7_CO_block_start, $type7_CO_block_end, $type7_CO_affected_spores, \%recombination_events, \%filtered_linkage_blocks);
			}
		    }
		    merge_nearby_COs($d, \%recombination_events, \%filtered_linkage_blocks, \%CO_associated_GC, \@merging_log);
		}
	    }
	    # generate outputs
	    my $recombination_event_output = "$output_dir/$tetrad_id/$prefix.events.txt";
	    my $recombination_event_output_fh = write_file($recombination_event_output);
	    print $recombination_event_output_fh "tetrad_id\tevent_id\tevent_type\tevent_subtype\tchr\tmarker_pos_start\tmarker_pos_end\tadjusted_pos_start\tadjusted_pos_end\tadjusted_size\tmarker_raw_index_start\tmarker_raw_index_end\tmarker_effective_index_start\tmarker_effective_index_end\taffected_spores\n";
	    
	    my %event_type_count = ();
	    my $total_CO_count = 0;
	    my $total_GC_count = 0;
	    # foreach my $event_id (sort {(($a =~ /(\d+\.?\d+?)/)[0] || 0) <=> (($b =~ /(\d+\.?\d+?)/)[0] || 0)} keys %recombination_events) {
	    my %event_sort_by_coord = ();
	    foreach my $event_id (sort custom_sort keys %recombination_events) {
		my $type = $recombination_events{$event_id}{'type'};
		my $subtype = $recombination_events{$event_id}{'subtype'};
		my $chr = $recombination_events{$event_id}{'chr'};
		my $marker_pos_start = $recombination_events{$event_id}{'marker_pos_start'};
		my $marker_pos_end = $recombination_events{$event_id}{'marker_pos_end'};
		my $marker_raw_index_start = $recombination_events{$event_id}{'marker_raw_index_start'};
		my $marker_raw_index_end = $recombination_events{$event_id}{'marker_raw_index_end'};
		my $marker_effective_index_start = $recombination_events{$event_id}{'marker_effective_index_start'};
		my $marker_effective_index_end = $recombination_events{$event_id}{'marker_effective_index_end'};
		my $adjusted_pos_start;
		my $adjusted_pos_end;
		if ($type eq "CO") {
		    $total_CO_count++;
		    $adjusted_pos_start = ($marker_pos_start + $marker_pos_end)/2;
		    $adjusted_pos_end = ($marker_pos_start + $marker_pos_end)/2;
		} else {
		    # type is GC
		    $total_GC_count++;
		    if ((exists $markers{'effective'}{$marker_effective_index_start - 1}) and ($markers{'effective'}{$marker_effective_index_start - 1}{'chr'} eq $chr)) {
			# the previous marker is an internal marker on the same chromosome
			my $previous_marker_pos = $markers{'effective'}{$marker_effective_index_start-1}{'pos'};
			$adjusted_pos_start = ($previous_marker_pos + $marker_pos_start)/2;
		    } else {
			# the GC starting marker is the first marker on the chromosome
			$adjusted_pos_start = $marker_pos_start;
		    }
		    if ((exists $markers{'effective'}{$marker_effective_index_end + 1}) and ($markers{'effective'}{$marker_effective_index_end + 1}{'chr'} eq $chr)) {
			# the next marker is an internal marker on the same chromosome
			my $next_marker_pos = $markers{'effective'}{$marker_effective_index_end + 1}{'pos'};
			$adjusted_pos_end = ($marker_pos_end + $next_marker_pos)/2;
		    } else {
			# the GC ending marker is the last marker on the chromosome
			$adjusted_pos_end = $marker_pos_end;
		    }
		}
		my $adjusted_size = $adjusted_pos_end - $adjusted_pos_start + 1; 
		my $affected_spores = $recombination_events{$event_id}{'affected_spores'};
		$event_sort_by_coord{$chr}{$marker_pos_start}{$marker_pos_end}{$type}{$event_id} =  "$tetrad_id\t$event_id\t$type\t$subtype\t$chr\t$marker_pos_start\t$marker_pos_end\t$adjusted_pos_start\t$adjusted_pos_end\t$adjusted_size\t$marker_raw_index_start\t$marker_raw_index_end\t$marker_effective_index_start\t$marker_effective_index_end\t$affected_spores\n";
		
		if (exists $event_type_count{$type}{$subtype}) {
		    $event_type_count{$type}{$subtype}++;
		} else {
		    $event_type_count{$type}{$subtype} = 1;
		}
	    }

	    foreach my $chr (sort keys %event_sort_by_coord) {
		foreach my $s (sort {$a <=> $b} keys %{$event_sort_by_coord{$chr}}) {
		    foreach my $e (sort {$a <=> $b} keys %{$event_sort_by_coord{$chr}{$s}}) {
			foreach my $t (sort keys %{$event_sort_by_coord{$chr}{$s}{$e}}) {
			    foreach my $id (sort keys %{$event_sort_by_coord{$chr}{$s}{$e}{$t}}) {
				print $recombination_event_output_fh $event_sort_by_coord{$chr}{$s}{$e}{$t}{$id}
			    }
			}
		    }
		}
	    }

	    # event type count summary output
	    my $event_type_count_out = "$output_dir/$tetrad_id/$prefix.event_type_count.txt";
	    my $event_type_count_out_fh = write_file($event_type_count_out);
	    print $event_type_count_out_fh "total CO events: $total_CO_count\n";
	    print $event_type_count_out_fh "total GC events: $total_GC_count\n";
	    print $event_type_count_out_fh "################################\n";
	    foreach my $type (sort keys %event_type_count) {
		my %subtype_index = ();
		foreach my $subtype (sort keys %{$event_type_count{$type}}) {
		    my ($subtype_index) = ($subtype =~ /Type(\d+)_/);
		    $subtype_index{$subtype} = $subtype_index;
		}
		foreach my $subtype (sort {$subtype_index{$a} <=> $subtype_index{$b}} keys %subtype_index) {
		    print $event_type_count_out_fh "$type\t$subtype\t$event_type_count{$type}{$subtype}\n";
		}
	    }
	    print $event_type_count_out_fh "################################\n";
	    # output merging record
	    my $merging_log_out = "$output_dir/$tetrad_id/$prefix.merging_log.txt";
	    my $merging_log_out_fh = write_file($merging_log_out);
	    foreach my $line (@merging_log) {
		print $merging_log_out_fh "$line\n";
	    }
	    # output CO_associated GC output
	    my $CO_associated_GC_out = "$output_dir/$tetrad_id/$prefix.co_associated_gc.txt";
	    my $CO_associated_GC_out_fh = write_file($CO_associated_GC_out);
	    print $CO_associated_GC_out_fh "event_id_for_CO\tCO_subtype\tevent_id_for_associated_GC\tassociated_GC_subtype\n";
	    foreach my $CO_id (sort keys %CO_associated_GC) {
		my $CO_subtype = $recombination_events{$CO_id}{'subtype'};
		my @GC_id = sort @{$CO_associated_GC{$CO_id}};
		if ((scalar @GC_id) > 0) {
		    foreach my $GC_id (@GC_id) {
			# print "CO_id: $CO_id, GC_id: $GC_id\n";
			my $GC_subtype = $recombination_events{$GC_id}{'subtype'};
			print $CO_associated_GC_out_fh "$CO_id\t$CO_subtype\t$GC_id\t$GC_subtype\n";
		    }
		}
	    }
	    # output marker record
	    my $markers_out = "$output_dir/$tetrad_id/$prefix.markers.txt";
	    my $markers_out_fh = write_file($markers_out);
	    print $markers_out_fh "tetrad_id\tchr\tpos\traw_index\teffective_index\tgenotype_pattern\tsegregation_pattern\n";
	    foreach my $raw_index (sort {$a<=>$b} keys %{$markers{'raw'}}) {
		my $chr = $markers{'raw'}{$raw_index}{'chr'};
		my $pos = $markers{'raw'}{$raw_index}{'pos'};
		my $genotype_pattern = $markers{'raw'}{$raw_index}{'genotype_pattern'};
		my $segregation_pattern = $markers{'raw'}{$raw_index}{'segregation_pattern'};
		my $effective_index = $markers{'raw'}{$raw_index}{'effective_index'};
		print $markers_out_fh "$tetrad_id\t$chr\t$pos\t$raw_index\t$effective_index\t$genotype_pattern\t$segregation_pattern\n";
	    }
	    # output preliminary linkage blocks
	    my $preliminary_linkage_blocks_out = "$output_dir/$tetrad_id/$prefix.preliminary_linkage_blocks.txt";
	    my $preliminary_linkage_blocks_out_fh = write_file($preliminary_linkage_blocks_out);
	    print $preliminary_linkage_blocks_out_fh "tetrad_id\tchr:preliminary_block_id\tmarker_pos_start-marker_pos_end\tmarker_raw_index_start-marker_raw_index_end\tmarker_effective_index_start-marker_effective_index_end\tgenotype_pattern\tsegregation_pattern\n";
	    
	    foreach my $chr (sort keys %preliminary_linkage_blocks) {
		foreach my $block_id (sort {$a<=>$b} keys %{$preliminary_linkage_blocks{$chr}}) {
		    my $marker_pos_start = $preliminary_linkage_blocks{$chr}{$block_id}{'marker_pos_start'};
		    my $marker_pos_end = $preliminary_linkage_blocks{$chr}{$block_id}{'marker_pos_end'};
		    my $marker_raw_index_start = $preliminary_linkage_blocks{$chr}{$block_id}{'marker_raw_index_start'};
		    my $marker_raw_index_end = $preliminary_linkage_blocks{$chr}{$block_id}{'marker_raw_index_end'};
		    my $marker_effective_index_start = $preliminary_linkage_blocks{$chr}{$block_id}{'marker_effective_index_start'};
		    my $marker_effective_index_end = $preliminary_linkage_blocks{$chr}{$block_id}{'marker_effective_index_end'};
		    my $genotype_pattern = $preliminary_linkage_blocks{$chr}{$block_id}{'genotype_pattern'}; # genotype_a:genotype_b:genotype_c:genotype_d
		    my $segregation_pattern = $preliminary_linkage_blocks{$chr}{$block_id}{'segregation_pattern'};

		    print $preliminary_linkage_blocks_out_fh "$tetrad_id\t$chr:$block_id\t$marker_pos_start-$marker_pos_end\t$marker_raw_index_start-$marker_raw_index_end\t$marker_effective_index_start-$marker_effective_index_end\t$genotype_pattern\t$segregation_pattern\n";
		}
	    }
	    # output linkage blocks
	    my $linkage_blocks_out = "$output_dir/$tetrad_id/$prefix.linkage_blocks.txt";
	    my $linkage_blocks_out_fh = write_file($linkage_blocks_out);
	    print $linkage_blocks_out_fh "tetrad_id\tchr:block_id\tmarker_pos_start-marker_pos_end\tmarker_raw_index_start-marker_raw_index_end\tmarker_effective_index_start-marker_effective_index_end\tgenotype_pattern\tsegregation_pattern\n";
	    
	    foreach my $chr (sort keys %linkage_blocks) {
		foreach my $block_id (sort {$a<=>$b} keys %{$linkage_blocks{$chr}}) {
		    my $marker_pos_start = $linkage_blocks{$chr}{$block_id}{'marker_pos_start'};
		    my $marker_pos_end = $linkage_blocks{$chr}{$block_id}{'marker_pos_end'};
		    my $marker_raw_index_start = $linkage_blocks{$chr}{$block_id}{'marker_raw_index_start'};
		    my $marker_raw_index_end = $linkage_blocks{$chr}{$block_id}{'marker_raw_index_end'};
		    my $marker_effective_index_start = $linkage_blocks{$chr}{$block_id}{'marker_effective_index_start'};
		    my $marker_effective_index_end = $linkage_blocks{$chr}{$block_id}{'marker_effective_index_end'};
		    my $genotype_pattern = $linkage_blocks{$chr}{$block_id}{'genotype_pattern'}; # genotype_a:genotype_b:genotype_c:genotype_d
		    my $segregation_pattern = $linkage_blocks{$chr}{$block_id}{'segregation_pattern'};
		    print $linkage_blocks_out_fh "$tetrad_id\t$chr:$block_id\t$marker_pos_start-$marker_pos_end\t$marker_raw_index_start-$marker_raw_index_end\t$marker_effective_index_start-$marker_effective_index_end\t$genotype_pattern\t$segregation_pattern\n";
		}
	    }
	    # output filtered linkage blocks
	    my $filtered_linkage_blocks_out = "$output_dir/$tetrad_id/$prefix.filtered_linkage_blocks.txt";
	    my $filtered_linkage_blocks_out_fh = write_file($filtered_linkage_blocks_out);
	    print $filtered_linkage_blocks_out_fh "tetrad_id\tchr:block_id\tmarker_pos_start-marker_pos_end\tmarker_raw_index_start-marker_raw_index_end\tmarker_effective_index_start-marker_effective_index_end\tgenotype_pattern\tsegregation_pattern\n";	    
	    foreach my $chr (sort keys %filtered_linkage_blocks) {
		foreach my $block_id (sort {$a<=>$b} keys %{$filtered_linkage_blocks{$chr}}) {
		    my $marker_pos_start = $filtered_linkage_blocks{$chr}{$block_id}{'marker_pos_start'};
		    my $marker_pos_end = $filtered_linkage_blocks{$chr}{$block_id}{'marker_pos_end'};
		    my $marker_raw_index_start = $filtered_linkage_blocks{$chr}{$block_id}{'marker_raw_index_start'};
		    my $marker_raw_index_end = $filtered_linkage_blocks{$chr}{$block_id}{'marker_raw_index_end'};
		    my $marker_effective_index_start = $filtered_linkage_blocks{$chr}{$block_id}{'marker_effective_index_start'};
		    my $marker_effective_index_end = $filtered_linkage_blocks{$chr}{$block_id}{'marker_effective_index_end'};
		    my $genotype_pattern = $filtered_linkage_blocks{$chr}{$block_id}{'genotype_pattern'}; # genotype_a:genotype_b:genotype_c:genotype_d
		    my $segregation_pattern = $filtered_linkage_blocks{$chr}{$block_id}{'segregation_pattern'};
		    print $filtered_linkage_blocks_out_fh "$tetrad_id\t$chr:$block_id\t$marker_pos_start-$marker_pos_end\t$marker_raw_index_start-$marker_raw_index_end\t$marker_effective_index_start-$marker_effective_index_end\t$genotype_pattern\t$segregation_pattern\n";
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

sub calculate_spore_availibility {
    my $tetrads_hashref = shift @_;
    my @spore_index_array = qw(a b c d);
    foreach my $tetrad_id (sort keys %$tetrads_hashref) {
	my %spore_availibility_hash = (
	    'a' => 0,
	    'b' => 0,
	    'c' => 0,
	    'd' => 0
	    );
	my @spore_availibility_array = ();
	foreach my $spore_index (@{$$tetrads_hashref{$tetrad_id}{'spore_index'}}) {
	    $spore_availibility_hash{$spore_index} = 1;
	}
	foreach my $spore_index (@spore_index_array) {
	    push @spore_availibility_array, $spore_availibility_hash{$spore_index};
	}
	$$tetrads_hashref{$tetrad_id}{'spore_availibility'} = join ":", @spore_availibility_array;
    }
}

sub filter_tetrad_with_inviable_spores {
    my $tetrads_hashref = shift @_;
    foreach my $tetrad_id (sort keys %$tetrads_hashref) {
	my $viable_spore_count = scalar @{$$tetrads_hashref{$tetrad_id}{'spore_index'}};
	if ($viable_spore_count < 4) {
	    delete $$tetrads_hashref{$tetrad_id};
	}
    }
}


sub compare_genotype_patterns {
    my ($genotype_pattern1, $genotype_pattern2) = @_;
    my @genotype_pattern1 = split /:/, $genotype_pattern1;
    my @genotype_pattern2 = split /:/, $genotype_pattern2;
    my $genotype_diff_count = 0;
    my @affected_spores = (0, 0, 0, 0);
    for (my $i = 0; $i < 4; $i++) {
	if ($genotype_pattern1[$i] ne $genotype_pattern2[$i]) {
	  $genotype_diff_count++;
	  $affected_spores[$i] = 1;
	}
    }
    my $affected_spores = join ":", @affected_spores;
    return ($genotype_diff_count, $affected_spores);
}

sub compare_genotype_patterns_with_NA {
    my ($genotype_pattern1, $genotype_pattern2) = @_;
    my @genotype_pattern1 = split /:/, $genotype_pattern1;
    my @genotype_pattern2 = split /:/, $genotype_pattern2;
    my $genotype_diff_count = 0;
    my $genotype_diff_by_NA_count = 0;
    my @affected_spores = (0, 0, 0, 0);
    for (my $i = 0; $i < 4; $i++) {
	if ($genotype_pattern1[$i] ne $genotype_pattern2[$i]) {
	  $genotype_diff_count++;
	  if (($genotype_pattern1[$i] eq "NA") or ($genotype_pattern2[$i] eq "NA")) {
	      $genotype_diff_by_NA_count++;
	  } else {
	      $affected_spores[$i] = 1;
	  }
	}
    }
    my $affected_spores = join ":", @affected_spores;
    return ($genotype_diff_count, $genotype_diff_by_NA_count, $affected_spores);
}



sub parse_tetrad_genotypes_file {
    my ($fh, $parent1, $parent2, $ref, $spore_availibility, $tetrad_genotypes_hashref, $linkage_blocks_hashref, $markers_hashref) = @_;
    my $marker_raw_index = 0;
    my $marker_effective_index = 0;
    my $block_id = 0;
    my $genotype_pattern_just_saw;
    my @spore_availibility_array = split ":", $spore_availibility;
    while (<$fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	my ($chr, $pos, $reftag, $genotype_a, $genotype_b, $genotype_c, $genotype_d) = split /\t/, $_;
	my $parent1_genotype_count = 0;
	my $parent2_genotype_count = 0;
	my $NA_genotype_count = 0;
	my $missing_genotype_count = 0;
	my @genotypes = ($genotype_a, $genotype_b, $genotype_c, $genotype_d);
	for (my $i = 0; $i < (scalar @genotypes); $i++) {
	    if ($genotypes[$i] eq $parent1) {
		    $parent1_genotype_count++;
	    } elsif ($genotypes[$i] eq $parent2) {
		    $parent2_genotype_count++;
	    } else {
		$NA_genotype_count++;
		if ($spore_availibility_array[$i] == 1) {
		    $missing_genotype_count++;
		}
	    }

	}
	my $segregation_pattern = "$parent1_genotype_count:$parent2_genotype_count:$NA_genotype_count";
	my $genotype_pattern = "$genotype_a:$genotype_b:$genotype_c:$genotype_d";
	$marker_raw_index++;
	$$markers_hashref{'raw'}{$marker_raw_index}{'chr'} = $chr;
	$$markers_hashref{'raw'}{$marker_raw_index}{'pos'} = $pos;
	$$markers_hashref{'raw'}{$marker_raw_index}{'reftag'} = $reftag;
	$$markers_hashref{'raw'}{$marker_raw_index}{'genotype_pattern'} = $genotype_pattern;
	$$markers_hashref{'raw'}{$marker_raw_index}{'segregation_pattern'} = $segregation_pattern;
	$$markers_hashref{'raw'}{$marker_raw_index}{'effective_index'} = "";
	##################
	#ignore markers with missing genotypes (for viable spores) for now
	if ($missing_genotype_count > 0) {
	   next;
	} else {
	   # only consider and index non-ambiguous markers, this will help calculating inter-marker boundaries
	    $marker_effective_index++;
	    $$markers_hashref{'raw'}{$marker_raw_index}{'effective_index'} = $marker_effective_index;
	    $$markers_hashref{'effective'}{$marker_effective_index}{'chr'} = $chr;
	    $$markers_hashref{'effective'}{$marker_effective_index}{'pos'} = $pos;
	    $$markers_hashref{'effective'}{$marker_effective_index}{'reftag'} = $reftag;
	    $$markers_hashref{'effective'}{$marker_effective_index}{'genotype_pattern'} = $genotype_pattern;
	    $$markers_hashref{'effective'}{$marker_effective_index}{'segregation_pattern'} = $segregation_pattern;
	    $$markers_hashref{'effective'}{$marker_effective_index}{'raw_index'} = $marker_raw_index;
	}
	$$tetrad_genotypes_hashref{$chr}{$pos}{'marker_raw_index'} = $marker_raw_index;
	$$tetrad_genotypes_hashref{$chr}{$pos}{'marker_effective_index'} = $marker_effective_index;
	$$tetrad_genotypes_hashref{$chr}{$pos}{'genotype_a'} = $genotype_a;
	$$tetrad_genotypes_hashref{$chr}{$pos}{'genotype_b'} = $genotype_b;
	$$tetrad_genotypes_hashref{$chr}{$pos}{'genotype_c'} = $genotype_c;
	$$tetrad_genotypes_hashref{$chr}{$pos}{'genotype_d'} = $genotype_d;
	$$tetrad_genotypes_hashref{$chr}{$pos}{'segregation_pattern'} = $segregation_pattern;
	$$tetrad_genotypes_hashref{$chr}{$pos}{'genotype_pattern'} = $genotype_pattern;
	if (not exists $$linkage_blocks_hashref{$chr}) {
	    # this is the first marker on the chromosome
	    $block_id++;
	    @{$$linkage_blocks_hashref{$chr}{$block_id}{'markers'}} = ($pos);
	    $$linkage_blocks_hashref{$chr}{$block_id}{'marker_pos_start'} = $pos;
	    $$linkage_blocks_hashref{$chr}{$block_id}{'marker_pos_end'} = $pos;
	    $$linkage_blocks_hashref{$chr}{$block_id}{'marker_raw_index_start'} = $marker_raw_index;
	    $$linkage_blocks_hashref{$chr}{$block_id}{'marker_raw_index_end'} = $marker_raw_index;
	    $$linkage_blocks_hashref{$chr}{$block_id}{'marker_effective_index_start'} = $marker_effective_index;
	    $$linkage_blocks_hashref{$chr}{$block_id}{'marker_effective_index_end'} = $marker_effective_index;
	    $$linkage_blocks_hashref{$chr}{$block_id}{'segregation_pattern'} = $segregation_pattern;
	    $$linkage_blocks_hashref{$chr}{$block_id}{'genotype_pattern'} = $genotype_pattern;
	    $$linkage_blocks_hashref{$chr}{$block_id}{'genotype_switch_pattern'} = "0:0:0:0"; # genotype switches relative to the previous linkeage block
	} elsif ($genotype_pattern eq $genotype_pattern_just_saw) {
	    # markers with the same genotype pattern: genotype_a:genotype_b:genotype_c:genotype_d 
	    push @{$$linkage_blocks_hashref{$chr}{$block_id}{'markers'}}, $pos;
	    $$linkage_blocks_hashref{$chr}{$block_id}{'marker_pos_end'} = $pos;
	    $$linkage_blocks_hashref{$chr}{$block_id}{'marker_raw_index_end'} = $marker_raw_index;
	    $$linkage_blocks_hashref{$chr}{$block_id}{'marker_effective_index_end'} = $marker_effective_index;
	} else {
	    # markers with different genotype pattern: genotype_a:genotype_b:genotype_c:genotype_d 
	    $block_id++;
	    @{$$linkage_blocks_hashref{$chr}{$block_id}{'markers'}} = ($pos);
	    $$linkage_blocks_hashref{$chr}{$block_id}{'marker_pos_start'} = $pos;
	    $$linkage_blocks_hashref{$chr}{$block_id}{'marker_pos_end'} = $pos;
	    $$linkage_blocks_hashref{$chr}{$block_id}{'marker_raw_index_start'} = $marker_raw_index;
	    $$linkage_blocks_hashref{$chr}{$block_id}{'marker_raw_index_end'} = $marker_raw_index;
	    $$linkage_blocks_hashref{$chr}{$block_id}{'marker_effective_index_start'} = $marker_effective_index;
	    $$linkage_blocks_hashref{$chr}{$block_id}{'marker_effective_index_end'} = $marker_effective_index;
	    $$linkage_blocks_hashref{$chr}{$block_id}{'segregation_pattern'} = $segregation_pattern;
	    $$linkage_blocks_hashref{$chr}{$block_id}{'genotype_pattern'} = $genotype_pattern;
	    # genotype switches relative to the previous linkage block
	    my ($genotype_pattern_total_diff, $genotype_pattern_diff_by_NA, $affected_spores) = compare_genotype_patterns_with_NA($genotype_pattern_just_saw, $genotype_pattern);
 	    $$linkage_blocks_hashref{$chr}{$block_id}{'genotype_switch_pattern'} = $affected_spores;
	}
	$genotype_pattern_just_saw = $genotype_pattern;
    }
}


sub count_NA_genotype {
    my $genotype_pattern = shift @_;
    my $NA_genotype_count = 0;
    my @genotypes = split /:/, $genotype_pattern;
    foreach my $gt (@genotypes) {
	if ($gt eq "NA") {
	    $NA_genotype_count++;
	}
    }
    return $NA_genotype_count;
}

sub adjust_linkage_blocks {
    my ($n, $l, $preliminary_linkage_blocks_hashref) = @_;
    my %new_linkage_blocks = ();
    foreach my $chr (sort keys %$preliminary_linkage_blocks_hashref) {
	my $new_block_id = 0;
	foreach my $block_id (sort {$a <=> $b} keys %{$$preliminary_linkage_blocks_hashref{$chr}}) {
	    my $block_marker_count = scalar @{$$preliminary_linkage_blocks_hashref{$chr}{$block_id}{'markers'}};
	    my $min_block_size = $$preliminary_linkage_blocks_hashref{$chr}{$block_id}{'marker_pos_end'} - $$preliminary_linkage_blocks_hashref{$chr}{$block_id}{'marker_pos_start'} + 1;
	    my $marker_count_by_index = $$preliminary_linkage_blocks_hashref{$chr}{$block_id}{'marker_effective_index_end'} - $$preliminary_linkage_blocks_hashref{$chr}{$block_id}{'marker_effective_index_start'} + 1;
	    # print "checkpoint 1: chr=$chr, block_id=$block_id, block_marker_count=$block_marker_count, marker_count_by_index=$marker_count_by_index, min_block_size=$min_block_size\n";
	    # for debug
	    if (($block_marker_count >= $n) and ($min_block_size >= $l)) {
		# qualified block
		if ($new_block_id == 0) {
		    # this is the first new block on the chromosome
		    $new_block_id++;
		    %{$new_linkage_blocks{$chr}{$new_block_id}} = %{$$preliminary_linkage_blocks_hashref{$chr}{$block_id}};
		} else {
		    # check if merging with the previous qualified block is needed
		    my $current_block_genotype_pattern = $$preliminary_linkage_blocks_hashref{$chr}{$block_id}{'genotype_pattern'};
		    my $previous_new_block_genotype_pattern = $new_linkage_blocks{$chr}{$new_block_id}{'genotype_pattern'};
		    my ($genotype_diff_count, $genotype_switch_pattern) = compare_genotype_patterns($current_block_genotype_pattern, $previous_new_block_genotype_pattern);
		    if ($genotype_diff_count == 0) {
			# this block has the same genotype as the previous qualified block, merging the two blocks
			push @{$new_linkage_blocks{$chr}{$new_block_id}{'markers'}}, @{$$preliminary_linkage_blocks_hashref{$chr}{$block_id}{'markers'}};
			$new_linkage_blocks{$chr}{$new_block_id}{'marker_pos_end'} = $$preliminary_linkage_blocks_hashref{$chr}{$block_id}{'marker_pos_end'};
			$new_linkage_blocks{$chr}{$new_block_id}{'marker_raw_index_end'} = $$preliminary_linkage_blocks_hashref{$chr}{$block_id}{'marker_raw_index_end'};
			$new_linkage_blocks{$chr}{$new_block_id}{'marker_effective_index_end'} = $$preliminary_linkage_blocks_hashref{$chr}{$block_id}{'marker_effective_index_end'};
		    } else {
			# this block has different genotype as the previous qualified block, initialize a new block
			$new_block_id++;
			%{$new_linkage_blocks{$chr}{$new_block_id}} = %{$$preliminary_linkage_blocks_hashref{$chr}{$block_id}};
			$new_linkage_blocks{$chr}{$new_block_id}{'genotype_switch_pattern'} = $genotype_switch_pattern;
		    }
		}
	    }
	}
    }
    return %new_linkage_blocks;
}


sub filter_linkage_blocks {
    my ($raw_linkage_blocks_hashref) = @_;
    my %filtered_linkage_blocks = ();
    foreach my $chr (sort keys %$raw_linkage_blocks_hashref) {
	my $filtered_block_id = 0;
	foreach my $block_id (sort {$a <=> $b} keys %{$$raw_linkage_blocks_hashref{$chr}}) {
	    my $segregation_pattern = $$raw_linkage_blocks_hashref{$chr}{$block_id}{'segregation_pattern'};
	    if ($segregation_pattern !~ /(3\:1\:0|1\:3\:0|3\:0\:1|0\:3\:1|4\:0\:0|0\:4\:0)/) {
		# qualified block
		if ($filtered_block_id == 0) {
		    # this is the first new block on the chromosome
		    $filtered_block_id++;
		    %{$filtered_linkage_blocks{$chr}{$filtered_block_id}} = %{$$raw_linkage_blocks_hashref{$chr}{$block_id}};
		} else {
		    # check if merging with the previous qualified block is needed
		    my $current_block_genotype_pattern = $$raw_linkage_blocks_hashref{$chr}{$block_id}{'genotype_pattern'};
		    my $previous_filtered_block_genotype_pattern = $filtered_linkage_blocks{$chr}{$filtered_block_id}{'genotype_pattern'};
		    my ($genotype_pattern_total_diff, $genotype_pattern_diff_by_NA, $genotype_switch_pattern) = compare_genotype_patterns_with_NA($current_block_genotype_pattern, $previous_filtered_block_genotype_pattern);
		    my $genotype_pattern_actual_diff = $genotype_pattern_total_diff - $genotype_pattern_diff_by_NA;
		    if ($genotype_pattern_actual_diff == 0) {
			# this block has the same genotype as the previous qualified block, merging the two blocks
			push @{$filtered_linkage_blocks{$chr}{$filtered_block_id}{'markers'}}, @{$$raw_linkage_blocks_hashref{$chr}{$block_id}{'markers'}};
			$filtered_linkage_blocks{$chr}{$filtered_block_id}{'marker_pos_end'} = $$raw_linkage_blocks_hashref{$chr}{$block_id}{'marker_pos_end'};
			$filtered_linkage_blocks{$chr}{$filtered_block_id}{'marker_raw_index_end'} = $$raw_linkage_blocks_hashref{$chr}{$block_id}{'marker_raw_index_end'};
			$filtered_linkage_blocks{$chr}{$filtered_block_id}{'marker_effective_index_end'} = $$raw_linkage_blocks_hashref{$chr}{$block_id}{'marker_effective_index_end'};
		    } else {
			# this block has different genotype as the previous qualified block, initialize a new block
			$filtered_block_id++;
			%{$filtered_linkage_blocks{$chr}{$filtered_block_id}} = %{$$raw_linkage_blocks_hashref{$chr}{$block_id}};
			$filtered_linkage_blocks{$chr}{$filtered_block_id}{'genotype_switch_pattern'} = $genotype_switch_pattern;
		    }
		}
	    }
	}
    }
    return %filtered_linkage_blocks;
}

sub check_genotype_switches {
    # this function will also consider the interstitial blocks between the two input blocks if there is any
    my ($block_start, $block_end, $linkage_blocks_on_chr_hashref) = @_;
    my @spores = qw(a b c d);
    my %genotype_switches = ();
    foreach my $spore (@spores) {
	$genotype_switches{$spore}{'count'} = 0;
	@{$genotype_switches{$spore}{'block_boundary'}} = ();
    }
    for (my $i = $block_start; $i < $block_end; $i++) {
	my @current_block_genotypes = split /:/, $$linkage_blocks_on_chr_hashref{$i}{'genotype_pattern'};
	my @next_block_genotypes = split /:/, $$linkage_blocks_on_chr_hashref{$i + 1}{'genotype_pattern'};
	for (my $j = 0; $j < scalar @spores; $j++) {
	    if ($current_block_genotypes[$j] ne $next_block_genotypes[$j]) {
		$genotype_switches{$spores[$j]}{'count'}++;
		my $next_block_id = $i + 1;
		push @{$genotype_switches{$spores[$j]}{'block_boundary'}}, "$i|$next_block_id";
	    }
	}
    }
    my @genotype_switches_count = ();
    my %uniq_single_switches = (); # a non-redundant set of block boundaries with a single spore-specific genotype switch
    my @genotype_switch_boundaries = ();
    foreach my $spore (@spores) {
	push @genotype_switches_count,  $genotype_switches{$spore}{'count'};
	my $spore_genotype_switch_boundaries = join ",", @{$genotype_switches{$spore}{'block_boundary'}};
	push @genotype_switch_boundaries, $spore_genotype_switch_boundaries;
	if ($genotype_switches{$spore}{'count'} == 1) {
	    # single switch
	    my $switch_boundary = ${$genotype_switches{$spore}{'block_boundary'}}[0];
	    if (not exists $uniq_single_switches{$switch_boundary}) {
		$uniq_single_switches{$switch_boundary} = 1;
	    } else {
		$uniq_single_switches{$switch_boundary}++;
	    }
	}
    }
    
    my $genotype_switches_count_pattern = join ":", @genotype_switches_count;
    my @genotype_switches_count_sorted = sort {$a <=> $b} @genotype_switches_count;
    my $genotype_switches_count_pattern_sorted = join ":", @genotype_switches_count_sorted;
    my $uniq_single_switch_count = scalar keys %uniq_single_switches;
    
    my $genotype_switch_boundaries = join ":", @genotype_switch_boundaries;

    return ($genotype_switches_count_pattern, $genotype_switches_count_pattern_sorted, $uniq_single_switch_count, $genotype_switch_boundaries);
}


sub record_recombination_event {
    my ($id, $type, $subtype, $chr, $block_start, $block_end, $affected_spores, $events_hashref, $linkage_blocks_hashref) = @_;
    # print "record_event: id=$id, type=$type, subtype=$subtype, chr=$chr, block_start=$block_start, block_end=$block_end, affected_spores=$affected_spores\n";
    $$events_hashref{$id}{'type'} = $type;
    $$events_hashref{$id}{'subtype'} = $subtype;
    $$events_hashref{$id}{'chr'} = $chr;
    $$events_hashref{$id}{'block_start'} = $block_start;
    $$events_hashref{$id}{'block_end'} = $block_end;
    if ($type eq "CO") {
	# if event is CO
	$$events_hashref{$id}{'marker_pos_start'} = $$linkage_blocks_hashref{$chr}{$block_start}{'marker_pos_end'};
	$$events_hashref{$id}{'marker_pos_end'} = $$linkage_blocks_hashref{$chr}{$block_end}{'marker_pos_start'};
	$$events_hashref{$id}{'marker_raw_index_start'} = $$linkage_blocks_hashref{$chr}{$block_start}{'marker_raw_index_end'};
	$$events_hashref{$id}{'marker_raw_index_end'} = $$linkage_blocks_hashref{$chr}{$block_end}{'marker_raw_index_start'};
	$$events_hashref{$id}{'marker_effective_index_start'} = $$linkage_blocks_hashref{$chr}{$block_start}{'marker_effective_index_end'};
	$$events_hashref{$id}{'marker_effective_index_end'} = $$linkage_blocks_hashref{$chr}{$block_end}{'marker_effective_index_start'};
    } else {
	# if event is GC
	$$events_hashref{$id}{'marker_pos_start'} = $$linkage_blocks_hashref{$chr}{$block_start}{'marker_pos_start'};
	$$events_hashref{$id}{'marker_pos_end'} = $$linkage_blocks_hashref{$chr}{$block_end}{'marker_pos_end'};
	$$events_hashref{$id}{'marker_raw_index_start'} = $$linkage_blocks_hashref{$chr}{$block_start}{'marker_raw_index_start'};
	$$events_hashref{$id}{'marker_raw_index_end'} = $$linkage_blocks_hashref{$chr}{$block_end}{'marker_raw_index_end'};
	$$events_hashref{$id}{'marker_effective_index_start'} = $$linkage_blocks_hashref{$chr}{$block_start}{'marker_effective_index_start'};
	$$events_hashref{$id}{'marker_effective_index_end'} = $$linkage_blocks_hashref{$chr}{$block_end}{'marker_effective_index_end'};
    }	
    $$events_hashref{$id}{'affected_spores'} = $affected_spores;
}


sub merge_nearby_COs {
    my ($d, $events_hashref, $linkage_blocks_hashref, $CO_associated_GC_hashref, $merging_log_arrayref) = @_;
    my $merging_log;
    my %stack = ();
    $stack{'id'} = "";
    $stack{'type'} = "";
    $stack{'subtype'} = "";
    $stack{'chr'} = "";
    $stack{'block_start'} = "";
    $stack{'block_end'} = "";
    $stack{'marker_pos_start'} = "";
    $stack{'marker_pos_end'} = "";
    $stack{'marker_raw_index_start'} = "";
    $stack{'marker_raw_index_end'} = "";
    $stack{'marker_effective_index_start'} = "";
    $stack{'marker_effective_index_end'} = "";
    $stack{'affected_spores'} = "";

    my @raw_event_ids = sort custom_sort keys %$events_hashref;
    my $raw_event_count = scalar @raw_event_ids;
    if ($raw_event_count > 0) {
	foreach my $id (@raw_event_ids) {
	    # if this event has not been deleted during previous merging
	    if (exists $$events_hashref{$id}) { 
		# print "for debug: id=$id\n";
		my $type = $$events_hashref{$id}{'type'};
		my $subtype = $$events_hashref{$id}{'subtype'};
		my $chr = $$events_hashref{$id}{'chr'};
		my $block_start = $$events_hashref{$id}{'block_start'};
		my $block_end = $$events_hashref{$id}{'block_end'};
		my $marker_pos_start = $$events_hashref{$id}{'marker_pos_start'};
		my $marker_pos_end = $$events_hashref{$id}{'marker_pos_end'};
		my $marker_raw_index_start = $$events_hashref{$id}{'marker_raw_index_start'};
		my $marker_raw_index_end = $$events_hashref{$id}{'marker_raw_index_end'};
		my $marker_effective_index_start = $$events_hashref{$id}{'marker_effective_index_start'};
		my $marker_effective_index_end = $$events_hashref{$id}{'marker_effective_index_end'};
		my $affected_spores = $$events_hashref{$id}{'affected_spores'};
		# print "for debug: id=$id, type=$type, subtype=$subtype, chr=$chr, marker_pos_start=$marker_pos_start, marker_pos_end=$marker_pos_end\n";
		if ($type eq "GC") {
		    # no merging for GCs
		} elsif ($stack{'id'} eq "") {
		    # the current event is single CO
		    # initialize the stack
		    $stack{'id'} = $id;
		    $stack{'type'} = $type;
		    $stack{'subtype'} = $subtype;
		    $stack{'chr'} = $chr;
		    $stack{'block_start'} = $block_start;
		    $stack{'block_end'} = $block_end;
		    $stack{'marker_pos_start'} = $marker_pos_start;
		    $stack{'marker_pos_end'} = $marker_pos_end;
		    $stack{'marker_raw_index_start'} = $marker_raw_index_start;
		    $stack{'marker_raw_index_end'} = $marker_raw_index_end;
		    $stack{'marker_effective_index_start'} = $marker_effective_index_start;
		    $stack{'marker_effective_index_end'} = $marker_effective_index_end;
		    $stack{'affected_spores'} = $affected_spores;
		} elsif ($stack{'chr'} ne $chr) {
		    # no merging but reset stack
		    $stack{'id'} = $id;
		    $stack{'type'} = $type;
		    $stack{'subtype'} = $subtype;
		    $stack{'chr'} = $chr;
		    $stack{'block_start'} = $block_start;
		    $stack{'block_end'} = $block_end;
		    $stack{'marker_pos_start'} = $marker_pos_start;
		    $stack{'marker_pos_end'} = $marker_pos_end;
		    $stack{'marker_raw_index_start'} = $marker_raw_index_start;
		    $stack{'marker_raw_index_end'} = $marker_raw_index_end;
		    $stack{'marker_effective_index_start'} = $marker_effective_index_start;
		    $stack{'marker_effective_index_end'} = $marker_effective_index_end;
		    $stack{'affected_spores'} = $affected_spores;
		} else {
		    # for single COs on the same chromosome, we need to check the exact genomic location.
		    if ($marker_pos_start - $stack{'marker_pos_end'} > $d) {
			# distance > d
			# no merging but reset stack
			$stack{'id'} = $id;
			$stack{'type'} = $type;
			$stack{'subtype'} = $subtype;
			$stack{'chr'} = $chr;
			$stack{'block_start'} = $block_start;
			$stack{'block_end'} = $block_end;
			$stack{'marker_pos_start'} = $marker_pos_start;
			$stack{'marker_pos_end'} = $marker_pos_end;
			$stack{'marker_raw_index_start'} = $marker_raw_index_start;
			$stack{'marker_raw_index_end'} = $marker_raw_index_end;
			$stack{'marker_effective_index_start'} = $marker_effective_index_start;
			$stack{'marker_effective_index_end'} = $marker_effective_index_end;
			$stack{'affected_spores'} = $affected_spores;
		    } else {
			# distance <= d, merging might be needed.
			# check how many chromatids are involved in CO exchange of the two merging COs
			# print "merging expected for $stack{'id'} and $id between $stack{'block_start'} and $block_end\n";
			my @stack_CO_affected_spores = split /:/, $stack{'affected_spores'};
			my @current_CO_affected_spores = split /:/, $affected_spores;
			my $total_CO_affected_spores_count = 0;
			my @total_CO_affected_spores = (0, 0, 0, 0);
			for (my $i = 0; $i < 4; $i++) {
			    if (($stack_CO_affected_spores[$i] == 1) or ($current_CO_affected_spores[$i] == 1)) {
				$total_CO_affected_spores_count++;
				$total_CO_affected_spores[$i] = 1;
			    }
			}
			my $total_CO_affected_spores = join ":", @total_CO_affected_spores;
			my ($genotype_switches_count_pattern, $genotype_switches_count_pattern_sorted, $uniq_single_switch_count, $genotype_switch_boundaries) = check_genotype_switches($stack{'block_start'}, $block_end, \%{$$linkage_blocks_hashref{$chr}});

			my @genotype_switches_count_pattern = split /:/, $genotype_switches_count_pattern;
			my @genotype_switch_boundaries = split /:/, $genotype_switch_boundaries;
			if (($total_CO_affected_spores_count == 2) and ($genotype_switches_count_pattern_sorted eq "0:0:2:2")) {
			    # recategorize as Type6_GC # double NCO			
			    my $type6_GC_id = "$stack{'id'}.1";
			    my $type6_GC_type = "GC";
			    my $type6_GC_subtype = "Type6_GC"; # double NCO or apparent double CO
			    $merging_log = "$stack{'id'} $stack{'subtype'} $stack{'chr'} $stack{'block_start'} $stack{'block_end'} $stack{'marker_pos_start'} $stack{'marker_pos_end'}\n";
			    $merging_log .= "+\n";
			    $merging_log .= "$id $subtype $chr $block_start $block_end $marker_pos_start $marker_pos_end\n";
			    $merging_log .= "=> new $type6_GC_id ";
			    my $type6_GC_chr = $chr;
			    my $type6_GC_block_start = $stack{'block_start'} + 1;
			    my $type6_GC_block_end = $block_end - 1;
			    # for debug: $type6_GC_block_start and $type6_GC_block_end should be the same
			    my $type6_GC_marker_pos_start = $$linkage_blocks_hashref{$type6_GC_chr}{$type6_GC_block_start}{'marker_pos_start'};
			    my $type6_GC_marker_pos_end = $$linkage_blocks_hashref{$type6_GC_chr}{$type6_GC_block_end}{'marker_pos_end'};
			    my $type6_GC_marker_raw_index_start = $$linkage_blocks_hashref{$type6_GC_chr}{$type6_GC_block_start}{'marker_raw_index_start'};
			    my $type6_GC_marker_raw_index_end = $$linkage_blocks_hashref{$type6_GC_chr}{$type6_GC_block_end}{'marker_raw_index_end'};
			    my $type6_GC_marker_effective_index_start = $$linkage_blocks_hashref{$type6_GC_chr}{$type6_GC_block_start}{'marker_effective_index_start'};
			    my $type6_GC_marker_effective_index_end = $$linkage_blocks_hashref{$type6_GC_chr}{$type6_GC_block_end}{'marker_effective_index_end'};
			    my $type6_GC_affected_spores = $affected_spores;
			    $merging_log .= "$type6_GC_subtype $type6_GC_chr $type6_GC_block_start $type6_GC_block_end $type6_GC_marker_pos_start $type6_GC_marker_pos_end\n\n";
			    push @$merging_log_arrayref, $merging_log;
			    # add the newly merged event
			    $$events_hashref{$type6_GC_id}{'type'} = $type6_GC_type;
			    $$events_hashref{$type6_GC_id}{'subtype'} = $type6_GC_subtype;
			    $$events_hashref{$type6_GC_id}{'chr'} = $type6_GC_chr;
			    $$events_hashref{$type6_GC_id}{'block_start'} = $type6_GC_block_start;
			    $$events_hashref{$type6_GC_id}{'block_end'} = $type6_GC_block_end;
			    $$events_hashref{$type6_GC_id}{'marker_pos_start'} = $type6_GC_marker_pos_start;
			    $$events_hashref{$type6_GC_id}{'marker_pos_end'} = $type6_GC_marker_pos_end;
			    $$events_hashref{$type6_GC_id}{'marker_raw_index_start'} = $type6_GC_marker_raw_index_start;
			    $$events_hashref{$type6_GC_id}{'marker_raw_index_end'} = $type6_GC_marker_raw_index_end;
			    $$events_hashref{$type6_GC_id}{'marker_effective_index_start'} = $type6_GC_marker_effective_index_start;
			    $$events_hashref{$type6_GC_id}{'marker_effective_index_end'} = $type6_GC_marker_effective_index_end;
			    $$events_hashref{$type6_GC_id}{'affected_spores'} = $type6_GC_affected_spores;
			    
			    # print "delete CO $stack{'id'}\n";
			    delete $$events_hashref{$stack{'id'}};
			    # print "delete CO $id\n";
			    delete $$events_hashref{$id};
			    foreach my $GC_id (@{$$CO_associated_GC_hashref{$stack{'id'}}}) {
				# print "deleted GC associated with the CO $stack{'id'}: $GC_id\n";
				delete $$events_hashref{$GC_id};
			    }
			    foreach my $GC_id (@{$$CO_associated_GC_hashref{$id}}) {
				# print "deleted GC associated with the CO $id: $GC_id\n";
				delete $$events_hashref{$GC_id};
			    }
			    delete $$CO_associated_GC_hashref{$stack{'id'}};
			    delete $$CO_associated_GC_hashref{$id};
			    # reset stack
			    $stack{'id'} = "";
			    $stack{'type'} = "";
			    $stack{'subtype'} = "";
			    $stack{'chr'} = "";
			    $stack{'block_start'} = "";
			    $stack{'block_end'} = "";
			    $stack{'marker_pos_start'} = "";
			    $stack{'marker_pos_end'} = "";
			    $stack{'marker_raw_index_start'} = "";
			    $stack{'marker_raw_index_end'} = "";
			    $stack{'marker_effective_index_start'} = "";
			    $stack{'marker_effective_index_end'} = "";
			    $stack{'affected_spores'} = "";
			} elsif (($total_CO_affected_spores_count == 2) and ($genotype_switches_count_pattern_sorted ne "0:0:2:2")) {
			    # print "total_CO_affected_spores_count=$total_CO_affected_spores_count, genotype_switches_count_pattern_sorted=$genotype_switches_count_pattern_sorted\n";
			    # print "stack_id=$stack{'id'}, stack_subtype=$stack{'subtype'}, stack_chr=$stack{'chr'}, stack_affected_spores=$stack{'affected_spores'}\n";
			    # print "id=$id, subtype=$subtype, chr=$chr, affected_spores=$affected_spores\n"; 
			    # recategorize as Type17_GC
			    my $type17_GC_id = "$stack{'id'}.1";
			    my $type17_GC_type = "GC";
			    my $type17_GC_subtype = "Type17_GC"; # complex CO and GC combination
			    my $type17_GC_chr = $chr;
			    my $type17_GC_block_start = $stack{'block_start'} + 1;
			    my $type17_GC_block_end = $block_end - 1;
			    my $type17_GC_marker_pos_start = $$linkage_blocks_hashref{$type17_GC_chr}{$type17_GC_block_start}{'marker_pos_start'};
			    my $type17_GC_marker_pos_end = $$linkage_blocks_hashref{$type17_GC_chr}{$type17_GC_block_end}{'marker_pos_end'};
			    my $type17_GC_marker_raw_index_start = $$linkage_blocks_hashref{$type17_GC_chr}{$type17_GC_block_start}{'marker_raw_index_start'};
			    my $type17_GC_marker_raw_index_end = $$linkage_blocks_hashref{$type17_GC_chr}{$type17_GC_block_end}{'marker_raw_index_end'};
			    my $type17_GC_marker_effective_index_start = $$linkage_blocks_hashref{$type17_GC_chr}{$type17_GC_block_start}{'marker_effective_index_start'};
			    my $type17_GC_marker_effective_index_end = $$linkage_blocks_hashref{$type17_GC_chr}{$type17_GC_block_end}{'marker_effective_index_end'};

			    my $type17_GC_affected_spores = $total_CO_affected_spores;

			    $merging_log = "$stack{'id'} $stack{'subtype'} $stack{'chr'} $stack{'block_start'} $stack{'block_end'} $stack{'marker_pos_start'} $stack{'marker_pos_end'}\n";
			    $merging_log .= "+\n";
			    $merging_log .= "$id $subtype $chr $block_start $block_end $marker_pos_start $marker_pos_end\n";
			    $merging_log .= "=> new $type17_GC_id ";
			    $merging_log .= "$type17_GC_subtype $type17_GC_chr $type17_GC_block_start $type17_GC_block_end $type17_GC_marker_pos_start $type17_GC_marker_pos_end\n\n";
			    
			    push @$merging_log_arrayref, $merging_log;
			    
			    # add the newly merged event
			    $$events_hashref{$type17_GC_id}{'type'} = $type17_GC_type;
			    $$events_hashref{$type17_GC_id}{'subtype'} = $type17_GC_subtype;
			    $$events_hashref{$type17_GC_id}{'chr'} = $type17_GC_chr;
			    $$events_hashref{$type17_GC_id}{'block_start'} = $type17_GC_block_start;
			    $$events_hashref{$type17_GC_id}{'block_end'} = $type17_GC_block_end;
			    $$events_hashref{$type17_GC_id}{'marker_pos_start'} = $type17_GC_marker_pos_start;
			    $$events_hashref{$type17_GC_id}{'marker_pos_end'} = $type17_GC_marker_pos_end;
			    $$events_hashref{$type17_GC_id}{'marker_raw_index_start'} = $type17_GC_marker_raw_index_start;
			    $$events_hashref{$type17_GC_id}{'marker_raw_index_end'} = $type17_GC_marker_raw_index_end;
			    $$events_hashref{$type17_GC_id}{'marker_effective_index_start'} = $type17_GC_marker_effective_index_start;
			    $$events_hashref{$type17_GC_id}{'marker_effective_index_end'} = $type17_GC_marker_effective_index_end;
			    $$events_hashref{$type17_GC_id}{'affected_spores'} = $type17_GC_affected_spores;
			    
			    # print "delete CO $stack{'id'}\n";
			    delete $$events_hashref{$stack{'id'}};
			    # print "delete CO $id\n";
			    delete $$events_hashref{$id};
			    foreach my $GC_id (@{$$CO_associated_GC_hashref{$stack{'id'}}}) {
				# print "deleted GC associated with the CO $stack{'id'}: $GC_id\n";
				delete $$events_hashref{$GC_id};
			    }
			    foreach my $GC_id (@{$$CO_associated_GC_hashref{$id}}) {
				# print "deleted GC associated with the CO $id: $GC_id\n";
				delete $$events_hashref{$GC_id};
			    }
			    delete $$CO_associated_GC_hashref{$stack{'id'}};
			    delete $$CO_associated_GC_hashref{$id};
			    
			    # reset stack
			    $stack{'id'} = "";
			    $stack{'type'} = "";
			    $stack{'subtype'} = "";
			    $stack{'chr'} = "";
			    $stack{'block_start'} = "";
			    $stack{'block_end'} = "";
			    $stack{'marker_pos_start'} = "";
			    $stack{'marker_pos_end'} = "";
			    $stack{'marker_raw_index_start'} = "";
			    $stack{'marker_raw_index_end'} = "";
			    $stack{'marker_effective_index_start'} = "";
			    $stack{'marker_effective_index_end'} = "";
			    $stack{'affected_spores'} = "";
			} elsif (($total_CO_affected_spores_count == 3) and ($genotype_switches_count_pattern_sorted eq "0:1:1:2")) {
			    # recategorize as Type5_CO + Type10_GC + Type11_GC
			    my $type5_CO_id = "$stack{'id'}.1";
			    my $type5_CO_type = "CO";
			    my $type5_CO_subtype = "Type5_CO";
			    my $type10_GC_id = "$stack{'id'}.2";
			    my $type10_GC_type = "GC";
			    my $type10_GC_subtype = "Type10_GC";
			    my $type11_GC_id = "$stack{'id'}.3";
			    my $type11_GC_type = "GC";
			    my $type11_GC_subtype = "Type11_GC";
			    @{$$CO_associated_GC_hashref{$type5_CO_id}} = ();
			    
			    my $type5_CO_chr = $chr;
			    my $type5_CO_block_start;
			    my $type5_CO_block_end;
			    my $type5_CO_marker_pos_start;
			    my $type5_CO_marker_pos_end;
			    my $type5_CO_marker_raw_index_start;
			    my $type5_CO_marker_raw_index_end;
			    my $type5_CO_marker_effective_index_start;
			    my $type5_CO_marker_effective_index_end;
			    
			    my $type10_GC_chr = $chr;
			    my $type10_GC_block_start;
			    my $type10_GC_block_end;
			    my $type10_GC_marker_pos_start;
			    my $type10_GC_marker_pos_end;
			    my $type10_GC_marker_raw_index_start;
			    my $type10_GC_marker_raw_index_end;
			    my $type10_GC_marker_effective_index_start;
			    my $type10_GC_marker_effective_index_end;
			    
			    my $type11_GC_chr = $chr;
			    my $type11_GC_block_start;
			    my $type11_GC_block_end;
			    my $type11_GC_marker_pos_start;
			    my $type11_GC_marker_pos_end;
			    my $type11_GC_marker_raw_index_start;
			    my $type11_GC_marker_raw_index_end;
			    my $type11_GC_marker_effective_index_start;
			    my $type11_GC_marker_effective_index_end;
			    
			    my @type5_CO_affected_spores = (0, 0, 0, 0);
			    my @type10_GC_affected_spores = (0, 0, 0, 0);
			    my @type11_GC_affected_spores = (0, 0, 0, 0);
			    
			    my @type5_CO_boundary_blocks = ();
			    my @type10_GC_boundary_blocks = ();
			    
			    for (my $i = 0; $i < 4; $i++) {
				# print "i=$i, genotype_switch_count_pattern=$genotype_switches_count_pattern[$i]\n";
				if ($genotype_switches_count_pattern[$i] == 2) {
				    # Type11_GC
				    $type11_GC_affected_spores[$i] = 1;
				    # print "i=$i, genotype_switch_boundaries=$genotype_switch_boundaries[$i]\n";
				    my ($b1, $b2) = split /,/, $genotype_switch_boundaries[$i];
				    # print "b1=$b1, b2=$b2\n";
				    ($type11_GC_block_start) = ($b1 =~ /\|(\S+)/);
				    ($type11_GC_block_end) = ($b2 =~ /(\S+)\|/);
				    $type11_GC_marker_pos_start = $$linkage_blocks_hashref{$type11_GC_chr}{$type11_GC_block_start}{'marker_pos_start'};
				    $type11_GC_marker_pos_end = $$linkage_blocks_hashref{$type11_GC_chr}{$type11_GC_block_end}{'marker_pos_end'};
				    $type11_GC_marker_raw_index_start = $$linkage_blocks_hashref{$type11_GC_chr}{$type11_GC_block_start}{'marker_raw_index_start'};
				    $type11_GC_marker_raw_index_end = $$linkage_blocks_hashref{$type11_GC_chr}{$type11_GC_block_end}{'marker_raw_index_end'};
				    $type11_GC_marker_effective_index_start = $$linkage_blocks_hashref{$type11_GC_chr}{$type11_GC_block_start}{'marker_effective_index_start'};
				    $type11_GC_marker_effective_index_end = $$linkage_blocks_hashref{$type11_GC_chr}{$type11_GC_block_end}{'marker_effective_index_end'};
				    # add the newly merged event
				    $$events_hashref{$type11_GC_id}{'type'} = $type11_GC_type;
				    $$events_hashref{$type11_GC_id}{'subtype'} = $type11_GC_subtype;
				    $$events_hashref{$type11_GC_id}{'chr'} = $type11_GC_chr;
				    $$events_hashref{$type11_GC_id}{'block_start'} = $type11_GC_block_start;
				    $$events_hashref{$type11_GC_id}{'block_end'} = $type11_GC_block_end;
				    $$events_hashref{$type11_GC_id}{'marker_pos_start'} = $type11_GC_marker_pos_start;
				    $$events_hashref{$type11_GC_id}{'marker_pos_end'} = $type11_GC_marker_pos_end;
				    $$events_hashref{$type11_GC_id}{'marker_raw_index_start'} = $type11_GC_marker_raw_index_start;
				    $$events_hashref{$type11_GC_id}{'marker_raw_index_end'} = $type11_GC_marker_raw_index_end;
				    $$events_hashref{$type11_GC_id}{'marker_effective_index_start'} = $type11_GC_marker_effective_index_start;
				    $$events_hashref{$type11_GC_id}{'marker_effective_index_end'} = $type11_GC_marker_effective_index_end;
				    $$events_hashref{$type11_GC_id}{'affected_spores'} = join ":", @type11_GC_affected_spores;
				    # print "type11_GC_id: $type11_GC_id, affected_spores: $$events_hashref{$type11_GC_id}{'affected_spores'}\n";
				    push @{$$CO_associated_GC_hashref{$type5_CO_id}}, $type11_GC_id;
				} elsif ($genotype_switches_count_pattern[$i] == 1) {
				    # Type5_CO and Type10_GC
				    $type5_CO_affected_spores[$i] = 1;
				    my ($b1, $b2) = split /\|/, $genotype_switch_boundaries[$i];
				    push @type5_CO_boundary_blocks, $b1;
				    push @type5_CO_boundary_blocks, $b2;
				    # test if Type10_GC exists or not
				    if ($uniq_single_switch_count > 1) {
					$type10_GC_affected_spores[$i] = 1;
					push @type10_GC_boundary_blocks, $b1;
					push @type10_GC_boundary_blocks, $b2;
				    } else {
					die "Type10_GC should exist. Exit when this is not the case."
				    }
				}
			    }
			    # add the outter boundries of type11_GC for comparison:
			    # push @type5_CO_boundary_blocks, $type11_GC_block_start - 1;
			    # push @type5_CO_boundary_blocks, $type11_GC_block_end + 1;
			    
			    my @type5_CO_boundary_blocks_sorted = sort {$a <=> $b} @type5_CO_boundary_blocks;
			    # print "type5_CO_boundary_blocks_sorted: @type5_CO_boundary_blocks_sorted\n";
			    # extract the two outtermost boundaries 
			    $type5_CO_block_start = $type5_CO_boundary_blocks_sorted[0];
			    $type5_CO_block_end = $type5_CO_boundary_blocks_sorted[-1];
			    # print "type5_CO_block_start=$type5_CO_block_start, type5_CO_block_end=$type5_CO_block_end\n";
			    $type5_CO_marker_pos_start = $$linkage_blocks_hashref{$type5_CO_chr}{$type5_CO_block_start}{'marker_pos_end'};
			    $type5_CO_marker_pos_end = $$linkage_blocks_hashref{$type5_CO_chr}{$type5_CO_block_end}{'marker_pos_start'};
			    # print "type5_CO_marker_pos_start=$type5_CO_marker_pos_start, type5_CO_marker_pos_end=$type5_CO_marker_pos_end\n";
			    $type5_CO_marker_raw_index_start = $$linkage_blocks_hashref{$type5_CO_chr}{$type5_CO_block_start}{'marker_raw_index_end'};
			    $type5_CO_marker_raw_index_end = $$linkage_blocks_hashref{$type5_CO_chr}{$type5_CO_block_end}{'marker_raw_index_start'};
			    $type5_CO_marker_effective_index_start = $$linkage_blocks_hashref{$type5_CO_chr}{$type5_CO_block_start}{'marker_effective_index_end'};
			    $type5_CO_marker_effective_index_end = $$linkage_blocks_hashref{$type5_CO_chr}{$type5_CO_block_end}{'marker_effective_index_start'};
			    # add the newly merged event
			    $$events_hashref{$type5_CO_id}{'type'} = $type5_CO_type;
			    $$events_hashref{$type5_CO_id}{'subtype'} = $type5_CO_subtype;
			    $$events_hashref{$type5_CO_id}{'chr'} = $type5_CO_chr;
			    $$events_hashref{$type5_CO_id}{'block_start'} = $type5_CO_block_start;
			    $$events_hashref{$type5_CO_id}{'block_end'} = $type5_CO_block_end;
			    $$events_hashref{$type5_CO_id}{'marker_pos_start'} = $type5_CO_marker_pos_start;
			    $$events_hashref{$type5_CO_id}{'marker_pos_end'} = $type5_CO_marker_pos_end;
			    $$events_hashref{$type5_CO_id}{'marker_raw_index_start'} = $type5_CO_marker_raw_index_start;
			    $$events_hashref{$type5_CO_id}{'marker_raw_index_end'} = $type5_CO_marker_raw_index_end;
			    $$events_hashref{$type5_CO_id}{'marker_effective_index_start'} = $type5_CO_marker_effective_index_start;
			    $$events_hashref{$type5_CO_id}{'marker_effective_index_end'} = $type5_CO_marker_effective_index_end;
			    $$events_hashref{$type5_CO_id}{'affected_spores'} = join ":", @type5_CO_affected_spores;
			    
			    if ($uniq_single_switch_count > 1) {
				# type10 GC does exist
				my @type10_GC_boundary_blocks_sorted = sort {$a <=> $b} @type10_GC_boundary_blocks;
				$type10_GC_block_start = $type10_GC_boundary_blocks_sorted[1];
				$type10_GC_block_end = $type10_GC_boundary_blocks_sorted[2];
				$type10_GC_marker_pos_start = $$linkage_blocks_hashref{$type10_GC_chr}{$type10_GC_block_start}{'marker_pos_start'}; 
				$type10_GC_marker_pos_end = $$linkage_blocks_hashref{$type10_GC_chr}{$type10_GC_block_end}{'marker_pos_end'}; 
				$type10_GC_marker_raw_index_start = $$linkage_blocks_hashref{$type10_GC_chr}{$type10_GC_block_start}{'marker_raw_index_start'}; 
				$type10_GC_marker_raw_index_end = $$linkage_blocks_hashref{$type10_GC_chr}{$type10_GC_block_end}{'marker_raw_index_end'}; 
				$type10_GC_marker_effective_index_start = $$linkage_blocks_hashref{$type10_GC_chr}{$type10_GC_block_start}{'marker_effective_index_start'}; 
				$type10_GC_marker_effective_index_end = $$linkage_blocks_hashref{$type10_GC_chr}{$type10_GC_block_end}{'marker_effective_index_end'}; 
				# add the newly merged event
				$$events_hashref{$type10_GC_id}{'type'} = $type10_GC_type;
				$$events_hashref{$type10_GC_id}{'subtype'} = $type10_GC_subtype;
				$$events_hashref{$type10_GC_id}{'chr'} = $type10_GC_chr;
				$$events_hashref{$type10_GC_id}{'block_start'} = $type10_GC_block_start;
				$$events_hashref{$type10_GC_id}{'block_end'} = $type10_GC_block_end;
				$$events_hashref{$type10_GC_id}{'marker_pos_start'} = $type10_GC_marker_pos_start;
				$$events_hashref{$type10_GC_id}{'marker_pos_end'} =  $type10_GC_marker_pos_end;
				$$events_hashref{$type10_GC_id}{'marker_raw_index_start'} = $type10_GC_marker_raw_index_start;
				$$events_hashref{$type10_GC_id}{'marker_raw_index_end'} =  $type10_GC_marker_raw_index_end;
				$$events_hashref{$type10_GC_id}{'marker_effective_index_start'} = $type10_GC_marker_effective_index_start;
				$$events_hashref{$type10_GC_id}{'marker_effective_index_end'} =  $type10_GC_marker_effective_index_end;
				$$events_hashref{$type10_GC_id}{'affected_spores'} = join ":", @type10_GC_affected_spores;
				push @{$$CO_associated_GC_hashref{$type5_CO_id}}, $type10_GC_id;
			    }
			    $merging_log = "$stack{'id'} $stack{'subtype'} $stack{'chr'} $stack{'block_start'} $stack{'block_end'} $stack{'marker_pos_start'} $stack{'marker_pos_end'}\n";
			    $merging_log .= "+\n";
			    $merging_log .= "$id $subtype $chr $block_start $block_end $marker_pos_start $marker_pos_end\n";
			    $merging_log .= "=> new $type5_CO_id ";
			    $merging_log .= "$type5_CO_subtype $type5_CO_chr $type5_CO_block_start $type5_CO_block_end $type5_CO_marker_pos_start $type5_CO_marker_pos_end\n";
			    if ($uniq_single_switch_count != 1) {
				$merging_log .= "+ new $type10_GC_id ";
				$merging_log .= "$type10_GC_subtype $type10_GC_chr $type10_GC_block_start $type10_GC_block_end $type10_GC_marker_pos_start $type10_GC_marker_pos_end\n";
			    }
			    $merging_log .= "+ new $type11_GC_id ";
			    $merging_log .= "$type11_GC_subtype $type11_GC_chr $type11_GC_block_start $type11_GC_block_end $type11_GC_marker_pos_start $type11_GC_marker_pos_end\n";
			    push @$merging_log_arrayref, $merging_log;
			    
			    # print "delete CO $stack{'id'}\n";
			    delete $$events_hashref{$stack{'id'}};
			    # print "delete CO $id\n";
			    delete $$events_hashref{$id};
			    foreach my $GC_id (@{$$CO_associated_GC_hashref{$stack{'id'}}}) {
				# print "deleted GC associated with the CO $stack{'id'}: $GC_id\n";
				delete $$events_hashref{$GC_id};
			    }
			    foreach my $GC_id (@{$$CO_associated_GC_hashref{$id}}) {
				# print "deleted GC associated with the CO $id: $GC_id\n";
				delete $$events_hashref{$GC_id};
			    }
			    delete $$CO_associated_GC_hashref{$stack{'id'}};
			    delete $$CO_associated_GC_hashref{$id};
			    
			    # reset stack
			    $stack{'id'} = $type5_CO_id;
			    $stack{'type'} = $type5_CO_type;
			    $stack{'subtype'} = $type5_CO_subtype;
			    $stack{'chr'} = $type5_CO_chr;
			    $stack{'block_start'} = $type5_CO_block_start;
			    $stack{'block_end'} = $type5_CO_block_end;
			    $stack{'marker_pos_start'} = $type5_CO_marker_pos_start;
			    $stack{'marker_pos_end'} = $type5_CO_marker_pos_end;
			    $stack{'marker_raw_index_start'} = $type5_CO_marker_raw_index_start;
			    $stack{'marker_raw_index_end'} = $type5_CO_marker_raw_index_end;
			    $stack{'marker_effective_index_start'} = $type5_CO_marker_effective_index_start;
			    $stack{'marker_effective_index_end'} = $type5_CO_marker_effective_index_end;
			    $stack{'affected_spores'} = join ":", @type5_CO_affected_spores;
			} elsif (($total_CO_affected_spores_count == 3) and ($genotype_switches_count_pattern_sorted ne "0:1:1:2")) {
			    # print "total_CO_affected_spores_count=2, genotype_switches_count_pattern_sorted=$genotype_switches_count_pattern_sorted\n";
			    # print "stack_id=$stack{'id'}, stack_subtype=$stack{'subtype'}, stack_chr=$stack{'chr'}\n";
			    # print "stack_block_start=$stack{'block_start'}, stack_block_end=$stack{'block_end'}, stack_affected_spores=$stack{'affected_spores'}\n";
			    # print "id=$id, subtype=$subtype, chr=$chr\n";
			    # print "block_start=$block_start, block_end=$block_end, affected_spores=$affected_spores\n";
			    # print "recategorize as Type6_CO\n";
			    # recategorize as Type6_CO
			    my $type6_CO_id = "$stack{'id'}.1";
			    my $type6_CO_type = "CO";
			    my $type6_CO_subtype = "Type6_CO"; # single CO with complex GC pattern
			    my $type6_CO_chr = $chr;
			    my $type6_CO_block_start = $stack{'block_start'};
			    my $type6_CO_block_end = $block_end;
			    my $type6_CO_marker_pos_start = $$linkage_blocks_hashref{$type6_CO_chr}{$type6_CO_block_start}{'marker_pos_end'};
			    my $type6_CO_marker_pos_end = $$linkage_blocks_hashref{$type6_CO_chr}{$type6_CO_block_end}{'marker_pos_start'};
			    my $type6_CO_marker_raw_index_start = $$linkage_blocks_hashref{$type6_CO_chr}{$type6_CO_block_start}{'marker_raw_index_end'};
			    my $type6_CO_marker_raw_index_end = $$linkage_blocks_hashref{$type6_CO_chr}{$type6_CO_block_end}{'marker_raw_index_start'};
			    my $type6_CO_marker_effective_index_start = $$linkage_blocks_hashref{$type6_CO_chr}{$type6_CO_block_start}{'marker_effective_index_end'};
			    my $type6_CO_marker_effective_index_end = $$linkage_blocks_hashref{$type6_CO_chr}{$type6_CO_block_end}{'marker_effective_index_start'};
			    my @type6_CO_affected_spores = (0, 0, 0, 0);
			    @{$$CO_associated_GC_hashref{$type6_CO_id}} = ();
			    my @genotype_switches_count_pattern = split /:/, $genotype_switches_count_pattern;
			    my $type15_GC_affected_spore_index;
			    my $type6_CO_affected_spores_count = 0;
			    for (my $i = 0; $i < 4; $i++) {
				if ($genotype_switches_count_pattern[$i] == 1) {
				    $type6_CO_affected_spores[$i] = 1;
				    $type6_CO_affected_spores_count++;
				} elsif ($genotype_switches_count_pattern[$i] > 1) {
				    $type15_GC_affected_spore_index = $i;
				}
			    }

			    if ($type6_CO_affected_spores_count < 2) {
				$type6_CO_affected_spores[$type15_GC_affected_spore_index] = 1;
				$type6_CO_affected_spores_count++;
			    }
			    my $type6_CO_affected_spores = join ":", @type6_CO_affected_spores;
			    
			    my $type15_GC_id = "$stack{'id'}.2";
			    my $type15_GC_type = "GC";
			    my $type15_GC_subtype = "Type15_GC"; # complex GC associated with Type6_CO
			    push @{$$CO_associated_GC_hashref{$type6_CO_id}}, $type15_GC_id;
			    my $type15_GC_chr = $type6_CO_chr;

			    
			    my @type15_GC_block_boundaries = ();
			    my $type15_GC_block_start = $type6_CO_block_start + 1;
			    my $type15_GC_block_end = $type6_CO_block_end - 1;
			    my @type15_GC_affected_spores = @type6_CO_affected_spores;

			    for (my $i = 0; $i < 4; $i++) {
				if ($genotype_switches_count_pattern[$i] > 1) {
				    $type15_GC_affected_spores[$i] = "1*";
				}
			    }
			    my $type15_GC_affected_spores = join ":", @type15_GC_affected_spores;
			    my $type15_GC_marker_pos_start = $$linkage_blocks_hashref{$type15_GC_chr}{$type15_GC_block_start}{'marker_pos_start'};
			    my $type15_GC_marker_pos_end = $$linkage_blocks_hashref{$type15_GC_chr}{$type15_GC_block_end}{'marker_pos_end'};
			    my $type15_GC_marker_raw_index_start = $$linkage_blocks_hashref{$type15_GC_chr}{$type15_GC_block_start}{'marker_raw_index_start'};
			    my $type15_GC_marker_raw_index_end = $$linkage_blocks_hashref{$type15_GC_chr}{$type15_GC_block_end}{'marker_raw_index_end'};
			    my $type15_GC_marker_effective_index_start = $$linkage_blocks_hashref{$type15_GC_chr}{$type15_GC_block_start}{'marker_effective_index_start'};
			    my $type15_GC_marker_effective_index_end = $$linkage_blocks_hashref{$type15_GC_chr}{$type15_GC_block_end}{'marker_effective_index_end'};
			    $merging_log = "$stack{'id'} $stack{'subtype'} $stack{'chr'} $stack{'block_start'} $stack{'block_end'} $stack{'marker_pos_start'} $stack{'marker_pos_end'}\n";
			    $merging_log .= "+\n";
			    $merging_log .= "$id $subtype $chr $block_start $block_end $marker_pos_start $marker_pos_end\n";
			    $merging_log .= "=> new $type6_CO_id ";
			    
			    $merging_log .= "$type6_CO_subtype $type6_CO_chr $type6_CO_block_start $type6_CO_block_end $type6_CO_marker_pos_start $type6_CO_marker_pos_end\n\n";
			    $merging_log .= "\+ new $type15_GC_id ";
			    $merging_log .= "$type15_GC_subtype $type15_GC_chr $type15_GC_block_start $type15_GC_block_end $type15_GC_marker_pos_start $type15_GC_marker_pos_end\n\n";
			    push @$merging_log_arrayref, $merging_log;
			    
			    # add the newly merged event
			    $$events_hashref{$type6_CO_id}{'type'} = $type6_CO_type;
			    $$events_hashref{$type6_CO_id}{'subtype'} = $type6_CO_subtype;
			    $$events_hashref{$type6_CO_id}{'chr'} = $type6_CO_chr;
			    $$events_hashref{$type6_CO_id}{'block_start'} = $type6_CO_block_start;
			    $$events_hashref{$type6_CO_id}{'block_end'} = $type6_CO_block_end;
			    $$events_hashref{$type6_CO_id}{'marker_pos_start'} = $type6_CO_marker_pos_start;
			    $$events_hashref{$type6_CO_id}{'marker_pos_end'} = $type6_CO_marker_pos_end;
			    $$events_hashref{$type6_CO_id}{'marker_raw_index_start'} = $type6_CO_marker_raw_index_start;
			    $$events_hashref{$type6_CO_id}{'marker_raw_index_end'} = $type6_CO_marker_raw_index_end;
			    $$events_hashref{$type6_CO_id}{'marker_effective_index_start'} = $type6_CO_marker_effective_index_start;
			    $$events_hashref{$type6_CO_id}{'marker_effective_index_end'} = $type6_CO_marker_effective_index_end;
			    $$events_hashref{$type6_CO_id}{'affected_spores'} = $type6_CO_affected_spores;
			    
			    $$events_hashref{$type15_GC_id}{'type'} = $type15_GC_type;
			    $$events_hashref{$type15_GC_id}{'subtype'} = $type15_GC_subtype;
			    $$events_hashref{$type15_GC_id}{'chr'} = $type15_GC_chr;
			    $$events_hashref{$type15_GC_id}{'block_start'} = $type15_GC_block_start;
			    $$events_hashref{$type15_GC_id}{'block_end'} = $type15_GC_block_end;
			    $$events_hashref{$type15_GC_id}{'marker_pos_start'} = $type15_GC_marker_pos_start;
			    $$events_hashref{$type15_GC_id}{'marker_pos_end'} = $type15_GC_marker_pos_end;
			    $$events_hashref{$type15_GC_id}{'marker_raw_index_start'} = $type15_GC_marker_raw_index_start;
			    $$events_hashref{$type15_GC_id}{'marker_raw_index_end'} = $type15_GC_marker_raw_index_end;
			    $$events_hashref{$type15_GC_id}{'marker_effective_index_start'} = $type15_GC_marker_effective_index_start;
			    $$events_hashref{$type15_GC_id}{'marker_effective_index_end'} = $type15_GC_marker_effective_index_end;
			    $$events_hashref{$type15_GC_id}{'affected_spores'} = $type15_GC_affected_spores;
			    
			    # print "delete CO $stack{'id'}\n";
			    delete $$events_hashref{$stack{'id'}};
			    # print "delete CO $id\n";
			    delete $$events_hashref{$id};
			    foreach my $GC_id (@{$$CO_associated_GC_hashref{$stack{'id'}}}) {
				# print "deleted GC associated with the CO $stack{'id'}: $GC_id\n";
				delete $$events_hashref{$GC_id};
			    }
			    foreach my $GC_id (@{$$CO_associated_GC_hashref{$id}}) {
				# print "deleted GC associated with the CO $id: $GC_id\n";
				delete $$events_hashref{$GC_id};
			    }
			    delete $$CO_associated_GC_hashref{$stack{'id'}};
			    delete $$CO_associated_GC_hashref{$id};
			    
			    # reset stack
			    $stack{'id'} = $type6_CO_id;
			    $stack{'type'} = $type6_CO_type;
			    $stack{'subtype'} = $type6_CO_subtype;
			    $stack{'chr'} = $type6_CO_chr;
			    $stack{'block_start'} = $type6_CO_block_start;
			    $stack{'block_end'} = $type6_CO_block_end;
			    $stack{'marker_pos_start'} = $type6_CO_marker_pos_start;
			    $stack{'marker_pos_end'} = $type6_CO_marker_pos_end;
			    $stack{'marker_raw_index_start'} = $type6_CO_marker_raw_index_start;
			    $stack{'marker_raw_index_end'} = $type6_CO_marker_raw_index_end;
			    $stack{'marker_effective_index_start'} = $type6_CO_marker_effective_index_start;
			    $stack{'marker_effective_index_end'} = $type6_CO_marker_effective_index_end;
			    $stack{'affected_spores'} = $type6_CO_affected_spores;
			    
			} elsif (($total_CO_affected_spores_count == 4) and ($genotype_switches_count_pattern_sorted eq "1:1:1:1")) {
			    # recategorize as Type8_CO + Type12_GC
			    my $type8_CO_id = "$stack{'id'}.1";
			    my $type8_CO_type = "CO";
			    my $type8_CO_subtype = "Type8_CO"; # double CO with one or two associated GC
			    my $type8_CO_chr = $chr;
			    my $type8_CO_block_start = $stack{'block_start'};
			    my $type8_CO_block_end = $block_end;
			    my $type8_CO_marker_pos_start = $$linkage_blocks_hashref{$type8_CO_chr}{$type8_CO_block_start}{'marker_pos_end'};
			    my $type8_CO_marker_pos_end = $$linkage_blocks_hashref{$type8_CO_chr}{$type8_CO_block_end}{'marker_pos_start'};
			    my $type8_CO_marker_raw_index_start = $$linkage_blocks_hashref{$type8_CO_chr}{$type8_CO_block_start}{'marker_raw_index_end'};
			    my $type8_CO_marker_raw_index_end = $$linkage_blocks_hashref{$type8_CO_chr}{$type8_CO_block_end}{'marker_raw_index_start'};
			    my $type8_CO_marker_effective_index_start = $$linkage_blocks_hashref{$type8_CO_chr}{$type8_CO_block_start}{'marker_effective_index_end'};
			    my $type8_CO_marker_effective_index_end = $$linkage_blocks_hashref{$type8_CO_chr}{$type8_CO_block_end}{'marker_effective_index_start'};
			    my $type8_CO_affected_spores = "1:1:1:1";
			    @{$$CO_associated_GC_hashref{$type8_CO_id}} = ();
			    
			    my $type12_GC_id = "$stack{'id'}.2";
			    my $type12_GC_type = "GC";
			    my $type12_GC_subtype = "Type12_GC"; # complex GC associated with Type8_CO
			    push @{$$CO_associated_GC_hashref{$type8_CO_id}}, $type12_GC_id;
			    my $type12_GC_block_start = $type8_CO_block_start + 1;
			    my $type12_GC_block_end = $type8_CO_block_end - 1;
			    my $type12_GC_chr = $type8_CO_chr;
			    my $type12_GC_marker_pos_start = $$linkage_blocks_hashref{$type12_GC_chr}{$type12_GC_block_start}{'marker_pos_start'};
			    my $type12_GC_marker_pos_end = $$linkage_blocks_hashref{$type12_GC_chr}{$type12_GC_block_end}{'marker_pos_end'};
			    my $type12_GC_marker_raw_index_start = $$linkage_blocks_hashref{$type12_GC_chr}{$type12_GC_block_start}{'marker_raw_index_start'};
			    my $type12_GC_marker_raw_index_end = $$linkage_blocks_hashref{$type12_GC_chr}{$type12_GC_block_end}{'marker_raw_index_end'};
			    my $type12_GC_marker_effective_index_start = $$linkage_blocks_hashref{$type12_GC_chr}{$type12_GC_block_start}{'marker_effective_index_start'};
			    my $type12_GC_marker_effective_index_end = $$linkage_blocks_hashref{$type12_GC_chr}{$type12_GC_block_end}{'marker_effective_index_end'};
			    my $type12_GC_affected_spores = "1:1:1:1";

			    $merging_log = "$stack{'id'} $stack{'subtype'} $stack{'chr'} $stack{'block_start'} $stack{'block_end'} $stack{'marker_pos_start'} $stack{'marker_pos_end'}\n";
			    $merging_log .= "+\n";
			    $merging_log .= "$id $subtype $chr $block_start $block_end $marker_pos_start $marker_pos_end\n";
			    $merging_log .= "=> new $type8_CO_id ";
			    
			    $merging_log .= "$type8_CO_subtype $type8_CO_chr $type8_CO_block_start $type8_CO_block_end $type8_CO_marker_pos_start $type8_CO_marker_pos_end\n\n";
			    
			    $merging_log .= "\+ new $type12_GC_id ";
			    $merging_log .= "$type12_GC_subtype $type12_GC_chr $type12_GC_block_start $type12_GC_block_end $type12_GC_marker_pos_start $type12_GC_marker_pos_end\n\n";
			    
			    push @$merging_log_arrayref, $merging_log;
			    
			    # add the newly merged event
			    $$events_hashref{$type8_CO_id}{'type'} = $type8_CO_type;
			    $$events_hashref{$type8_CO_id}{'subtype'} = $type8_CO_subtype;
			    $$events_hashref{$type8_CO_id}{'chr'} = $type8_CO_chr;
			    $$events_hashref{$type8_CO_id}{'block_start'} = $type8_CO_block_start;
			    $$events_hashref{$type8_CO_id}{'block_end'} = $type8_CO_block_end;
			    $$events_hashref{$type8_CO_id}{'marker_pos_start'} = $type8_CO_marker_pos_start;
			    $$events_hashref{$type8_CO_id}{'marker_pos_end'} = $type8_CO_marker_pos_end;
			    $$events_hashref{$type8_CO_id}{'marker_raw_index_start'} = $type8_CO_marker_raw_index_start;
			    $$events_hashref{$type8_CO_id}{'marker_raw_index_end'} = $type8_CO_marker_raw_index_end;
			    $$events_hashref{$type8_CO_id}{'marker_effective_index_start'} = $type8_CO_marker_effective_index_start;
			    $$events_hashref{$type8_CO_id}{'marker_effective_index_end'} = $type8_CO_marker_effective_index_end;
			    $$events_hashref{$type8_CO_id}{'affected_spores'} = $type8_CO_affected_spores;
			    
			    $$events_hashref{$type12_GC_id}{'type'} = $type12_GC_type;
			    $$events_hashref{$type12_GC_id}{'subtype'} = $type12_GC_subtype;
			    $$events_hashref{$type12_GC_id}{'chr'} = $type12_GC_chr;
			    $$events_hashref{$type12_GC_id}{'block_start'} = $type12_GC_block_start;
			    $$events_hashref{$type12_GC_id}{'block_end'} = $type12_GC_block_end;
			    $$events_hashref{$type12_GC_id}{'marker_pos_start'} = $type12_GC_marker_pos_start;
			    $$events_hashref{$type12_GC_id}{'marker_pos_end'} = $type12_GC_marker_pos_end;
			    $$events_hashref{$type12_GC_id}{'marker_raw_index_start'} = $type12_GC_marker_raw_index_start;
			    $$events_hashref{$type12_GC_id}{'marker_raw_index_end'} = $type12_GC_marker_raw_index_end;
			    $$events_hashref{$type12_GC_id}{'marker_effective_index_start'} = $type12_GC_marker_effective_index_start;
			    $$events_hashref{$type12_GC_id}{'marker_effective_index_end'} = $type12_GC_marker_effective_index_end;
			    $$events_hashref{$type12_GC_id}{'affected_spores'} = $type12_GC_affected_spores;
			    
			    # print "delete CO $stack{'id'}\n";
			    delete $$events_hashref{$stack{'id'}};
			    # print "delete CO $id\n";
			    delete $$events_hashref{$id};
			    foreach my $GC_id (@{$$CO_associated_GC_hashref{$stack{'id'}}}) {
				# print "deleted GC associated with the CO $stack{'id'}: $GC_id\n";
				delete $$events_hashref{$GC_id};
			    }
			    foreach my $GC_id (@{$$CO_associated_GC_hashref{$id}}) {
				# print "deleted GC associated with the CO $id: $GC_id\n";
				delete $$events_hashref{$GC_id};
			    }
			    delete $$CO_associated_GC_hashref{$stack{'id'}};
			    delete $$CO_associated_GC_hashref{$id};
			    
			    $stack{'id'} = $type8_CO_id;
			    $stack{'type'} = $type8_CO_type;
			    $stack{'subtype'} = $type8_CO_subtype;
			    $stack{'chr'} = $type8_CO_chr;
			    $stack{'block_start'} = $type8_CO_block_start;
			    $stack{'block_end'} = $type8_CO_block_end;
			    $stack{'marker_pos_start'} = $type8_CO_marker_pos_start;
			    $stack{'marker_pos_end'} = $type8_CO_marker_pos_end;
			    $stack{'marker_raw_index_start'} = $type8_CO_marker_raw_index_start;
			    $stack{'marker_raw_index_end'} = $type8_CO_marker_raw_index_end;
			    $stack{'marker_effective_index_start'} = $type8_CO_marker_effective_index_start;
			    $stack{'marker_effective_index_end'} = $type8_CO_marker_effective_index_end;
			    $stack{'affected_spores'} = $type8_CO_affected_spores;
			} elsif (($total_CO_affected_spores_count == 4) and ($genotype_switches_count_pattern_sorted ne "1:1:1:1")) {
			    # print "total_CO_affected_spores_count=$total_CO_affected_spores_count, genotype_switches_count_pattern_sorted=$genotype_switches_count_pattern_sorted\n";
			    # recategorize as Type9_CO + Type16_GC
			    my $type9_CO_id = "$stack{'id'}.1";
			    my $type9_CO_type = "CO";
			    my $type9_CO_subtype = "Type9_CO"; # double CO with complex GC pattern
			    my $type9_CO_chr = $chr;
			    my $type9_CO_block_start = $stack{'block_start'};
			    my $type9_CO_block_end = $block_end;
			    my $type9_CO_marker_pos_start = $$linkage_blocks_hashref{$type9_CO_chr}{$type9_CO_block_start}{'marker_pos_end'};
			    my $type9_CO_marker_pos_end = $$linkage_blocks_hashref{$type9_CO_chr}{$type9_CO_block_end}{'marker_pos_start'};
			    my $type9_CO_marker_raw_index_start = $$linkage_blocks_hashref{$type9_CO_chr}{$type9_CO_block_start}{'marker_raw_index_end'};
			    my $type9_CO_marker_raw_index_end = $$linkage_blocks_hashref{$type9_CO_chr}{$type9_CO_block_end}{'marker_raw_index_start'};
			    my $type9_CO_marker_effective_index_start = $$linkage_blocks_hashref{$type9_CO_chr}{$type9_CO_block_start}{'marker_effective_index_end'};
			    my $type9_CO_marker_effective_index_end = $$linkage_blocks_hashref{$type9_CO_chr}{$type9_CO_block_end}{'marker_effective_index_start'};
			    my $type9_CO_affected_spores = "1:1:1:1";
			    print "type9_CO_id=$type9_CO_id, type9_CO_chr=$type9_CO_chr, type9_CO_block_start=$type9_CO_block_start, type9_CO_marker_pos_start=$type9_CO_marker_pos_start\n";
			    @{$$CO_associated_GC_hashref{$type9_CO_id}} = ();
			    
			    my $type16_GC_id = "$stack{'id'}.2";
			    my $type16_GC_type = "GC";
			    my $type16_GC_subtype = "Type16_GC"; # complex GC associated with Type9_CO
			    push @{$$CO_associated_GC_hashref{$type9_CO_id}}, $type16_GC_id;
			    my $type16_GC_block_start = $type9_CO_block_start + 1;
			    my $type16_GC_block_end = $type9_CO_block_end - 1;
			    my $type16_GC_chr = $type9_CO_chr;
			    my $type16_GC_marker_pos_start = $$linkage_blocks_hashref{$type16_GC_chr}{$type16_GC_block_start}{'marker_pos_start'};
			    my $type16_GC_marker_pos_end = $$linkage_blocks_hashref{$type16_GC_chr}{$type16_GC_block_end}{'marker_pos_end'};
			    my $type16_GC_marker_raw_index_start = $$linkage_blocks_hashref{$type16_GC_chr}{$type16_GC_block_start}{'marker_raw_index_start'};
			    my $type16_GC_marker_raw_index_end = $$linkage_blocks_hashref{$type16_GC_chr}{$type16_GC_block_end}{'marker_raw_index_end'};
			    my $type16_GC_marker_effective_index_start = $$linkage_blocks_hashref{$type16_GC_chr}{$type16_GC_block_start}{'marker_effective_index_start'};
			    my $type16_GC_marker_effective_index_end = $$linkage_blocks_hashref{$type16_GC_chr}{$type16_GC_block_end}{'marker_effective_index_end'};
			    my @type16_GC_affected_spores = (1, 1, 1, 1);
			    my @genotype_switches_count_pattern = split /:/, $genotype_switches_count_pattern;
			    for (my $i = 0; $i < 4; $i++) {
				if ($genotype_switches_count_pattern[$i] > 1) {
				    $type16_GC_affected_spores[$i] = "1*";
				}
			    }
			    my $type16_GC_affected_spores = join ":", @type16_GC_affected_spores;
			    $merging_log = "$stack{'id'} $stack{'subtype'} $stack{'chr'} $stack{'block_start'} $stack{'block_end'} $stack{'marker_pos_start'} $stack{'marker_pos_end'}\n";
			    $merging_log .= "+\n";
			    $merging_log .= "$id $subtype $chr $block_start $block_end $marker_pos_start $marker_pos_end\n";
			    $merging_log .= "=> new $type9_CO_id ";
			    
			    $merging_log .= "$type9_CO_subtype $type9_CO_chr $type9_CO_block_start $type9_CO_block_end $type9_CO_marker_pos_start $type9_CO_marker_pos_end\n\n";
			    
			    $merging_log .= "\+ new $type16_GC_id ";
			    $merging_log .= "$type16_GC_subtype $type16_GC_chr $type16_GC_block_start $type16_GC_block_end $type16_GC_marker_pos_start $type16_GC_marker_pos_end\n\n";
			    
			    push @$merging_log_arrayref, $merging_log;
			    
			    # add the newly merged event
			    $$events_hashref{$type9_CO_id}{'type'} = $type9_CO_type;
			    $$events_hashref{$type9_CO_id}{'subtype'} = $type9_CO_subtype;
			    $$events_hashref{$type9_CO_id}{'chr'} = $type9_CO_chr;
			    $$events_hashref{$type9_CO_id}{'block_start'} = $type9_CO_block_start;
			    $$events_hashref{$type9_CO_id}{'block_end'} = $type9_CO_block_end;
			    $$events_hashref{$type9_CO_id}{'marker_pos_start'} = $type9_CO_marker_pos_start;
			    $$events_hashref{$type9_CO_id}{'marker_pos_end'} = $type9_CO_marker_pos_end;
			    $$events_hashref{$type9_CO_id}{'marker_raw_index_start'} = $type9_CO_marker_raw_index_start;
			    $$events_hashref{$type9_CO_id}{'marker_raw_index_end'} = $type9_CO_marker_raw_index_end;
			    $$events_hashref{$type9_CO_id}{'marker_effective_index_start'} = $type9_CO_marker_effective_index_start;
			    $$events_hashref{$type9_CO_id}{'marker_effective_index_end'} = $type9_CO_marker_effective_index_end;
			    $$events_hashref{$type9_CO_id}{'affected_spores'} = $type9_CO_affected_spores;
			    
			    $$events_hashref{$type16_GC_id}{'type'} = $type16_GC_type;
			    $$events_hashref{$type16_GC_id}{'subtype'} = $type16_GC_subtype;
			    $$events_hashref{$type16_GC_id}{'chr'} = $type16_GC_chr;
			    $$events_hashref{$type16_GC_id}{'block_start'} = $type16_GC_block_start;
			    $$events_hashref{$type16_GC_id}{'block_end'} = $type16_GC_block_end;
			    $$events_hashref{$type16_GC_id}{'marker_pos_start'} = $type16_GC_marker_pos_start;
			    $$events_hashref{$type16_GC_id}{'marker_pos_end'} = $type16_GC_marker_pos_end;
			    $$events_hashref{$type16_GC_id}{'marker_raw_index_start'} = $type16_GC_marker_raw_index_start;
			    $$events_hashref{$type16_GC_id}{'marker_raw_index_end'} = $type16_GC_marker_raw_index_end;
			    $$events_hashref{$type16_GC_id}{'marker_effective_index_start'} = $type16_GC_marker_effective_index_start;
			    $$events_hashref{$type16_GC_id}{'marker_effective_index_end'} = $type16_GC_marker_effective_index_end;
			    $$events_hashref{$type16_GC_id}{'affected_spores'} = $type16_GC_affected_spores;
			    
			    # print "delete CO $stack{'id'}\n";
			    delete $$events_hashref{$stack{'id'}};
			    # print "delete CO $id\n";
			    delete $$events_hashref{$id};
			    foreach my $GC_id (@{$$CO_associated_GC_hashref{$stack{'id'}}}) {
				# print "deleted GC associated with the CO $stack{'id'}: $GC_id\n";
				delete $$events_hashref{$GC_id};
			    }
			    foreach my $GC_id (@{$$CO_associated_GC_hashref{$id}}) {
				# print "deleted GC associated with the CO $id: $GC_id\n";
				delete $$events_hashref{$GC_id};
			    }
			    delete $$CO_associated_GC_hashref{$stack{'id'}};
			    delete $$CO_associated_GC_hashref{$id};
			    
			    # reset stack
			    $stack{'id'} = $type9_CO_id;
			    $stack{'type'} = $type9_CO_type;
			    $stack{'subtype'} = $type9_CO_subtype;
			    $stack{'chr'} = $type9_CO_chr;
			    $stack{'block_start'} = $type9_CO_block_start;
			    $stack{'block_end'} = $type9_CO_block_end;
			    $stack{'marker_pos_start'} = $type9_CO_marker_pos_start;
			    $stack{'marker_pos_end'} = $type9_CO_marker_pos_end;
			    $stack{'marker_raw_index_start'} = $type9_CO_marker_raw_index_start;
			    $stack{'marker_raw_index_end'} = $type9_CO_marker_raw_index_end;
			    $stack{'marker_effective_index_start'} = $type9_CO_marker_effective_index_start;
			    $stack{'marker_effective_index_end'} = $type9_CO_marker_effective_index_end;
			    $stack{'affected_spores'} = $type9_CO_affected_spores;
			}
		    }
		}
	    }
	}
    }
}


sub find_closest_event {
    my ($query_type, $query_chr, $query_marker_pos_start, $query_marker_pos_end, $events_hashref, $linkage_blocks_hashref) = @_;
    my %closest_event = ();
    $closest_event{'distance_to_query'} = 9999999999;
    foreach my $id (sort custom_sort keys %$events_hashref) {
    # foreach my $id (sort {(($a =~ /(\d+\.?\d+?)/)[0] || 0) <=> (($b =~ /(\d+\.?\d+?)/)[0] || 0)} keys %$events_hashref) {
        my $type = $$events_hashref{$id}{'type'};
        my $subtype = $$events_hashref{$id}{'subtype'};
        my $chr = $$events_hashref{$id}{'chr'};
        my $block_start = $$events_hashref{$id}{'block_start'};
        my $block_end = $$events_hashref{$id}{'block_end'};
        my $marker_pos_start = $$events_hashref{$id}{'marker_pos_start'};
        my $marker_pos_end = $$events_hashref{$id}{'marker_pos_end'};
        my $affected_spores = $$events_hashref{$id}{'affected_spores'};
        if (($chr eq $query_chr) and ($type eq $query_type)) {
	    if (($marker_pos_start <= $query_marker_pos_start) and ($marker_pos_end >= $query_marker_pos_end)) {
		# the query region is contained in the current event
		$closest_event{'id'} = $id;
		$closest_event{'type'} = $type;
		$closest_event{'subtype'} = $subtype;
		$closest_event{'chr'} = $chr;
		$closest_event{'block_start'} = $block_start;
		$closest_event{'block_end'} = $block_end;
		$closest_event{'marker_pos_start'} = $marker_pos_start;
		$closest_event{'marker_pos_end'} = $marker_pos_end;
		$closest_event{'affected_spores'} = $affected_spores;
		$closest_event{'distance_to_query'} = -1;
		last;
	    } elsif ($marker_pos_end <= $query_marker_pos_start) {
		# query is on the right
		my $distance_to_query = $query_marker_pos_start - $marker_pos_end;
		if($distance_to_query < $closest_event{'distance_to_query'}) {
		    $closest_event{'id'} = $id;
		    $closest_event{'type'} = $type;
		    $closest_event{'subtype'} = $subtype;
		    $closest_event{'chr'} = $chr;
		    $closest_event{'block_start'} = $block_start;
		    $closest_event{'block_end'} = $block_end;
		    $closest_event{'marker_pos_start'} = $marker_pos_start;
		    $closest_event{'marker_pos_end'} = $marker_pos_end;
		    $closest_event{'affected_spores'} = $affected_spores;
		    $closest_event{'distance_to_query'} = $distance_to_query;
		}
	    } elsif ($marker_pos_start >=  $query_marker_pos_end) {
		# query is on the left
		my $distance_to_query = $marker_pos_start - $query_marker_pos_end;
		if ($distance_to_query < $closest_event{'distance_to_query'}) {
		    $closest_event{'id'} = $id;
		    $closest_event{'type'} = $type;
		    $closest_event{'subtype'} = $subtype;
		    $closest_event{'chr'} = $chr;
		    $closest_event{'block_start'} = $block_start;
		    $closest_event{'block_end'} = $block_end;
		    $closest_event{'marker_pos_start'} = $marker_pos_start;
		    $closest_event{'marker_pos_end'} = $marker_pos_end;
		    $closest_event{'affected_spores'} = $affected_spores;
		    $closest_event{'distance_to_query'} = $distance_to_query;
		}
	    } else {
		# query has partial overlap with tested event
		$closest_event{'id'} = $id;
		$closest_event{'type'} = $type;
		$closest_event{'subtype'} = $subtype;
		$closest_event{'chr'} = $chr;
		$closest_event{'block_start'} = $block_start;
		$closest_event{'block_end'} = $block_end;
		$closest_event{'marker_pos_start'} = $marker_pos_start;
		$closest_event{'marker_pos_end'} = $marker_pos_end;
		$closest_event{'affected_spores'} = $affected_spores;
		$closest_event{'distance_to_query'} = -2;
	    }
	}
    }
    return %closest_event;
} 



sub custom_sort {
    my @a = split /\./, $a;
    my @b = split /\./, $b;
    my $a_array_length = scalar @a;
    my $b_array_length = scalar @b;
    my $array_length;
    if ($a_array_length >= $b_array_length) {
	$array_length = $a_array_length;
    } else {
	$array_length = $b_array_length;
    }

    for (my $i = 0; $i < $array_length; $i++) {
	if (not exists $a[$i]) {
	    $a[$i] = 0;
	}
	if (not exists $b[$i]) {
	    $b[$i] = 0;
	}
	if ($a[$i] > $b[$i]) {
	    return 1;
	} elsif ($a[$i] < $b[$i]) {
            return -1;
	} else {
	    next;
	}
    }
}


