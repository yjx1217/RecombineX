#!/usr/bin/perl
use warnings FATAL => 'all';
use strict;
use Getopt::Long;
use Math::Random qw(:all);
use Math::Round;
use List::Util qw(shuffle);
use Data::Dumper;

##############################################################
#  script: simulate_recombinant_tetrad.pl 
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2020.01.19
#  description: simulating tetrad genotype after meiotic recombination based on input reference genome and marker table
#  example: perl simulate_recombinant_tetrad.pl -g ref.genome.raw.relabel.fa(.gz) -g1 S288C -g2 SK1 -m S288C-SK1.ref.final.SNP.markers.txt(.gz) -co_num 90 -nco_num 60 -gc_by_co_ratio 1 -co_gc_size_mean 2500 -co_gc_size_stdev 2000 -nco_gc_size_mean 2200 -nco_gc_size_stdev 2200 -d 10000 -l linked_regions.forSim.txt -m 0 -p Sim 
##############################################################

################
# input parameters
################
my $coordinate_genome; # The input genome file for set up the coordinate system (FASTA format)
my $parent1_tag = "P1"; # The genome tag for simulated parent1 genome (P1). Default = "P1".
my $parent2_tag = "P2"; # The genome tag for simulated parent2 genome (P2). Default = "P2".
my $marker_table_file; # SNP markers based on the input genome.
my $co_num = 90; # The number of crossover (CO) per tetrad. Default = 90. # (mean = 90.5 for S. cerevisiae based on Mancera et al. 2008 Nature paper).
my $nco_num = 60; # The number of noncrossover (NCO) per tetrad. Default = 60. # (mean = 66.1 for S. cerevisiae based on Mancera et al. 2008 Nature paper).
my $min_interevent_distance = 10000; # The minimal genomic distance in basepairs (bp) allowed between two independent recombination events. Default = 10000. (i.e. 10 kb). 
my $linked_region_bed; # genomic regions for complete linkeage (i.e. no recombination) in 3-column BED format (without header).
my $gc_by_co_ratio = 1; # the ratio of the number of crossover with associated gene conversions (GC) among all crossover events. 
my $co_gc_size_mean = 2500; # The mean of crossover-associated gene conversion (CO-GC) size in basepairs (bp). Default = 2500. ("2461.836" for S. cerevisiae based on Mancera et al. 2008 Nature paper).
my $co_gc_size_stdev = 2000; # The standard deviation of crossover-associated gene conversion (CO-GC) size in basepair (bp). Default = 2000. ("2057.673" for S. cerevisiae based on Mancera et al. 2008 Nature paper).
my $nco_gc_size_mean = 2250; # The mean of noncrossover-associated gene conversion (NCO-GC) size in basepairs (bp). Default = 2250. ("2247.716" for NCO-GC for S. cerevisiae based on Mancera et al. 2008 Nature paper).
my $nco_gc_size_stdev=2200; # The standard deviation of noncrossover-associated gene conversion (NCO-GC) size in basepair (bp). Default = 2200. ("2173.46" for NCO-GC for S. cerevisiae based on Mancera et al. 2008 Nature paper).
my $min_gc_size=100; # Safe lower bound gene conversion (GC) size in basepair (bp). Default = 100. 
my $max_gc_size=5000; # Safe upper bound gene conversion (GC) size in basepair (bp). Default = 5000. 
my $output_prefix = "simulated_tetrad"; # prefix for outputs
my $random_seed = "simulated_tetrad"; # random seed phrase

my $mut_num = 0; # mutation number per tetrad


GetOptions('cg|coordinate_genome_fasta:s' => \$coordinate_genome,
	   'p1|parent1_tag:s' => \$parent1_tag,
	   'p2|parent2_tag:s' => \$parent2_tag,
	   'm|markers:s' => \$marker_table_file,
	   'co|co_num:i' => \$co_num,
	   'nco|nco_num:i' => \$nco_num,
	   'gc_by_co_ratio|gc_by_co_ratio:f' => \$gc_by_co_ratio,
	   'co_gc_size_mean|co_gc_size_mean:f' => \$co_gc_size_mean,
	   'co_gc_size_stdev|co_gc_size_stdev:f' => \$co_gc_size_stdev,
	   'nco_gc_size_mean|nco_gc_size_mean:f' => \$nco_gc_size_mean,
	   'nco_gc_size_stdev|nco_gc_size_stdev:f' => \$nco_gc_size_stdev,
	   'min_gc_size|min_gc_size:f' => \$min_gc_size,
	   'max_gc_size|max_gc_size:f' => \$max_gc_size,
	   'd|min_interevent_distance:i' => \$min_interevent_distance,
	   'l|linked_region_bed:s' => \$linked_region_bed,
	   'mut_num|mut_num:i' => \$mut_num,
	   'p|prefix:s' => \$output_prefix,
	   's|seed:s' => \$random_seed);

# initialize the seed for random number generator
random_set_seed_from_phrase($random_seed);

print "\nThis simulation use the random seed: $random_seed\n\n";

my $coordinate_genome_fh = read_file($coordinate_genome);
my %coordinate_genome = ();
my @coordinate_genome = ();
parse_fasta_file($coordinate_genome_fh, \%coordinate_genome, \@coordinate_genome);
close $coordinate_genome_fh;

my $markers_fh = read_file($marker_table_file);
my %markers = ();
parse_markers_table_file($markers_fh, \%markers);
close $markers_fh;

# reconstruct the genome of the two parents based on SNP markers
my %parent1_genome = %coordinate_genome;
my %parent2_genome = %coordinate_genome;
modify_genome_based_on_markers(\%coordinate_genome, \%parent1_genome, \%parent2_genome, \%markers);

my $tetrad_id = "simulated_tetrad";
my @spores = qw(a b c d);
my %initial_genotype_pattern = (
    'a' => $parent1_tag, 
    'b' => $parent1_tag, 
    'c' => $parent2_tag, 
    'd' => $parent2_tag,
    );
my %spore_id2index = (
    'a' => 0,
    'b' => 1,
    'c' => 2,
    'd' => 3,
    );
my %spore_index2id = (
    0 => 'a',
    1 => 'b',
    2 => 'c',
    3 => 'd',
    );

my @initial_genotype_pattern = ($initial_genotype_pattern{'a'}, $initial_genotype_pattern{'b'}, $initial_genotype_pattern{'c'}, $initial_genotype_pattern{'d'});

# create genome space
my %coordinate_genome_space = create_genome_space(\%coordinate_genome);
my $total_coordinate_genome_size = 0;
foreach my $chr (@coordinate_genome) {
    if ($coordinate_genome_space{$chr}{'end'} > $total_coordinate_genome_size) {
	$total_coordinate_genome_size = $coordinate_genome_space{$chr}{'end'};
    }
}

# get infomration for complete linked regions
my %invalid_recomb_regions = ();
if ($linked_region_bed ne "") {
    my $linked_region_bed_fh = read_file($linked_region_bed);
    parse_linked_regions_file($linked_region_bed_fh, \%invalid_recomb_regions, \%coordinate_genome_space);
}

# simulate recombination
my %recomb_history = ();
my %recomb_breakpoints = ();
foreach my $chr (@coordinate_genome) {
    my $chr_start = 1;
    my $chr_end = length $coordinate_genome{$chr};
    $recomb_breakpoints{$chr}{$chr_start}{'note'} = "chr_start";
    $recomb_breakpoints{$chr}{$chr_end}{'note'} = "chr_end";
}

# simulate CO
my %CO = ();
# sample genomic regions
for (my $i = 1; $i <= $co_num; $i++) {
    my $co_size;
    # determine whether to sample GC tract based on the gc_by_co_ratio
    if (random_uniform(1, 0, 1) <= $gc_by_co_ratio) {
	# having CO-associated GC 
	my $co_gc_size = $min_gc_size - 1;
	while (($co_gc_size < $min_gc_size) or ($co_gc_size > $max_gc_size)) {
	    $co_gc_size = random_normal(1, $co_gc_size_mean, $co_gc_size_stdev);
	    # my $co_gc_size = random_gamma(1, $co_gc_size_gamma_rate, $co_gc_size_gamma_shape) * $co_gc_size_scale;
	}
	$co_gc_size = round($co_gc_size);
	$co_size += $co_gc_size;
    } else {
	# no CO-associated GC
	$co_size = 1;
    }
    my $chr_sampled = sample_chr($total_coordinate_genome_size, \%coordinate_genome_space, \%invalid_recomb_regions);
    my ($interval_start_sampled, $interval_end_sampled) = sample_genomic_interval($co_size, $chr_sampled, \%coordinate_genome_space, \%invalid_recomb_regions);
    print "CO:$i > sampled interval with size of $co_size: $chr_sampled:$interval_start_sampled-$interval_end_sampled\n";
    $CO{$chr_sampled}{$interval_start_sampled}{'end'} = $interval_end_sampled;
    my $key = "${chr_sampled}:${interval_start_sampled}-${interval_end_sampled}";
    $invalid_recomb_regions{$key}{'chr'} = $chr_sampled;
    if (($interval_start_sampled - $min_interevent_distance) > 1) {
	$invalid_recomb_regions{$key}{'start'} = $interval_start_sampled - $min_interevent_distance;
    } else {
	$invalid_recomb_regions{$key}{'start'} = 1;
    }
    if (($interval_end_sampled + $min_interevent_distance) < (length $coordinate_genome{$chr_sampled})) {
	$invalid_recomb_regions{$key}{'end'} = $interval_end_sampled + $min_interevent_distance;
    } else {
	$invalid_recomb_regions{$key}{'end'} = length $coordinate_genome{$chr_sampled};
    }

    my $start_in_genome_space = $coordinate_genome_space{$chr_sampled}{'start'} + $invalid_recomb_regions{$key}{'start'} - 1;
    my $end_in_genome_space = $coordinate_genome_space{$chr_sampled}{'start'} + $invalid_recomb_regions{$key}{'end'} - 1;
    $invalid_recomb_regions{$key}{'start_in_genome_space'} = $start_in_genome_space;
    $invalid_recomb_regions{$key}{'end_in_genome_space'} = $end_in_genome_space;

    $invalid_recomb_regions{$key}{'note'} = "existing CO";

    if ($interval_start_sampled == $interval_end_sampled) {
	$recomb_breakpoints{$chr_sampled}{$interval_start_sampled}{'note'} = "CO_interval_start_and_end";
    } else {
	$recomb_breakpoints{$chr_sampled}{$interval_start_sampled}{'note'} = "CO_interval_start";
	$recomb_breakpoints{$chr_sampled}{$interval_end_sampled}{'note'} = "CO_interval_end";
    }
}

my %preliminary_linkage_blocks_by_coord = ();
my %preliminary_linkage_blocks_by_id = ();
my $preliminary_block_id = 0;
foreach my $chr (sort keys %recomb_breakpoints) {
    my $offset = 1;
    foreach my $bp (sort {$a <=> $b} keys %{$recomb_breakpoints{$chr}}) {
	if ( $bp == 1) {
	    next;
	} else {
	    $preliminary_block_id++;
	    # for simplicity, here we will use closed intervals to represent each linkage blocks (i.e. not including the exact breakpoints)                                                 
	    my $preliminary_block_start = $offset + 1;
	    my $preliminary_block_end = $bp - 1;
	    my $preliminary_block_size = $preliminary_block_end - $preliminary_block_start + 1;
	    if ($offset == 1) {
		$preliminary_linkage_blocks_by_coord{$chr}{$offset}{'genotype_pattern'} = join ":", @initial_genotype_pattern;
		$preliminary_linkage_blocks_by_id{$preliminary_block_id}{'genotype_pattern'} = join ":", @initial_genotype_pattern;
	    }
	    
	    $preliminary_linkage_blocks_by_coord{$chr}{$offset}{'block_id'} = $preliminary_block_id;
	    $preliminary_linkage_blocks_by_coord{$chr}{$offset}{'block_start'} = $preliminary_block_start;
	    $preliminary_linkage_blocks_by_coord{$chr}{$offset}{'chr'} = $chr;
	    $preliminary_linkage_blocks_by_coord{$chr}{$offset}{'block_end'} = $preliminary_block_end;
	    $preliminary_linkage_blocks_by_coord{$chr}{$offset}{'block_size'} = $preliminary_block_size;
	    $preliminary_linkage_blocks_by_coord{$chr}{$offset}{'interval_start'} = $offset;
	    $preliminary_linkage_blocks_by_coord{$chr}{$offset}{'interval_end'} = $bp;
	    
	    
	    $preliminary_linkage_blocks_by_id{$preliminary_block_id}{'chr'} = $chr;
	    $preliminary_linkage_blocks_by_id{$preliminary_block_id}{'block_start'} = $preliminary_block_start;
	    $preliminary_linkage_blocks_by_id{$preliminary_block_id}{'block_end'} = $preliminary_block_end;
	    $preliminary_linkage_blocks_by_id{$preliminary_block_id}{'block_size'} = $preliminary_block_size;
	    $preliminary_linkage_blocks_by_id{$preliminary_block_id}{'interval_start'} = $offset;
	    $preliminary_linkage_blocks_by_id{$preliminary_block_id}{'interval_end'} = $bp;
	    
	    
	    if ($offset eq $bp) {
		$preliminary_linkage_blocks_by_coord{$chr}{$offset}{'interval_type'} = "CO-GC_interval";
		$preliminary_linkage_blocks_by_id{$preliminary_block_id}{'interval_type'} = "CO-GC_interval";
	    } else {
		$preliminary_linkage_blocks_by_coord{$chr}{$offset}{'interval_type'} = "CO_interval";
		$preliminary_linkage_blocks_by_id{$preliminary_block_id}{'interval_type'} = "CO_interval";
	    }
	    $offset = $bp;
	}
    }
}

# simulate CO associated GC
my $CO_id = 0;
my $GC_id = 0;
foreach my $chr (sort keys %CO) {
    my %genotype_pattern_before_CO = %initial_genotype_pattern;
    my @CO_interval_start_sampled = sort {$a <=> $b} keys %{$CO{$chr}};
    foreach my $CO_interval_start (@CO_interval_start_sampled) {
	my $CO_interval_end = $CO{$chr}{$CO_interval_start}{'end'};
	if ($CO_interval_start eq $CO_interval_end) {
	    # CO with no GC tract
	    my $preliminary_block_id_after_CO = $preliminary_linkage_blocks_by_coord{$chr}{$CO_interval_start}{'block_id'};
	    my $preliminary_block_id_before_CO = $preliminary_block_id_after_CO - 1;
	    my @genotype_pattern_before_CO = split /:/, $preliminary_linkage_blocks_by_id{$preliminary_block_id_before_CO}{'genotype_pattern'};
	    my %genotype_pattern_before_CO = ();
	    for(my $i = 0; $i < 4; $i++) {
		my $spore = $spores[$i];
		$genotype_pattern_before_CO{$spore} = $genotype_pattern_before_CO[$i];
	    }
	    # simulate genotype switch
	    my ($donor_spore) = random_sample_from_array(1, \@spores);
	    my $donor_spore_index = $spore_id2index{$donor_spore};
	    print "donor_spore: $donor_spore, donor_spore_index: $donor_spore_index\n";
	    my $receipent_spore;
	    my @available_recipient_spores = ();
	    foreach my $spore (@spores) {
		if (($spore ne $donor_spore) and ($genotype_pattern_before_CO{$spore} ne $genotype_pattern_before_CO{$donor_spore})) {
		    push @available_recipient_spores, $spore;
		}
	    }
	    my ($recipient_spore) = random_sample_from_array(1, \@available_recipient_spores);
	    my $recipient_spore_index = $spore_id2index{$recipient_spore};
	    print "recipient_spore: $recipient_spore, recipient_spore_index: $recipient_spore_index\n";

	    my @co_affected_spores = qw(0 0 0 0);
	    $co_affected_spores[$donor_spore_index] = 1;
	    $co_affected_spores[$recipient_spore_index] = 1;
	    my $co_affected_spores = join ":", @co_affected_spores;

	    my %genotype_pattern_after_CO = ();
	    # switch genotypes for the donor and recipient
	    my @genotype_pattern_after_CO = @genotype_pattern_before_CO;
	    ($genotype_pattern_after_CO{$donor_spore}, $genotype_pattern_after_CO{$recipient_spore}) = ($genotype_pattern_before_CO{$recipient_spore}, $genotype_pattern_before_CO{$donor_spore});
	    $genotype_pattern_after_CO[$donor_spore_index] = $genotype_pattern_after_CO{$donor_spore};
	    $genotype_pattern_after_CO[$recipient_spore_index] = $genotype_pattern_after_CO{$recipient_spore};
	    $preliminary_linkage_blocks_by_id{$preliminary_block_id_after_CO}{'genotype_pattern'} = join ":", @genotype_pattern_after_CO;
	    $preliminary_linkage_blocks_by_coord{$chr}{$CO_interval_start}{'genotype_pattern'} = join ":", @genotype_pattern_after_CO;
	    
	    $CO_id++;
	    $recomb_history{'CO'}{$CO_id}{'co_affected_spores'} = $co_affected_spores;
	    $recomb_history{'CO'}{$CO_id}{'recomb_desc'} = "CO,Type1_CO,$co_affected_spores";
	    $recomb_history{'CO'}{$CO_id}{'chr'} = $chr;
	    $recomb_history{'CO'}{$CO_id}{'start'} = $CO_interval_start;
	    $recomb_history{'CO'}{$CO_id}{'end'} = $CO_interval_end;
	    $recomb_history{'CO'}{$CO_id}{'subtype'} = "Type1_CO";

	} else {
	    # CO with GC tract
	    print "chr=$chr, CO_interval_start=$CO_interval_start\n";
	    my $preliminary_block_id_at_CO = $preliminary_linkage_blocks_by_coord{$chr}{$CO_interval_start}{'block_id'};
	    my $preliminary_block_id_before_CO = $preliminary_block_id_at_CO - 1;
	    my $preliminary_block_id_after_CO = $preliminary_block_id_at_CO + 1;

	    print "preliminary_block_id_before_CO=$preliminary_block_id_before_CO\n";
	    print "preliminary_block_id_at_CO=$preliminary_block_id_at_CO\n";
	    print "preliminary_block_id_after_CO=$preliminary_block_id_after_CO\n";

	    my $genotype_pattern_before_CO = $preliminary_linkage_blocks_by_id{$preliminary_block_id_before_CO}{'genotype_pattern'};
	    my @genotype_pattern_before_CO = split /:/, $genotype_pattern_before_CO;
	    my %genotype_pattern_before_CO = ();
	    print "genotype_pattern_before_CO=$genotype_pattern_before_CO\n";

	    for(my $i = 0; $i < 4; $i++) {
		my $spore = $spores[$i];
		$genotype_pattern_before_CO{$spore} = $genotype_pattern_before_CO[$i];
	    }

	    # simulate genotype switch
	    my ($donor_spore) = random_sample_from_array(1, \@spores);
	    print "donor_spore: $donor_spore\n";
	    my $donor_spore_index = $spore_id2index{$donor_spore};
	    print "donor_spore: $donor_spore, donor_spore_index: $donor_spore_index\n";
	    my $receipent_spore;
	    my @available_recipient_spores = ();
	    foreach my $spore (@spores) {
		if (($spore ne $donor_spore) and ($genotype_pattern_before_CO{$spore} ne $genotype_pattern_before_CO{$donor_spore})) {
		    push @available_recipient_spores, $spore;
		}
	    }
	    my ($recipient_spore) = random_sample_from_array(1, \@available_recipient_spores);
	    my $recipient_spore_index = $spore_id2index{$recipient_spore};
	    print "recipient_spore: $recipient_spore, recipient_spore_index: $recipient_spore_index\n";

	    my @co_affected_spores = qw(0 0 0 0);
	    $co_affected_spores[$donor_spore_index] = 1;
	    $co_affected_spores[$recipient_spore_index] = 1;
	    my $co_affected_spores = join ":", @co_affected_spores;
	    my %genotype_pattern_after_CO = ();
	    # switch genotypes for the donor and recipient
	    ($genotype_pattern_after_CO{$donor_spore}, $genotype_pattern_after_CO{$recipient_spore}) = ($genotype_pattern_before_CO{$recipient_spore}, $genotype_pattern_before_CO{$donor_spore});
	    my @genotype_pattern_after_CO = @genotype_pattern_before_CO;
	    $genotype_pattern_after_CO[$donor_spore_index] = $genotype_pattern_after_CO{$donor_spore};
	    $genotype_pattern_after_CO[$recipient_spore_index] = $genotype_pattern_after_CO{$recipient_spore};
	    my $genotype_pattern_after_CO = join ":", @genotype_pattern_after_CO;
	    print "genotype_pattern_after_CO=$genotype_pattern_after_CO\n";
	    $preliminary_linkage_blocks_by_id{$preliminary_block_id_after_CO}{'genotype_pattern'} = $genotype_pattern_after_CO;
	    $preliminary_linkage_blocks_by_coord{$chr}{$CO_interval_end}{'genotype_pattern'} = $genotype_pattern_after_CO;

	    my @genotype_pattern_at_CO = @genotype_pattern_before_CO;
	    # flip a coin to decide which genotype to use for the GC tract
	    if (random_uniform(1, 0, 1) < 0.5) {
		$genotype_pattern_at_CO[$donor_spore_index] = $genotype_pattern_before_CO[$recipient_spore_index];
	    } else {
		$genotype_pattern_at_CO[$recipient_spore_index] = $genotype_pattern_before_CO[$donor_spore_index];
	    }
	    my $genotype_pattern_at_CO = join ":", @genotype_pattern_at_CO;
	    print "genotype_pattern_at_CO=$genotype_pattern_at_CO\n";
	    $preliminary_linkage_blocks_by_id{$preliminary_block_id_at_CO}{'genotype_pattern'} = $genotype_pattern_at_CO;
	    $preliminary_linkage_blocks_by_coord{$chr}{$CO_interval_start}{'genotype_pattern'} = $genotype_pattern_at_CO;


	    $CO_id++;
	    $recomb_history{'CO'}{$CO_id}{'co_affected_spores'} = $co_affected_spores;
	    $recomb_history{'CO'}{$CO_id}{'gc_affected_spores'} = $co_affected_spores;
	    $recomb_history{'CO'}{$CO_id}{'recomb_desc'} = "CO,Type2_CO,$co_affected_spores";
	    $recomb_history{'CO'}{$CO_id}{'chr'} = $chr;
	    $recomb_history{'CO'}{$CO_id}{'start'} = $CO_interval_start;
	    $recomb_history{'CO'}{$CO_id}{'end'} = $CO_interval_end;
	    $recomb_history{'CO'}{$CO_id}{'subtype'} = "Type2_CO";

	    $GC_id++;
	    $recomb_history{'GC'}{$GC_id}{'gc_affected_spores'} = $co_affected_spores;
	    $recomb_history{'GC'}{$GC_id}{'recomb_desc'} = "GC,Type2_GC,$co_affected_spores";
	    $recomb_history{'GC'}{$GC_id}{'chr'} = $chr;
	    $recomb_history{'GC'}{$GC_id}{'start'} = $CO_interval_start;
	    $recomb_history{'GC'}{$GC_id}{'end'} = $CO_interval_end;
	    $recomb_history{'GC'}{$GC_id}{'subtype'} = "Type2_GC";

	}
    }
}


# simulate NCO
my %NCO = ();
for (my $i = 1; $i <= $nco_num; $i++) {
    my $nco_size = $min_gc_size - 1;
    while (($nco_size < $min_gc_size) or ($nco_size > $max_gc_size)) {
	$nco_size = random_normal(1, $nco_gc_size_mean, $nco_gc_size_stdev);
	# my $nco_size = random_gamma(1, $nco_gc_size_gamma_rate, $nco_gc_size_gamma_shape) * $nco_gc_size_scale;
    }
    $nco_size = round($nco_size);
    my $chr_sampled = sample_chr($total_coordinate_genome_size, \%coordinate_genome_space, \%invalid_recomb_regions);
    my ($interval_start_sampled, $interval_end_sampled) = sample_genomic_interval($nco_size, $chr_sampled, \%coordinate_genome_space, \%invalid_recomb_regions);
    print "NCO:$i > sampled interval with size of $nco_size: $chr_sampled:$interval_start_sampled-$interval_end_sampled\n";
    $NCO{$chr_sampled}{$interval_start_sampled}{'end'} = $interval_end_sampled;
    my $key = "${chr_sampled}:${interval_start_sampled}-${interval_end_sampled}";
    $invalid_recomb_regions{$key}{'chr'} = $chr_sampled;
    $invalid_recomb_regions{$key}{'start'} = $interval_start_sampled - $min_interevent_distance;
    $invalid_recomb_regions{$key}{'end'} = $interval_end_sampled + $min_interevent_distance;

    my $start_in_genome_space = $coordinate_genome_space{$chr_sampled}{'start'} + $invalid_recomb_regions{$key}{'start'} - 1;
    my $end_in_genome_space = $coordinate_genome_space{$chr_sampled}{'start'} + $invalid_recomb_regions{$key}{'end'} - 1;
    $invalid_recomb_regions{$key}{'start_in_genome_space'} = $start_in_genome_space;
    $invalid_recomb_regions{$key}{'end_in_genome_space'} = $end_in_genome_space;

    $invalid_recomb_regions{$key}{'note'} = "existing NCO";
    
    $recomb_breakpoints{$chr_sampled}{$interval_start_sampled}{'note'} = "NCO_interval_start";
    $recomb_breakpoints{$chr_sampled}{$interval_end_sampled}{'note'} = "NCO_interval_end";

}

my $NCO_subblock_id = 0;
my %NCO_interval_end2CO_master_interval_map = ();
foreach my $chr (sort keys %NCO) {
    foreach my $NCO_interval_start (sort {$a <=> $b} keys %{$NCO{$chr}}) {
	my $NCO_interval_end = $NCO{$chr}{$NCO_interval_start}{'end'};
	my $NCO_block_start = $NCO_interval_start + 1;
	my $NCO_block_end = $NCO_interval_end - 1;
	my $NCO_block_size = $NCO_interval_end - $NCO_block_start + 1;
	foreach my $offset (sort {$a <=> $b} keys %{$preliminary_linkage_blocks_by_coord{$chr}}) {
	    if ($preliminary_linkage_blocks_by_coord{$chr}{$offset}{'interval_type'} eq "CO_interval") {
		my $preliminary_block_start = $preliminary_linkage_blocks_by_coord{$chr}{$offset}{'block_start'};
		my $preliminary_block_end = $preliminary_linkage_blocks_by_coord{$chr}{$offset}{'block_end'};
		my $genotype_pattern = $preliminary_linkage_blocks_by_coord{$chr}{$offset}{'genotype_pattern'};
		if ((($preliminary_block_start <= $NCO_block_start) and ($preliminary_block_end >= $NCO_block_start)) and (($preliminary_block_start <= $NCO_block_end) and ($preliminary_block_end >= $NCO_block_end))) {
		    my $preliminary_block_id = $preliminary_linkage_blocks_by_coord{$chr}{$offset}{'block_id'};
		    $NCO_subblock_id++; 
		    $preliminary_linkage_blocks_by_coord{$chr}{$offset}{'NCO_subblock'}{$NCO_subblock_id}{'block_id'} = $NCO_subblock_id;
		    $preliminary_linkage_blocks_by_coord{$chr}{$offset}{'NCO_subblock'}{$NCO_subblock_id}{'block_start'} = $NCO_block_start;
		    $preliminary_linkage_blocks_by_coord{$chr}{$offset}{'NCO_subblock'}{$NCO_subblock_id}{'block_end'} = $NCO_block_end;
		    $preliminary_linkage_blocks_by_coord{$chr}{$offset}{'NCO_subblock'}{$NCO_subblock_id}{'interval_start'} = $NCO_interval_start;
		    $preliminary_linkage_blocks_by_coord{$chr}{$offset}{'NCO_subblock'}{$NCO_subblock_id}{'interval_end'} = $NCO_interval_end;
		    $preliminary_linkage_blocks_by_coord{$chr}{$offset}{'NCO_subblock'}{$NCO_subblock_id}{'interval_type'} = "NCO_interval";
		    
		    $preliminary_linkage_blocks_by_coord{$chr}{$NCO_interval_start}{'block_id'} = "$preliminary_block_id.$NCO_subblock_id";
		    $preliminary_linkage_blocks_by_coord{$chr}{$NCO_interval_start}{'block_start'} = $NCO_block_start;
		    $preliminary_linkage_blocks_by_coord{$chr}{$NCO_interval_start}{'block_end'} = $NCO_block_end;
		    $preliminary_linkage_blocks_by_coord{$chr}{$NCO_interval_start}{'interval_start'} = $NCO_interval_start;
		    $preliminary_linkage_blocks_by_coord{$chr}{$NCO_interval_start}{'interval_end'} = $NCO_interval_end;
		    $preliminary_linkage_blocks_by_coord{$chr}{$NCO_interval_start}{'interval_type'} = "NCO_interval";
		    
		    $NCO_interval_end2CO_master_interval_map{$chr}{$NCO_interval_end}{'CO_masterblock'}{'block_id'} = $preliminary_block_id;
		    $NCO_interval_end2CO_master_interval_map{$chr}{$NCO_interval_end}{'CO_masterblock'}{'interval_start'} = $offset;
		    $NCO_interval_end2CO_master_interval_map{$chr}{$NCO_interval_end}{'CO_masterblock'}{'interval_end'} = $preliminary_linkage_blocks_by_coord{$chr}{$offset}{'interval_end'};
		    $NCO_interval_end2CO_master_interval_map{$chr}{$NCO_interval_end}{'CO_masterblock'}{'genotype_pattern'} = $genotype_pattern;
		    
		    my @genotype_pattern_before_NCO = split /:/, $genotype_pattern;
		    my %genotype_pattern_before_NCO = ();
		    my %genotype_pattern_at_NCO = ();
		    my @genotype_pattern_at_NCO = @genotype_pattern_before_NCO;
		    
		    for(my $i = 0; $i < 4; $i++) {
			my $spore = $spores[$i];
			$genotype_pattern_at_NCO{$spore} = $genotype_pattern_before_NCO[$i];
		    }
		    
		    # simulate genotype switch
		    my ($donor_spore) = random_sample_from_array(1, \@spores);
		    my $donor_spore_index = $spore_id2index{$donor_spore};
		    print "donor_spore: $donor_spore, donor_spore_index: $donor_spore_index\n";
		    my $receipent_spore;
		    my @available_recipient_spores = ();
		    foreach my $spore (@spores) {
			if (($spore ne $donor_spore) and ($genotype_pattern_at_NCO{$spore} ne $genotype_pattern_at_NCO{$donor_spore})) {
			    push @available_recipient_spores, $spore;
			}
		    }
		    my ($recipient_spore) = random_sample_from_array(1, \@available_recipient_spores);
		    my $recipient_spore_index = $spore_id2index{$recipient_spore};
		    print "recipient_spore: $recipient_spore, recipient_spore_index: $recipient_spore_index\n";
		    
		    my @gc_affected_spores = qw(0 0 0 0);
		    $gc_affected_spores[$recipient_spore_index] = 1;
		    my $gc_affected_spores = join ":", @gc_affected_spores;
		    # switch genotypes for the donor and recipient
		    $genotype_pattern_at_NCO{$recipient_spore} = $genotype_pattern_at_NCO{$donor_spore};
		    $genotype_pattern_at_NCO[$recipient_spore_index] = $genotype_pattern_at_NCO{$donor_spore};
		    $preliminary_linkage_blocks_by_coord{$chr}{$offset}{'NCO_subblock'}{$NCO_subblock_id}{'genotype_pattern'} = join ":", @genotype_pattern_at_NCO;
		    $preliminary_linkage_blocks_by_coord{$chr}{$NCO_interval_start}{'genotype_pattern'} = join ":", @genotype_pattern_at_NCO;
		    
		    $GC_id++;
		    $recomb_history{'GC'}{$GC_id}{'gc_affected_spores'} = $gc_affected_spores;
		    $recomb_history{'GC'}{$GC_id}{'recomb_desc'} = "GC,Type1_GC,$gc_affected_spores";
		    $recomb_history{'GC'}{$GC_id}{'chr'} = $chr;
		    $recomb_history{'GC'}{$GC_id}{'start'} = $NCO_interval_start;
		    $recomb_history{'GC'}{$GC_id}{'end'} = $NCO_interval_end;
		    $recomb_history{'GC'}{$GC_id}{'subtype'} = "Type1_GC";
		    last;
		}
	    }
	}
    }
}



# generate final linkage blocks
my $final_linkage_blocks = "$output_prefix.$tetrad_id.linkage_blocks.txt";
my $final_linkage_blocks_fh = write_file($final_linkage_blocks);
print $final_linkage_blocks_fh "tetrad_id\tblock_id\tchr\tpos_start-pos_end\tgenotype_pattern\tsegregation_pattern\n";

my %final_linkage_blocks_by_id = ();
my $final_block_id = 0;
my %tetrad_genomes = ();
foreach my $spore (@spores) {
    foreach my $chr (@coordinate_genome) {
	$tetrad_genomes{$spore}{$chr} = "";
    }
}

foreach my $chr (sort keys %recomb_breakpoints) {
    my $offset = 1;
    foreach my $bp (sort {$a <=> $b} keys %{$recomb_breakpoints{$chr}}) {
	if ($bp == 1) {
	    next;
	} else {
	    # for simplicity, here we will use closed intervals to represent each linkage blocks (i.e. not including the exact breakpoints) 
	    my $final_block_start = $offset + 1; 
	    my $final_block_end = $bp - 1;
	    my $final_block_size = $final_block_end - $final_block_start + 1;
	    $final_block_id++;
	    
	    $final_linkage_blocks_by_id{$final_block_id}{'chr'} = $chr;
	    $final_linkage_blocks_by_id{$final_block_id}{'interval_start'} = $offset;
	    $final_linkage_blocks_by_id{$final_block_id}{'interval_end'} = $bp;
	    print $final_linkage_blocks_fh "$tetrad_id\t$final_block_id\t$chr\t$offset-$bp";
	    my $parent1_genome_seq = substr $parent1_genome{$chr}, $final_block_start - 1, $final_block_size;
	    my $parent2_genome_seq = substr $parent2_genome{$chr}, $final_block_start - 1, $final_block_size;
	    my @genotype_pattern = ();
	    my $genotype_pattern = "";
	    my %genotype_pattern = ();
	    my %segregation_pattern = ();
	    my @segregation_pattern = (0, 0, 0);
	    my $segregation_pattern = "";
	    # print "genome1_seq=$parent1_genome_seq, genome2_seq=$parent2_genome_seq\n";
	    if ($offset == 1) {
		print "breakpoint: chr=$chr, offset=$offset\n";
		%genotype_pattern = %initial_genotype_pattern;
		@genotype_pattern = @initial_genotype_pattern;
		$genotype_pattern = join ":", @genotype_pattern;
		@segregation_pattern = (2, 2, 0);
		$segregation_pattern = join ":", @segregation_pattern;
		print $final_linkage_blocks_fh "\t$genotype_pattern\t$segregation_pattern\n";
		$final_linkage_blocks_by_id{$final_block_id}{'genotype_pattern'} = $genotype_pattern;
		$final_linkage_blocks_by_id{$final_block_id}{'segregation_pattern'} = $segregation_pattern;
	    } elsif ($recomb_breakpoints{$chr}{$offset}{'note'} =~ /^(CO_interval_start_and_end|CO_interval_start|CO_interval_end|NCO_interval_start)/) {
		print "breakpoint: chr=$chr, offset=$offset\n";
		print "breakpoint note: $recomb_breakpoints{$chr}{$offset}{'note'}\n";
		$final_linkage_blocks_by_id{$final_block_id}{'genotype_pattern'} = $preliminary_linkage_blocks_by_coord{$chr}{$offset}{'genotype_pattern'};
		@genotype_pattern = split ":", $final_linkage_blocks_by_id{$final_block_id}{'genotype_pattern'};
		$genotype_pattern = $final_linkage_blocks_by_id{$final_block_id}{'genotype_pattern'};
		for (my $i = 0; $i < 4; $i++) {
		    my $spore = $spores[$i];
		    my $genotype = $genotype_pattern[$i];
		    $genotype_pattern{$spore} = $genotype;
		}
		foreach my $genotype (@genotype_pattern) {
		    if ($genotype eq $parent1_tag) {
			$segregation_pattern[0]++;
		    } elsif ($genotype eq $parent2_tag) {
			$segregation_pattern[1]++;
		    } elsif ($genotype eq "NA") {
			$segregation_pattern[2]++;
		    } else {
			die "unexpected spore genotype: $genotype for the linkage block: $chr:$final_block_start-$final_block_end\n";
		    }
		}
		$segregation_pattern = join ":", @segregation_pattern;
		print $final_linkage_blocks_fh "\t$genotype_pattern\t$segregation_pattern\n";
		$final_linkage_blocks_by_id{$final_block_id}{'genotype_pattern'} = $genotype_pattern;
		$final_linkage_blocks_by_id{$final_block_id}{'segregation_pattern'} = $segregation_pattern;
	    } elsif ($recomb_breakpoints{$chr}{$offset}{'note'} =~ /NCO_interval_end/) {
		# NCO end
		print "breakpoint: chr=$chr, offset=$offset\n";
		print "breakpoint note: $recomb_breakpoints{$chr}{$offset}{'note'}\n";
		
		$genotype_pattern = $NCO_interval_end2CO_master_interval_map{$chr}{$offset}{'CO_masterblock'}{'genotype_pattern'};
		@genotype_pattern = split /:/, $genotype_pattern;
		for (my $i = 0; $i < 4; $i++) {
		    my $spore = $spores[$i];
		    my $genotype = $genotype_pattern[$i];
		    $genotype_pattern{$spore} = $genotype;
		}
		
		foreach my $genotype (@genotype_pattern) {
		    if ($genotype eq $parent1_tag) {
			$segregation_pattern[0]++;
		    } elsif ($genotype eq $parent2_tag) {
			$segregation_pattern[1]++;
		    } elsif ($genotype eq "NA") {
			$segregation_pattern[2]++;
		    } else {
			die "unexpected spore genotype: $genotype for the linkage block: $chr:$final_block_start-$final_block_end\n";
		    }
		}
		$segregation_pattern = join ":", @segregation_pattern;
		print $final_linkage_blocks_fh "\t$genotype_pattern\t$segregation_pattern\n";
		$final_linkage_blocks_by_id{$final_block_id}{'genotype_pattern'} = $genotype_pattern;
		$final_linkage_blocks_by_id{$final_block_id}{'segregation_pattern'} = $segregation_pattern;
	    } else {
		print "breakpoint: chr=$chr, offset=$offset\n";
		print "breakpoint note: $recomb_breakpoints{$chr}{$offset}{'note'}\n";
		die "unexpected breakpoints!!!\n";
	    }
	    print "tetrad genome interval: $chr:[$offset - $bp)\n"; 
	    foreach my $spore (@spores) {
		print "spore=$spore, genotype_pattern=$genotype_pattern{$spore}\n";

		my $spore_block_seq = "";
		if ($genotype_pattern{$spore} eq "$parent1_tag") {
		    $spore_block_seq = substr $parent1_genome{$chr}, $offset - 1, $bp - $offset;
		} else {
		    $spore_block_seq = substr $parent2_genome{$chr}, $offset - 1, $bp - $offset;
		}
		
		if (exists $tetrad_genomes{$spore}{$chr}) {
		    $tetrad_genomes{$spore}{$chr} .= $spore_block_seq;
		} else {
		    $tetrad_genomes{$spore}{$chr} = $spore_block_seq;
		}
	    }
	    $offset = $bp;
	}
    }
}
close $final_linkage_blocks_fh;

# generate recombination event log
my $recomb_log = "$output_prefix.$tetrad_id.recombination_events.txt";
my $recomb_log_fh = write_file($recomb_log);
my @recomb_types = qw(CO GC);
print $recomb_log_fh "tetrad_id\tevent_id\tevent_type\tevent_subtype\tchr\traw_start\traw_end\tadjusted_start\tadjusted_end\tadjusted_size\taffected_spores\tnote\n";
my $event_id = 0;
my %recomb_event_by_coord = ();
foreach my $recomb_type (@recomb_types) {
    foreach my $i (sort {$a <=> $b} keys %{$recomb_history{$recomb_type}}) {
	my $chr = $recomb_history{$recomb_type}{$i}{'chr'};
	my $start = $recomb_history{$recomb_type}{$i}{'start'};
	my $end = $recomb_history{$recomb_type}{$i}{'end'};
	my $adjusted_start = $start;
	my $adjusted_end = $end;
	if ($recomb_type eq "CO") {
	    $adjusted_start = ($start + $end)/2;
	    $adjusted_end = ($start + $end)/2;
	}
	my $size = $end - $start + 1;
	my $adjusted_size = $adjusted_end - $adjusted_start + 1;
	my $note = "";
	my $recomb_desc = $recomb_history{$recomb_type}{$i}{'recomb_desc'};
	my @recomb_desc;
	if ($recomb_desc =~ /;/) {
	    @recomb_desc = split /;/, $recomb_desc;
	} else {
	    @recomb_desc = ($recomb_desc);
	}
	foreach my $desc (@recomb_desc) {
	    my ($event_type, $event_subtype, $affected_spores) = split /,/, $desc;
	    $event_id++;
	    $recomb_event_by_coord{$chr}{$start}{$event_id} = "$tetrad_id\t$event_id\t$event_type\t$event_subtype\t$chr\t$start\t$end\t$adjusted_start\t$adjusted_end\t$adjusted_size\t$affected_spores\t$note";
	}
    }
}

foreach my $chr (sort keys %recomb_event_by_coord) {
    foreach my $start (sort {$a <=> $b} keys %{$recomb_event_by_coord{$chr}}) {
	foreach my $event_id (sort {$a <=> $b} keys %{$recomb_event_by_coord{$chr}{$start}}) {
	    print $recomb_log_fh "$recomb_event_by_coord{$chr}{$start}{$event_id}\n";
	}
    }
}

close $recomb_log_fh;



my %mutations = ();
if ($mut_num > 0) {
# overlay point mutations
# generate mutation log
    my $mut_log = "$output_prefix.$tetrad_id.mutation_events.txt";
    my $mut_log_fh = write_file($mut_log);
    my %mut_history = ();
    my $mut_type = "SNP";
    print $mut_log_fh "mutation_type\tmutation_id\tchr\tstart\tend\tsize\taffected_spores\tmutation_description\n";

    foreach (my $i = 1; $i <= $mut_num; $i++) {
	my %tmp = (); # place_holder for %invalid_recomb_regions
	my $chr_sampled = sample_chr($total_coordinate_genome_size, \%coordinate_genome_space, \%tmp);
	my $mut_size = 1;
	my ($interval_start_sampled, $interval_end_sampled) = sample_genomic_interval($mut_size, $chr_sampled, \%coordinate_genome_space, \%tmp);
	print "$mut_type:$i > sampled interval with size of $mut_size: $chr_sampled:$interval_start_sampled-$interval_end_sampled\n";
	process_mutation($mut_type, $i, $chr_sampled, $interval_start_sampled, $interval_end_sampled, \%tetrad_genomes, \%mut_history);
    }
    
    foreach my $i (sort {$a <=> $b} keys %{$mut_history{$mut_type}}) {
	my $chr = $mut_history{$mut_type}{$i}{'chr'};
	my $start = $mut_history{$mut_type}{$i}{'start'};
	my $end = $mut_history{$mut_type}{$i}{'end'};
	my $affected_spores = $mut_history{$mut_type}{$i}{'affected_spores'};
	my $size = $end - $start + 1;
	my $mut_desc = $mut_history{$mut_type}{$i}{'mut_desc'};
	$mutations{$chr}{$start}{'affected_spores'} = $affected_spores;
	$mutations{$chr}{$start}{'mut_desc'} = $mut_desc;
	print $mut_log_fh "$mut_type\t$i\t$chr\t$start\t$end\t$size\t$affected_spores\t$mut_desc\n";
    }
    close $mut_log_fh;
}


# output simulated parental genomes
my $parent1_genome_out = "$parent1_tag.genome.fa";
my $parent1_genome_out_fh = write_file($parent1_genome_out);
foreach my $chr (@coordinate_genome) {
    print $parent1_genome_out_fh ">$chr\n$parent1_genome{$chr}\n";
}
close $parent1_genome_out_fh;

my $parent2_genome_out = "$parent2_tag.genome.fa";
my $parent2_genome_out_fh = write_file($parent2_genome_out);
foreach my $chr (@coordinate_genome) {
    print $parent2_genome_out_fh ">$chr\n$parent2_genome{$chr}\n";
}
close $parent2_genome_out_fh;


# output recombined (and mutated) gamete genomes
foreach my $spore (@spores) {
    my $output_spore_genome = "$output_prefix.$tetrad_id.${spore}.genome.fa";
    my $output_spore_genome_fh = write_file($output_spore_genome);
    foreach my $chr (@coordinate_genome) {
	print $output_spore_genome_fh ">$chr\n$tetrad_genomes{$spore}{$chr}\n";
    }
    close $output_spore_genome_fh;
}


# generate tetrad genotype txt files across all marker sites

my $tetrad_genotype_txt = "$output_prefix.$tetrad_id.genotype.txt";
my $tetrad_genotype_txt_fh = write_file($tetrad_genotype_txt);
foreach my $chr (@coordinate_genome) {
    if (exists $markers{$chr}) {
	my @markers_chr_pos = sort {$a <=> $b} keys %{$markers{$chr}};
	foreach my $i (@markers_chr_pos) {
	    my $parent1_allele = substr $parent1_genome{$chr}, $i - 1, 1;
	    my $parent2_allele = $markers{$chr}{$i}{'query_allele'};
	    # print "chr=$chr i=$i coordinate_genome=$coordinate_genome\n";

	    print $tetrad_genotype_txt_fh "$chr\t$i\t$coordinate_genome";
	    my @genotype_pattern = ();
	    my @segregation_pattern = qw(0 0 0);
	    foreach my $spore (@spores) {
		my $spore_allele = substr $tetrad_genomes{$spore}{$chr}, $i - 1, 1;		
		if (($spore_allele eq $parent1_allele) and ($parent1_allele ne "NA")) {
		    push @genotype_pattern, $parent1_tag;
		    $segregation_pattern[0]++;
		} elsif (($spore_allele eq $parent2_allele) and ($parent2_allele ne "NA")) {
		    push @genotype_pattern, $parent2_tag;
		    $segregation_pattern[1]++;
		} else {
		    push @genotype_pattern, "NA";
		    $segregation_pattern[2]++;
		}
	    }
	    my $genotype_pattern = join "\t", @genotype_pattern;
	    my $segregation_pattern = join ":", @segregation_pattern;
	    # my @segregation_pattern_sorted = sort {$b <=> $a} @segregation_pattern;
	    # my $segregation_pattern_sorted = join ":", @segregation_pattern_sorted;
	    
	    print $tetrad_genotype_txt_fh "\t$genotype_pattern\n";
	}
    }
}

close $tetrad_genotype_txt_fh;

my $tetrad_genotype_input_fh = read_file($tetrad_genotype_txt);
my %genotypes = parse_genotype_file($tetrad_genotype_input_fh);
close $tetrad_genotype_input_fh;

my $tetrad_genotype_for_plotting_txt = "$output_prefix.$tetrad_id.genotype.for_genotype_plotting.txt";
my $tetrad_genotype_for_plotting_txt_fh = write_file($tetrad_genotype_for_plotting_txt);
print $tetrad_genotype_for_plotting_txt_fh "chr\tmarker_id\traw_marker_start\traw_marker_end\tadjusted_marker_start\tadjusted_marker_end\tspore_genotype\tspore_id\n";

foreach my $chr (sort keys %genotypes) {
    my @pos = sort {$a <=> $b} keys %{$genotypes{$chr}};
    my $first_marker_pos = $pos[0];
    my $last_marker_pos = $pos[-1];
    foreach my $spore_id (@spores) {
	my $marker_id = 1;
	for (my $i = 0; $i < scalar @pos; $i++) {
	    my $raw_marker_start = $pos[$i];
	    my $raw_marker_end = $pos[$i];
	    my $adjusted_marker_start;
	    my $adjusted_marker_end;
	    my $marker_genotype;
	    if ($pos[$i] == $first_marker_pos) {
		$adjusted_marker_start = $pos[$i];
		$adjusted_marker_end = ($pos[$i] + $pos[$i+1])/2;
	    } elsif ($pos[$i] == $last_marker_pos) {
		$adjusted_marker_start = ($pos[$i-1] + $pos[$i])/2;
		$adjusted_marker_end = $pos[$i];
	    } else {
		$adjusted_marker_start = ($pos[$i-1] + $pos[$i])/2;
		$adjusted_marker_end = ($pos[$i] + $pos[$i+1])/2;
	    }
	    $marker_genotype = $genotypes{$chr}{$pos[$i]}{$spore_id};
	    print $tetrad_genotype_for_plotting_txt_fh "$chr\t$marker_id\t$raw_marker_start\t$raw_marker_end\t$adjusted_marker_start\t$adjusted_marker_end\t$marker_genotype\t$spore_id\n";
	    $marker_id++;
	}
    }
}
close $tetrad_genotype_for_plotting_txt_fh;






###############
# subroutines (functions)
###############

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

sub parse_markers_table_file {
    my ($fh, $markers_hashref) = @_;
    while (<$fh>) {
        chomp;
        /^#/ and next;
	/^\s*$/ and next;
        /^chr\tstart\tend/ and next;
	/^ref_chr\tref_start\tref_end/ and next;
	/^parent1_chr\tparent1_start\tparent1_end/ and next;
	/^parent2_chr\tparent2_start\tparent2_end/ and next;
        my ($ref_chr, $ref_start, $ref_end, $ref_allele, $query_allele, $query_chr, $query_start, $query_end, $match_orientation) = split /\t/, $_;
        $$markers_hashref{$ref_chr}{$ref_start}{'ref_chr'} = $ref_chr;
        $$markers_hashref{$ref_chr}{$ref_start}{'ref_start'} = $ref_start;
        $$markers_hashref{$ref_chr}{$ref_start}{'ref_end'} = $ref_end;
        $$markers_hashref{$ref_chr}{$ref_start}{'ref_allele'} = $ref_allele;
        $$markers_hashref{$ref_chr}{$ref_start}{'query_chr'} = $query_chr;
        $$markers_hashref{$ref_chr}{$ref_start}{'query_start'} = $query_start;
        $$markers_hashref{$ref_chr}{$ref_start}{'query_end'} = $query_end;
        $$markers_hashref{$ref_chr}{$ref_start}{'query_allele'} = $query_allele;
        $$markers_hashref{$ref_chr}{$ref_start}{'match_orientation'} = $match_orientation;
    }
}

sub modify_genome_based_on_markers {
    my ($genome_hashref, $parent1_genome_hashref, $parent2_genome_hashref, $markers_hashref) = @_;
    my $count = 0;
    foreach my $chr (sort keys %$genome_hashref) {
	foreach my $pos (sort {$a <=> $b} keys %{$$markers_hashref{$chr}}) {
	    substr $$parent1_genome_hashref{$chr}, $pos - 1, 1, $$markers_hashref{$chr}{$pos}{'ref_allele'};
	    substr $$parent2_genome_hashref{$chr}, $pos - 1, 1, $$markers_hashref{$chr}{$pos}{'query_allele'};
	    $count++;
	}
    }
    print "reconstruct genome1 and genome2 by replace $count markers in the input genome\n";
}


# create genome space 
sub create_genome_space {
    my $genome_hashref = shift @_;
    my %genome_space = ();
    my $offset = 0;
    foreach my $chr (sort keys %$genome_hashref) {
	my $chr_length = length $$genome_hashref{$chr};
	my $start = $offset + 1;
	my $end = $start + $chr_length - 1;
	$genome_space{$chr}{"start"} = $start;
	$genome_space{$chr}{"end"} = $end;
	$genome_space{$chr}{"size"} = $chr_length;
	print "chr=$chr, chr_length=$chr_length, genome_space_start=$start, genome_space_end=$end\n";
	$offset = $end;
    }
    return %genome_space;
}

sub parse_linked_regions_file {
    my ($fh, $invalid_recomb_regions_hashref, $genome_space_hashref) = @_;
    while (<$fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	my ($chr, $start, $end) = split /\t/, $_;
	my $start_in_genome_space = $$genome_space_hashref{$chr}{'start'} + $start - 1;
	my $end_in_genome_space = $$genome_space_hashref{$chr}{'start'} + $end - 1;
	my $key = "${chr}:${start}-${end}";
	$$invalid_recomb_regions_hashref{$key}{'chr'} = $chr;
	$$invalid_recomb_regions_hashref{$key}{'start'} = $start;
	$$invalid_recomb_regions_hashref{$key}{'end'} = $end;
	$$invalid_recomb_regions_hashref{$key}{'start_in_genome_space'} = $start_in_genome_space;
	$$invalid_recomb_regions_hashref{$key}{'end_in_genome_space'} = $end_in_genome_space;
	if (($start == 1) and ($end == ($$genome_space_hashref{$chr}{'end'} - $$genome_space_hashref{$chr}{'start'} + 1))) {
	    $$invalid_recomb_regions_hashref{$key}{'note'} = "predefined_complete_chr";
	    print "complete linked chromosome: $chr\n";
	} else {
	    $$invalid_recomb_regions_hashref{$key}{'note'} = "predefined_region";
	}
	print "invalid recombination region: $chr:$start-$end in genome space: ${start_in_genome_space}-${end_in_genome_space}\n";
    }
}

sub genome_space_translator {
    my ($raw, $genome_space_hashref) = @_;
    my ($chr_translated, $pos_translated);
    foreach my $chr (sort keys %$genome_space_hashref) {
	my $chr_start = $$genome_space_hashref{$chr}{'start'};
	my $chr_end = $$genome_space_hashref{$chr}{'end'};
	if (($raw >= $chr_start) and ($raw <= $chr_end)) {
	    $chr_translated = $chr;
	    $pos_translated = $raw - $chr_start + 1;
	    last;
	}
    }
    return ($chr_translated, $pos_translated);
}

sub sample_chr {
    my ($total_genome_size, $genome_space_hashref, $invalid_recomb_regions_hashref) = @_;
    my $chr_sampling;
  SAMPLE:
    $chr_sampling = random_uniform(1, 1, $total_genome_size);
    print "\nchr_sampling = $chr_sampling\n";
    $chr_sampling = round($chr_sampling);
    # check for invalid chromosome
    foreach my $key (sort keys %$invalid_recomb_regions_hashref) {
	if ($$invalid_recomb_regions_hashref{$key}{'note'} eq "predefined_complete_chr") {
	    my $start_in_genome_space = $$invalid_recomb_regions_hashref{$key}{'start_in_genome_space'};
	    my $end_in_genome_space = $$invalid_recomb_regions_hashref{$key}{'end_in_genome_space'};
	    if (($chr_sampling >= $start_in_genome_space) and ($chr_sampling <= $end_in_genome_space)) {
		print "overlapped with invalid chromosome: $key, re-do sampling\n";
		goto SAMPLE;
	    }
	}
    }
    my ($chr_sampled, $start_sampled) = genome_space_translator($chr_sampling, $genome_space_hashref);
    print "chr_sampled = $chr_sampled\n";
    return $chr_sampled;
}


sub sample_genomic_interval {
    my ($interval_size, $chr, $genome_space_hashref, $invalid_recomb_regions_hashref) = @_;
    my $chr_length = $$genome_space_hashref{$chr}{'size'};
    my ($interval_start_sampled, $interval_end_sampled);
  SAMPLE:
    $interval_start_sampled = random_uniform(1, 1, $chr_length - $interval_size + 1);
    $interval_start_sampled = round($interval_start_sampled);
    $interval_end_sampled = $interval_start_sampled + $interval_size - 1;
    print "initially sampled interval: ${interval_start_sampled}-${interval_end_sampled}\n";
    # check for the overlap with invalid regions, if so, re-do sampling
    foreach my $key (sort keys %$invalid_recomb_regions_hashref) {
	if ($$invalid_recomb_regions_hashref{$key}{'chr'} eq $chr) {
	    my $invalid_region_start = $$invalid_recomb_regions_hashref{$key}{'start'};
	    my $invalid_region_end = $$invalid_recomb_regions_hashref{$key}{'end'};
            if (($interval_start_sampled <= $invalid_region_start) and ($interval_end_sampled >= $invalid_region_start)) {
		print "overlapped with ${invalid_region_start}-${invalid_region_end}, re-do sampling\n";
		goto SAMPLE;
            } elsif (($interval_end_sampled >= $invalid_region_start) and ($interval_start_sampled <= $invalid_region_end)) {
		print "overlapped with ${invalid_region_start}-${invalid_region_end}, re-do sampling\n";
		goto SAMPLE;
	    }
        }
    }
    print "final sampled interval: ${interval_start_sampled}-${interval_end_sampled}\n";
    return ($interval_start_sampled, $interval_end_sampled)
}

sub random_sample_from_array {
    my ($n, $array_ref) = @_;
    my $array_size = scalar @$array_ref;
    my @random_indices = random_permuted_index($array_size);
    my @samples = @$array_ref[@random_indices];
    return @samples;
}

sub process_mutation {
    my ($type, $id, $chr, $start, $end, $tetrad_genomes_hashref, $mut_history_hashref) = @_;
    my @spores = qw(a b c d);
    # randomly sample one spore for the mutation
    my ($spore_sampled) = random_sample_from_array(1, \@spores);
    print "spore sampled: $spore_sampled\n";
    my @affected_spores = qw(0 0 0 0);
    for (my $i = 0; $i < 4; $i++) {
	if ($spores[$i] eq $spore_sampled) {
	    $affected_spores[$i] = 1;
	    last;
	}
    }
    my $original_allele = substr $$tetrad_genomes_hashref{$spore_sampled}{$chr}, $start - 1, $end - $start + 1;
    my @alleles = qw(A T G C);
    my @possible_mut_alleles = grep {$_ ne $original_allele} @alleles;
    my ($mut_allele) = random_sample_from_array(1, \@possible_mut_alleles);
    substr $$tetrad_genomes_hashref{$spore_sampled}{$chr}, $start - 1, $end - $start + 1, $mut_allele;
    print "$original_allele -> $mut_allele\n";
    # record mutation event
    $$mut_history_hashref{$type}{$id}{'chr'} = $chr;
    $$mut_history_hashref{$type}{$id}{'start'} = $start;
    $$mut_history_hashref{$type}{$id}{'end'} = $end;
    $$mut_history_hashref{$type}{$id}{'mut_desc'} = "$original_allele->$mut_allele";
    $$mut_history_hashref{$type}{$id}{'affected_spores'} = join ":", @affected_spores;
}

sub parse_genotype_file {
    my $fh = shift @_;
    my %genotypes = ();
    while (<$fh>) {
        chomp;
        /^#/ and next;
	/^\s*$/ and next;
        my ($chr, $pos, $reftag, $genotype_a, $genotype_b, $genotype_c, $genotype_d) = split /\t/, $_;
        my @genotypes = ($genotype_a, $genotype_b, $genotype_c, $genotype_d);
        # print "chr=$chr, reftag=$reftag, pos=$pos, genotypes = @genotypes\n";
        my @spores = qw(a b c d);
        foreach my $spore_id (@spores) {
            my $gt = shift @genotypes;
            $genotypes{$chr}{$pos}{$spore_id} = $gt;
        }
    }
    return %genotypes;
}
