#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: batch_parental_allele_frequency_in_tetrads_profiling_by_reference_genome.pl 
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.08.16
#  description: calculate parental allele frequency in tetrads
#  example: perl batch_parental_allele_frequency_in_tetrads_profiling_by_reference_genome.pl -s Master_Sample_Table.txt -batch_id batch_id -q 30 -m raw -prefix prefix 
##############################################################

my ($sample_table, $batch_id, $qual_diff, $mode, $filter, $ignore_na, $prefix);
$qual_diff = 30;
$mode = "raw"; # "raw" or "inferred" 
$filter = "viable"; # "viable" or "all"
$ignore_na = "yes"; # "yes" or "no"

GetOptions('sample_table|s:s' => \$sample_table, # master sample list
	   'batch_id|b:s' => \$batch_id,
	   'qual_diff|q:s' => \$qual_diff, # quality cutoff used for genotyping,
	   'mode|m:s' => \$mode,
	   'filter|f:s' => \$filter,
	   'ignore_na|n:s' => \$ignore_na,
	   'prefix|p:s' => \$prefix);

my $sample_table_fh = read_file($sample_table);
my %tetrads = parse_sample_table_by_tetrad($sample_table_fh);
close $sample_table_fh;

my $genotype_dir = $batch_id;

my %marker_allele_freq = (); 
my %chr = ();
my @chr = ();

my @spore_index = qw(a b c d);
my %spore_index = ();
foreach my $s (@spore_index) {
    $spore_index{$s} = 0;
}

my $cross_pair = "NA";
my $genome1_tag = "NA";
my $genome2_tag = "NA";
foreach my $tetrad_id (sort keys %tetrads) {
    $cross_pair = $tetrads{$tetrad_id}{'cross_pair'};
    ($genome1_tag, $genome2_tag) = split "-", $cross_pair;
    my @viable_spore_index = @{$tetrads{$tetrad_id}{'spore_index'}};
    foreach my $s (@viable_spore_index) {
	$spore_index{$s} = 1;
    }
    print "tetrad: $tetrad_id, cross_pair: $cross_pair, parent1: $genome1_tag, parent2: $genome2_tag, spore_index = @viable_spore_index\n";
    my @ref = ("ref");
    foreach my $ref (@ref) {
	my @mode = ($mode);
	foreach my $mode (@mode) {
	    my $genotype_input = "$genotype_dir/$cross_pair.$tetrad_id.$ref.q${qual_diff}.genotype.lite.$mode.txt.gz";
	    my $genotype_input_fh = read_file($genotype_input);
	    while (<$genotype_input_fh>) {
		chomp;
		/^#/ and next;
		/^\s*$/ and next;
		my ($chr, $pos, $reftag, $genotype_a, $genotype_b, $genotype_c, $genotype_d) = split /\t/, $_;
		if (not exists $chr{$chr}) {
		    push @chr, $chr;
		    $chr{$chr} = 1;
		}
		my @genotypes = ($genotype_a, $genotype_b, $genotype_c, $genotype_d);
		my %genotypes = ();
		$genotypes{'a'} = $genotype_a;
		$genotypes{'b'} = $genotype_b;
		$genotypes{'c'} = $genotype_c;
		$genotypes{'d'} = $genotype_d;
		
		foreach my $s (@spore_index) {
		    my $genotype = $genotypes{$s};
		    if ($filter eq "viable") {
			if ($spore_index{$s} == 1) {
			    if (exists $marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}) {
				$marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}{'count'}++;
			    } else {
				$marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}{'count'} = 1;
			    }
			}
		    } elsif ($filter eq "all") {
			if (exists $marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}) {
			    $marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}{'count'}++;
			} else {
			    $marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}{'count'} = 1;
			}
		    }
		}
		if (not exists $marker_allele_freq{$reftag}{$chr}{$pos}{$genome1_tag}) {
		    $marker_allele_freq{$reftag}{$chr}{$pos}{$genome1_tag}{'count'} = 0;
		} elsif (not exists $marker_allele_freq{$reftag}{$chr}{$pos}{$genome2_tag}) {
                    $marker_allele_freq{$reftag}{$chr}{$pos}{$genome2_tag}{'count'} = 0;
		} elsif (not exists $marker_allele_freq{$reftag}{$chr}{$pos}{'NA'}) {
                    $marker_allele_freq{$reftag}{$chr}{$pos}{'NA'}{'count'} = 0;
		} elsif (not exists $marker_allele_freq{$reftag}{$chr}{$pos}{'heteroduplex'}) {
                    $marker_allele_freq{$reftag}{$chr}{$pos}{'heteroduplex'}{'count'} = 0;
		}
	    }
	}
    }
}

# count to frequency
my $output_by_ref = "$prefix.$cross_pair.ref.parental_allele_frequency.$mode.txt";
my $output_by_ref_fh = write_file($output_by_ref);
print $output_by_ref_fh "cross_pair\tref\tchr\tpos\tgenotype\ttotal_allele_count\tgenotype_match_allele_count\tfreq\n";

foreach my $reftag (sort keys %marker_allele_freq) {
    print "output cross_pair=$cross_pair, reftag = $reftag\n";
    foreach my $chr (@chr) {
	print "output chr = $chr\n";
	foreach my $pos (sort {$a <=> $b} keys %{$marker_allele_freq{$reftag}{$chr}}) {
	    my $total_allele_count = 0;
	    foreach my $genotype (sort keys %{$marker_allele_freq{$reftag}{$chr}{$pos}}) {
		if ($ignore_na eq "yes") { 
		    if ($genotype ne "NA") {
			$total_allele_count += $marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}{'count'};
		    }
		} else {
		    $total_allele_count += $marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}{'count'};
		}
	    }
	    foreach my $genotype (sort keys %{$marker_allele_freq{$reftag}{$chr}{$pos}}) {
		# print "reftag=$reftag, chr=$chr, pos=$pos\n";
		if ($total_allele_count > 0) {
		    my $freq = $marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}{'count'}/$total_allele_count;
		    $marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}{'freq'} = sprintf("%.2f", $freq);
		    if ($ignore_na eq "yes") {
			if ($genotype ne "NA") {
				print $output_by_ref_fh "$cross_pair\t$reftag\t$chr\t$pos\t$genotype\t$total_allele_count\t$marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}{'count'}\t$marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}{'freq'}\n";
			}
		    } else {
			print $output_by_ref_fh "$cross_pair\t$reftag\t$chr\t$pos\t$genotype\t$total_allele_count\t$marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}{'count'}\t$marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}{'freq'}\n";
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
    
