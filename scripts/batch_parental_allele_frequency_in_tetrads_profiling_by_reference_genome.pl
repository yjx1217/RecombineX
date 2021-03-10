#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: batch_parental_allele_frequency_in_tetrads_profiling_by_reference_genome.pl 
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.10.12
#  description: calculate parental allele frequency in tetrads
#  example: perl batch_parental_allele_frequency_in_tetrads_profiling_by_reference_genome.pl -s Master_Sample_Table.txt -batch_id batch_id -q 30 -m raw -output output.parental_allele_frequency.raw.txt.gz
##############################################################

my ($sample_table, $batch_id, $qual_diff, $mode, $output);
$output = "output.parental_allele_frequency.txt.gz";
$qual_diff = 30;
$mode = "raw";

GetOptions('sample_table|s:s' => \$sample_table, # master sample list
	   'batch_id|b:s' => \$batch_id,
	   'qual_diff|q:s' => \$qual_diff, # quality cutoff used for genotyping,
	   'mode|m:s' => \$mode,
	   'output|o:s' => \$output);

my $sample_table_fh = read_file($sample_table);
my %tetrads = parse_sample_table_by_tetrad($sample_table_fh);
close $sample_table_fh;

my $output_fh = write_file($output);
print $output_fh "ref\tchr\tpos\tgenotype\ttotal_allele_count\tgenotype_match_allele_count\tfreq\n";

my $genotype_dir = $batch_id;

my %marker_allele_freq = (); 
my %chr = ();
my @chr = ();

foreach my $tetrad_id (sort keys %tetrads) {
    my $cross_pair = $tetrads{$tetrad_id}{'cross_pair'};
    my ($genome1_tag, $genome2_tag) = split "-", $cross_pair;
    print "tetrad: $tetrad_id, cross_pair: $cross_pair, parent1: $genome1_tag, parent2: $genome2_tag\n";
    # my @ref = ($genome1_tag, $genome2_tag);
    my @ref = ("ref");
    foreach my $ref (@ref) {
	# my @mode = ("raw", "inferred");
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
		
		foreach my $genotype (@genotypes) {
		    if (($genotype ne $genome1_tag) and ($genotype ne $genome2_tag)) {
			$genotype = "NA";
		    }
		    if (exists $marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}) {
			$marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}{'count'}++;
		    } else {
			$marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}{'count'} = 1;
		    }
		}
		if (not exists $marker_allele_freq{$reftag}{$chr}{$pos}{$genome1_tag}) {
		    $marker_allele_freq{$reftag}{$chr}{$pos}{$genome1_tag}{'count'} = 0;
		} elsif (not exists $marker_allele_freq{$reftag}{$chr}{$pos}{$genome2_tag}) {
                    $marker_allele_freq{$reftag}{$chr}{$pos}{$genome2_tag}{'count'} = 0;
		} elsif (not exists $marker_allele_freq{$reftag}{$chr}{$pos}{'NA'}) {
                    $marker_allele_freq{$reftag}{$chr}{$pos}{'NA'}{'count'} = 0;
		}
	    }
	}
    }
}

# count to frequency

foreach my $reftag (sort keys %marker_allele_freq) {
    print "output reftag = $reftag\n";
    foreach my $chr (@chr) {
	print "output chr = $chr\n";
	foreach my $pos (sort {$a <=> $b} keys %{$marker_allele_freq{$reftag}{$chr}}) {
	    my $total_allele_count = 0;
	    foreach my $genotype (sort keys %{$marker_allele_freq{$reftag}{$chr}{$pos}}) {
		$total_allele_count += $marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}{'count'};
	    }
	    foreach my $genotype (sort keys %{$marker_allele_freq{$reftag}{$chr}{$pos}}) {
		my $freq = $marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}{'count'}/$total_allele_count;
		$marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}{'freq'} = sprintf("%.2f", $freq);
		print $output_fh "$reftag\t$chr\t$pos\t$genotype\t$total_allele_count\t$marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}{'count'}\t$marker_allele_freq{$reftag}{$chr}{$pos}{$genotype}{'freq'}\n";
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
    
