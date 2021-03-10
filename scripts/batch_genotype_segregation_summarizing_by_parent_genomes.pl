#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: batch_genotype_segregation_summarizing_by_parent_genomes.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2019.09.23
#  description: run batch summarizing the genotype segregation ratio for called tetrad genotype
#  example: perl batch_genotype_segregation_summarizing_by_parent_genomes.pl -s Master_Sample_Table.txt -genotype_dir genotype_dir -b batch_id -q qual_diff -output output.summary.txt
##############################################################

my ($sample_table, $genotype_dir, $batch_id, $qual_diff, $output);
$batch_id = "Batch_TEST";
$output = "$batch_id.segregation_summary.txt";

GetOptions('sample_table|s:s' => \$sample_table, # master sample list
	   'batch_id|b:s' => \$batch_id,
	   'genotype_dir|gt_dir:s' => \$genotype_dir,
	   'qual_diff|q:s' => \$qual_diff, # quality cutoff used for genotyping
	   'output|o:s' => \$output);

my $sample_table_fh = read_file($sample_table);
my %tetrads = parse_sample_table_by_tetrad($sample_table_fh);
close $sample_table_fh;

my $output_fh = write_file($output);
print $output_fh "tetrad_id\tcross_pair\treference\traw_2:2_pct\traw_3:1_pct\traw_1:3_pct\traw_4:0_pct\traw_0:4_pct\traw_NA_pct\tinferred_2:2_pct\tinferred_3:1_pct\tinferred_1:3_pct\tinferred_4:0_pct\tinferred_0:4_pct\tinferred_NA_pct\n";

foreach my $tetrad_id (sort keys %tetrads) {
    my $cross_pair = $tetrads{$tetrad_id}{'cross_pair'};
    my ($genome1_tag, $genome2_tag) = split "-", $cross_pair;
    print "tetrad: $tetrad_id, cross_pair: $cross_pair, parent1: $genome1_tag, parent2: $genome2_tag\n";
    my @ref = ($genome1_tag, $genome2_tag);
    foreach my $ref (@ref) {
	print $output_fh "$tetrad_id\t$cross_pair\t$ref";
	my @mode = ("raw", "inferred");
	foreach my $mode (@mode) {
	    my $genotype_input = "$genotype_dir/$batch_id/$cross_pair.$tetrad_id.$ref.q${qual_diff}.genotype.lite.$mode.txt.gz";
	    my $genotype_input_fh = read_file($genotype_input);
	    my @segregation_patterns = qw(2:2 3:1 1:3 4:0 0:4 NA); # parent1:parent2:NA
	    my %genotype_segregation = ();
	    foreach my $p (@segregation_patterns) {
		$genotype_segregation{$p}{'count'} = 0;
	    }
	    my $total_marker_count = 0;
	    while (<$genotype_input_fh>) {
		chomp;
		/^#/ and next;
		/^\s*$/ and next;
		my ($chr, $pos, $reftag, $genotype_a, $genotype_b, $genotype_c, $genotype_d) = split /\t/, $_;
		$total_marker_count++;
		my $parent1_genotype_count = 0;
		my $parent2_genotype_count = 0;
		my $NA_genotype_count = 0;
		my @genotypes = ($genotype_a, $genotype_b, $genotype_c, $genotype_d);
		foreach my $genotype (@genotypes) {
		    if ($genotype eq $genome1_tag) {
			    $parent1_genotype_count++;
		    } elsif ($genotype eq $genome2_tag) {
			    $parent2_genotype_count++;
		    } else {
			$NA_genotype_count++;
		    }
		}
		if ($NA_genotype_count > 0) {
		    $genotype_segregation{'NA'}{'count'}++;
		} elsif ("$parent1_genotype_count:$parent2_genotype_count" eq "2:2") {
		    $genotype_segregation{'2:2'}{'count'}++;
		} elsif ("$parent1_genotype_count:$parent2_genotype_count" eq "3:1") {
		    $genotype_segregation{'3:1'}{'count'}++;
		} elsif ("$parent1_genotype_count:$parent2_genotype_count" eq "1:3") {
		    $genotype_segregation{'1:3'}{'count'}++;
		} elsif ("$parent1_genotype_count:$parent2_genotype_count" eq "4:0") {
		    $genotype_segregation{'4:0'}{'count'}++;
		} elsif ("$parent1_genotype_count:$parent2_genotype_count" eq "0:4") {
		    $genotype_segregation{'0:4'}{'count'}++;
		}
	    }
	    foreach my $p (@segregation_patterns) {
		$genotype_segregation{$p}{'proportion'} = $genotype_segregation{$p}{'count'}/$total_marker_count ;
		$genotype_segregation{$p}{'pct'} = sprintf( "%.2f", 100 * $genotype_segregation{$p}{'proportion'});
		print $output_fh "\t$genotype_segregation{$p}{'pct'}";
	    }
	}
	print $output_fh "\n";
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
    
