#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: batch_recombinex_genotype2rqtl_genotype_by_parent_genomes.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2021.02.09
#  description: converting the genotping results from RecombineX to the csv genotype input for Rqtl.
#  example: perl batch_recombinex_genotype2rqtl_genotype_by_parent_genomes.pl -b Batch_S288C-SK1 -q 30 -g . -s Master_Sample_Table.S288C-SK1.txt
##############################################################

my $batch_id = "Batch_S288C-SK1";
my $sample_table = "Master_Sample_Table.$batch_id.txt";
my $genotype_dir = "";
my $qual_diff = 30;

GetOptions('sample_table|s:s' => \$sample_table,
	   'genotype_dir|g:s' => \$genotype_dir, # input directory for genotyping txt file
           'qual_diff|q:i' => \$qual_diff, # net quality difference used in tetrad genotyping, default = 30
           'batch_id|b:s' => \$batch_id); # batch id used for genotyping 
           

my $sample_table_fh = read_file($sample_table);
my %tetrads = parse_sample_table_by_tetrad($sample_table_fh);

foreach my $tetrad_id (sort keys %tetrads) {
    my $cross_pair = $tetrads{$tetrad_id}{'cross_pair'};
    my ($parent1, $parent2) = split "-", $cross_pair;
    print "tetrad: $tetrad_id, cross_pair: $cross_pair, parent1: $parent1, parent2: $parent2\n";
    my @ref = ($parent1, $parent2);
    #my @ref = ("ref");
    foreach my $ref (@ref) {
        my @mode = ("raw", "inferred");
        foreach my $mode (@mode) {
            my $genotype_input = "$genotype_dir/$batch_id/$cross_pair.$tetrad_id.$ref.q${qual_diff}.genotype.lite.$mode.txt.gz";
            my $genotype_input_fh = read_file($genotype_input);
	    my $output = "$genotype_dir/$batch_id/$cross_pair.$tetrad_id.$ref.q${qual_diff}.genotype.lite.$mode.for_Rqtl.csv.gz";
	    my $output_fh = write_file($output);
	    my @markers = ();
	    my %genotypes = (
		a => "$tetrad_id.a",
		b => "$tetrad_id.b",
		c => "$tetrad_id.c",
		d => "$tetrad_id.d",
		);
	    while (<$genotype_input_fh>) {
		chomp;
		/^#/ and next;
		/^\s*$/ and next;
		my ($chr, $start, $tag, $GT_a, $GT_b, $GT_c, $GT_d) = split /\t+/, $_;
		my $marker = "c${chr}m${start}";
		push @markers, $marker;
		if ($GT_a eq "NA") {
		    $GT_a = "-";
		} elsif ($GT_a eq $parent1) {
		    $GT_a = 1;
		} elsif ($GT_a eq $parent2) {
                    $GT_a = 2;
		} else {
		    die "unexpected genotype for gamete a: $GT_a\n";
		}

		if ($GT_b eq "NA") {
		    $GT_b = "-";
		} elsif($GT_b eq $parent1) {
                    $GT_b = 1;
		} elsif($GT_b eq $parent2) {
                    $GT_b = 2;
		} else {
                    die"unexpected genotype for gamete b: $GT_b\n";
		}

		if ($GT_c eq "NA") {
		    $GT_c = "-";
		} elsif($GT_c eq $parent1) {
                    $GT_c = 1;
                } elsif($GT_c eq $parent2) {
                    $GT_c = 2;
                } else {
                    die"unexpected genotype for gamete c: $GT_c\n";
                }

		if ($GT_d eq "NA") {
		    $GT_d = "-";
		} elsif($GT_d eq $parent1) {
                    $GT_d = 1;
                } elsif($GT_d eq $parent2) {
                    $GT_d = 2;
                } else {
                    die"unexpected genotype for gamete d: $GT_d\n";
                }

		$genotypes{'a'} .= ",$GT_a";
		$genotypes{'b'} .= ",$GT_b";
		$genotypes{'c'} .= ",$GT_c";
		$genotypes{'d'} .= ",$GT_d";
	    }
	    
	    print $output_fh "id";
	    foreach my $marker (@markers) {
		print $output_fh ",$marker";
	    }
	    print $output_fh "\n";
	    
	    foreach my $spore (sort keys %genotypes) {
		print $output_fh "$genotypes{$spore}\n";
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

