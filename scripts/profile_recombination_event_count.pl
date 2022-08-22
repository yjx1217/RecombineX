#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: profile_recombination_event_frequency.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.10.11
#  description: summarize genome-wide recombiantion event frequency by sliding-window
#  example: perl profile_recombination_event_frequency.pl -i batch_id.DBVPG6765.pooled_events.txt -p batch_id.DBVPG6765.pooled_events -r DBVPG6765.genome.raw.fa
##############################################################

my ($input, $refseq, $prefix);

GetOptions('input|i:s' => \$input,
	   'refseq|r:s' => \$refseq,
	   'prefix|p:s' => \$prefix);

my $input_fh = read_file($input);
my %events = parse_events_file($input_fh);

my $refseq_fh = read_file($refseq);
my %refseq = ();
my @refseq = ();
parse_fasta_file($refseq_fh, \%refseq, \@refseq);

my %event_by_tetrad = ();
my %event_by_chr = ();

foreach my $chr (@refseq) {
    my $chr_length = length $refseq{$chr};
    $event_by_chr{$chr}{'length'} = $chr_length;
}

my $tetrad_count = scalar keys %events;
foreach my $tetrad_id (sort keys %events) {
    foreach my $event_id (sort keys %{$events{$tetrad_id}}) {
	my $event_chr = $events{$tetrad_id}{$event_id}{'chr'};
	my $event_type = $events{$tetrad_id}{$event_id}{'type'};
	my $event_subtype = $events{$tetrad_id}{$event_id}{'subtype'};
	my $event_midpoint = $events{$tetrad_id}{$event_id}{'adjusted_pos_midpoint'};
	if (exists $event_by_chr{$event_chr}{'event_count'}{$event_type}) {
	    $event_by_chr{$event_chr}{'event_count'}{$event_type}++;
	} else {
	    $event_by_chr{$event_chr}{'event_count'}{$event_type} = 1;
	}

	if (exists $event_by_chr{$event_chr}{'event_count'}{$event_subtype}) {
	    $event_by_chr{$event_chr}{'event_count'}{$event_subtype}++;
	} else {
	    $event_by_chr{$event_chr}{'event_count'}{$event_subtype} = 1;
	}


	if (exists $event_by_tetrad{$tetrad_id}{'event_count'}{$event_type}) {
	    $event_by_tetrad{$tetrad_id}{'event_count'}{$event_type}++;
	} else {
	    $event_by_tetrad{$tetrad_id}{'event_count'}{$event_type} = 1;
	}

	if (exists $event_by_tetrad{$tetrad_id}{'event_count'}{$event_subtype}) {
	    $event_by_tetrad{$tetrad_id}{'event_count'}{$event_subtype}++;
	} else {
	    $event_by_tetrad{$tetrad_id}{'event_count'}{$event_subtype} = 1;
	}
    }
}

my @event_type_subtype_input = qw(CO GC Type1_CO Type2_CO Type3_CO Type4_CO Type5_CO Type6_CO Type7_CO Type8_CO Type9_CO Type1_GC Type2_GC Type3_GC Type4_GC Type5_GC Type6_GC Type7_GC Type8_GC Type9_GC Type10_GC Type11_GC Type12_GC Type13_GC Type14_GC Type15_GC Type16_GC Type17_GC);
my @event_type_subtype_output = qw(all_CO all_GC single_CO double_CO single_CO_associated_GC double_CO_associated_GC simple_NCO double_NCO Type1_CO Type2_CO Type2_GC Type3_GC Type4_GC Type5_GC);


foreach my $chr (@refseq) {
    foreach my $event_subtype (@event_type_subtype_input) {
	if (not exists $event_by_chr{$chr}{'event_count'}{$event_subtype}) {
	    $event_by_chr{$chr}{'event_count'}{$event_subtype} = 0;
	}
    }
    $event_by_chr{$chr}{'event_count'}{"all_CO"} = $event_by_chr{$chr}{'event_count'}{"CO"};
    $event_by_chr{$chr}{'event_count'}{"all_GC"} = $event_by_chr{$chr}{'event_count'}{"GC"};

    $event_by_chr{$chr}{'event_count'}{"single_CO"} = $event_by_chr{$chr}{'event_count'}{"Type1_CO"} + $event_by_chr{$chr}{'event_count'}{"Type2_CO"} + $event_by_chr{$chr}{'event_count'}{"Type3_CO"} + $event_by_chr{$chr}{'event_count'}{"Type4_CO"} + $event_by_chr{$chr}{'event_count'}{"Type5_CO"} + $event_by_chr{$chr}{'event_count'}{"Type6_CO"};

    $event_by_chr{$chr}{'event_count'}{"double_CO"} = $event_by_chr{$chr}{'event_count'}{"Type7_CO"} + $event_by_chr{$chr}{'event_count'}{"Type8_CO"} + $event_by_chr{$chr}{'event_count'}{"Type9_CO"};

    $event_by_chr{$chr}{'event_count'}{"single_CO_associated_GC"} = $event_by_chr{$chr}{'event_count'}{"Type2_GC"} + $event_by_chr{$chr}{'event_count'}{"Type7_GC"} + $event_by_chr{$chr}{'event_count'}{"Type8_GC"} + $event_by_chr{$chr}{'event_count'}{"Type9_GC"} + $event_by_chr{$chr}{'event_count'}{"Type10_GC"} + $event_by_chr{$chr}{'event_count'}{"Type11_GC"} + $event_by_chr{$chr}{'event_count'}{"Type13_GC"} + $event_by_chr{$chr}{'event_count'}{"Type14_GC"} + $event_by_chr{$chr}{'event_count'}{"Type15_GC"};

    $event_by_chr{$chr}{'event_count'}{"double_CO_associated_GC"} = $event_by_chr{$chr}{'event_count'}{"Type12_GC"} + $event_by_chr{$chr}{'event_count'}{"Type16_GC"};

    $event_by_chr{$chr}{'event_count'}{"simple_NCO"} = $event_by_chr{$chr}{'event_count'}{"Type1_GC"};
    $event_by_chr{$chr}{'event_count'}{"double_NCO"} = $event_by_chr{$chr}{'event_count'}{"Type6_GC"} + $event_by_chr{$chr}{'event_count'}{"Type17_GC"};
}


foreach my $tetrad_id (sort keys %events) {
    foreach my $event_subtype (@event_type_subtype_input) {
	if (not exists $event_by_tetrad{$tetrad_id}{'event_count'}{$event_subtype}) {
	    $event_by_tetrad{$tetrad_id}{'event_count'}{$event_subtype} = 0;
	}
    }
    $event_by_tetrad{$tetrad_id}{'event_count'}{"all_CO"} = $event_by_tetrad{$tetrad_id}{'event_count'}{"CO"};
    $event_by_tetrad{$tetrad_id}{'event_count'}{"all_GC"} = $event_by_tetrad{$tetrad_id}{'event_count'}{"GC"};

    $event_by_tetrad{$tetrad_id}{'event_count'}{"single_CO"} = $event_by_tetrad{$tetrad_id}{'event_count'}{"Type1_CO"} + $event_by_tetrad{$tetrad_id}{'event_count'}{"Type2_CO"} + $event_by_tetrad{$tetrad_id}{'event_count'}{"Type3_CO"} + $event_by_tetrad{$tetrad_id}{'event_count'}{"Type4_CO"} + $event_by_tetrad{$tetrad_id}{'event_count'}{"Type5_CO"} + $event_by_tetrad{$tetrad_id}{'event_count'}{"Type6_CO"};

    $event_by_tetrad{$tetrad_id}{'event_count'}{"double_CO"} = $event_by_tetrad{$tetrad_id}{'event_count'}{"Type7_CO"} + $event_by_tetrad{$tetrad_id}{'event_count'}{"Type8_CO"} + $event_by_tetrad{$tetrad_id}{'event_count'}{"Type9_CO"};

    $event_by_tetrad{$tetrad_id}{'event_count'}{"single_CO_associated_GC"} = $event_by_tetrad{$tetrad_id}{'event_count'}{"Type2_GC"} + $event_by_tetrad{$tetrad_id}{'event_count'}{"Type7_GC"} + $event_by_tetrad{$tetrad_id}{'event_count'}{"Type8_GC"} + $event_by_tetrad{$tetrad_id}{'event_count'}{"Type9_GC"} + $event_by_tetrad{$tetrad_id}{'event_count'}{"Type10_GC"} + $event_by_tetrad{$tetrad_id}{'event_count'}{"Type11_GC"} + $event_by_tetrad{$tetrad_id}{'event_count'}{"Type13_GC"} + $event_by_tetrad{$tetrad_id}{'event_count'}{"Type14_GC"} + $event_by_tetrad{$tetrad_id}{'event_count'}{"Type15_GC"};

    $event_by_tetrad{$tetrad_id}{'event_count'}{"double_CO_associated_GC"} = $event_by_tetrad{$tetrad_id}{'event_count'}{"Type12_GC"} + $event_by_tetrad{$tetrad_id}{'event_count'}{"Type16_GC"};

    $event_by_tetrad{$tetrad_id}{'event_count'}{"simple_NCO"} = $event_by_tetrad{$tetrad_id}{'event_count'}{"Type1_GC"};
    $event_by_tetrad{$tetrad_id}{'event_count'}{"double_NCO"} = $event_by_tetrad{$tetrad_id}{'event_count'}{"Type6_GC"} + $event_by_tetrad{$tetrad_id}{'event_count'}{"Type17_GC"};
}


my $event_by_chr_output = "$prefix.recombination_event_summary.event_by_chr.txt";
my $event_by_chr_output_fh = write_file($event_by_chr_output);
foreach my $chr (@refseq) {
    my $chr_size = length $refseq{$chr};
    if (exists $event_by_chr{$chr}) {
	foreach my $t (@event_type_subtype_output) {
	    if (not exists $event_by_chr{$chr}{'event_count'}{$t}) {
		$event_by_chr{$chr}{'event_count'}{$t} = 0;
	    }
	    $event_by_chr{$chr}{'event_count_per_tetrad'}{$t} = $event_by_chr{$chr}{'event_count'}{$t}/$tetrad_count;
	    print $event_by_chr_output_fh "$chr\t$chr_size\t$t\t$event_by_chr{$chr}{'event_count'}{$t}\t$event_by_chr{$chr}{'event_count_per_tetrad'}{$t}\t$tetrad_count\n";
	}
    } else {
	foreach my $t (@event_type_subtype_output) {
	    print $event_by_chr_output_fh "$chr\t$chr_size\t$t\t0\t0\t$tetrad_count\n";
	}
    }
}


my $event_by_tetrad_output = "$prefix.recombination_event_summary.event_by_tetrad.txt";
my $event_by_tetrad_output_fh = write_file($event_by_tetrad_output);
foreach my $tetrad_id (sort keys %event_by_tetrad) {
    foreach my $t (@event_type_subtype_output) {
	if (not exists $event_by_tetrad{$tetrad_id}{'event_count'}{$t}) {
	    $event_by_tetrad{$tetrad_id}{'event_count'}{$t} = 0;
	}
	print $event_by_tetrad_output_fh "$tetrad_id\t$t\t$event_by_tetrad{$tetrad_id}{'event_count'}{$t}\n";
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
    if ($file =~ /.gz$/) {
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

sub parse_events_file {
    my $fh = shift @_;
    my %events = ();
    my $index = 0;
    while (<$fh>) {
	chomp;
	/^\s*$/ and next;
	/^#/ and next;
	/^tetrad_id/ and next;
	my @line = split /\t/, $_;
	my $tetrad_id = $line[0];
	my $event_id = $line[1];
	$events{$tetrad_id}{$event_id}{'type'} = $line[2];
	$events{$tetrad_id}{$event_id}{'subtype'} = $line[3];
	$events{$tetrad_id}{$event_id}{'chr'} = $line[4];
	$events{$tetrad_id}{$event_id}{'adjusted_pos_start'} = $line[7];
	$events{$tetrad_id}{$event_id}{'adjusted_pos_end'} = $line[8];
	$events{$tetrad_id}{$event_id}{'adjusted_pos_midpoint'} = ($line[7] + $line[8])/2;
    }
    return %events;
}



sub cal_mean {
    my @data = @_;
    my $n = scalar @data;
    my $sum = 0;
    foreach my $d (@data){
        $sum += $d;
    }
    my $mean = $sum/$n;
    return $mean;
}

sub cal_stdev {
    my @data = @_;
    my $n = scalar @data;
    my $stdev;
    if ($n > 1) {
	my $mean = cal_mean(@data);
	my $var_sum = 0;
	foreach my $d (@data){
	    $var_sum += ($d - $mean) ** 2; 
	}
	$stdev = sqrt($var_sum/($n-1));
    } else {
	$stdev = 0;
    }
    return $stdev;
}
