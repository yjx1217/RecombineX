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

my ($input, $refseq, $window_size, $step_size, $prefix);
$window_size = 10000;
$step_size = 2000;

GetOptions('input|i:s' => \$input,
	   'refseq|r:s' => \$refseq,
	   'window|w:s' => \$window_size,
	   'step|s:s' => \$step_size,
	   'prefix|p:s' => \$prefix);

my $input_fh = read_file($input);
my %events = parse_events_file($input_fh);

my $refseq_fh = read_file($refseq);
my %refseq = ();
my @refseq = ();
parse_fasta_file($refseq_fh, \%refseq, \@refseq);

my %windows = ();
foreach my $chr (@refseq) {
    my $chr_length = length $refseq{$chr};
    my $num_of_window =  int (($chr_length - $window_size)/$step_size + 1);
    for (my $i = 0; $i < $num_of_window; $i++) {
	my $window_start = $step_size * $i + 1;
	my $window_end = $window_start + $window_size - 1;
	my $window_midpoint = ($window_start + $window_end)/2;
	$windows{$chr}{$i}{'start'} = $window_start;
	$windows{$chr}{$i}{'end'} = $window_end;
	$windows{$chr}{$i}{'midpoint'} = $window_midpoint;
	$windows{$chr}{$i}{'CO_count'} = 0;
	$windows{$chr}{$i}{'GC_count'} = 0;
	$windows{$chr}{$i}{'NCO_count'} = 0;
    }
}


foreach my $id (sort keys %events) {
    my $event_chr = $events{$id}{'chr'};
    my $event_type = $events{$id}{'type'};
    my $event_subtype = $events{$id}{'subtype'};
    my $event_midpoint = $events{$id}{'adjusted_pos_midpoint'};
    foreach my $i (sort {$a <=> $b} keys %{$windows{$event_chr}}) {
	if (($windows{$event_chr}{$i}{'start'} <= $event_midpoint) and ($windows{$event_chr}{$i}{'end'} >= $event_midpoint)) {
	    if ($event_type eq "CO") {
		$windows{$event_chr}{$i}{'CO_count'}++;
		last;
	    } elsif ($event_type eq "GC") {
		$windows{$event_chr}{$i}{'GC_count'}++;
		if ($event_subtype =~ /(Type1_GC|Type6_GC|Type17_GC)/) {
		    $windows{$event_chr}{$i}{'NCO_count'}++;
		}
		last;
	    }
	}
    }
}

my $output = "$prefix.recombination_frequency.W${window_size}S${step_size}.txt";
my $output_fh = write_file($output);
foreach my $chr (@refseq) {
    foreach my $i (sort {$a <=> $b} keys %{$windows{$chr}}) {
	print $output_fh "$chr\t$windows{$chr}{$i}{'midpoint'}\t$windows{$chr}{$i}{'CO_count'}\tCO\n";
	print $output_fh "$chr\t$windows{$chr}{$i}{'midpoint'}\t$windows{$chr}{$i}{'GC_count'}\tGC\n";
	print $output_fh "$chr\t$windows{$chr}{$i}{'midpoint'}\t$windows{$chr}{$i}{'NCO_count'}\tNCO\n";
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
	$index++;
	my $id = $index;
	$events{$id}{'type'} = $line[2];
	$events{$id}{'subtype'} = $line[3];
	$events{$id}{'chr'} = $line[4];
	$events{$id}{'adjusted_pos_start'} = $line[7];
	$events{$id}{'adjusted_pos_end'} = $line[8];
	$events{$id}{'adjusted_pos_midpoint'} = ($line[7] + $line[8])/2;
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
