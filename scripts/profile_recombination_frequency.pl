#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: profile_recombination_event_frequency.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.10.11
#  description: summarize genome-wide recombiantion event frequency by sliding-window
#  example: perl profile_recombination_event_frequency.pl -i batch_id.DBVPG6765.pooled_events.txt -p batch_id.DBVPG6765.pooled_events -r DBVPG6765.genome.raw.fa -n normalize no
##############################################################

my ($input, $refseq, $window_size, $step_size, $normalize, $prefix);
$window_size = 10000;
$step_size = 2000;
$normalize = "no";

GetOptions('input|i:s' => \$input,
	   'refseq|r:s' => \$refseq,
	   'window|w:s' => \$window_size,
	   'step|s:s' => \$step_size,
	   'normalize|n:s' => \$normalize,
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
	$windows{$chr}{$i}{'all_CO'}{'count'} = 0;
	$windows{$chr}{$i}{'all_GC'}{'count'} = 0;
	$windows{$chr}{$i}{'single_CO'}{'count'} = 0;
	$windows{$chr}{$i}{'double_CO'}{'count'} = 0;
	$windows{$chr}{$i}{'single_CO_associated_GC'}{'count'} = 0;
	$windows{$chr}{$i}{'double_CO_associated_GC'}{'count'} = 0;
	$windows{$chr}{$i}{'simple_NCO'}{'count'} = 0;
	$windows{$chr}{$i}{'double_NCO'}{'count'} = 0;
	$windows{$chr}{$i}{'Type1_CO'}{'count'} = 0;
	$windows{$chr}{$i}{'Type2_CO'}{'count'} = 0;
	$windows{$chr}{$i}{'Type1_GC'}{'count'} = 0;
	$windows{$chr}{$i}{'Type2_GC'}{'count'} = 0;
	$windows{$chr}{$i}{'Type3_GC'}{'count'} = 0;
	$windows{$chr}{$i}{'Type4_GC'}{'count'} = 0;
	$windows{$chr}{$i}{'Type5_GC'}{'count'} = 0;
    }
}

my $tetrad_count = scalar keys %events;
foreach my $tetrad_id (sort keys %events) {
    foreach my $event_id (sort keys %{$events{$tetrad_id}}) {
	my $event_chr = $events{$tetrad_id}{$event_id}{'chr'};
	my $event_type = $events{$tetrad_id}{$event_id}{'type'};
	my $event_subtype = $events{$tetrad_id}{$event_id}{'subtype'};
	my $event_midpoint = $events{$tetrad_id}{$event_id}{'adjusted_pos_midpoint'};
	foreach my $i (sort {$a <=> $b} keys %{$windows{$event_chr}}) {
	    if (($windows{$event_chr}{$i}{'start'} <= $event_midpoint) and ($windows{$event_chr}{$i}{'end'} >= $event_midpoint)) {
		if ($event_type eq "CO") {
		    $windows{$event_chr}{$i}{'all_CO'}{'count'}++;
		    if ($event_subtype =~ /(Type1_CO|Type2_CO|Type3_CO|Type4_CO|Type5_CO|Type6_CO)/) {
			$windows{$event_chr}{$i}{'single_CO'}{'count'}++;
		    }
		    if ($event_subtype eq "Type1_CO") {
			$windows{$event_chr}{$i}{'Type1_CO'}{'count'}++;
		    }

		    if ($event_subtype eq "Type2_CO") {
			$windows{$event_chr}{$i}{'Type2_CO'}{'count'}++;
		    }

		    if ($event_subtype =~ /(Type7_CO|Type8_CO|Type9_CO)/) {
			$windows{$event_chr}{$i}{'double_CO'}{'count'}++;
		    }
		    
		} elsif ($event_type eq "GC") {
		    $windows{$event_chr}{$i}{'all_GC'}{'count'}++;
		    if ($event_subtype =~ /(Type2_GC|Type7_GC|Type8_GC|Type9_GC|Type10_GC|Type11_GC|Type13_GC|Type14_GC|Type15_GC)/) {
			$windows{$event_chr}{$i}{'single_CO_associated_GC'}{'count'}++;
		    }

		    if ($event_subtype =~ /(Type12_GC|Type16_GC)/) {
			$windows{$event_chr}{$i}{'double_CO_associated_GC'}{'count'}++;
		    }

		    if ($event_subtype =~ /(Type6_GC|Type17_GC)/) {
			$windows{$event_chr}{$i}{'double_NCO'}{'count'}++;
		    }

		    if ($event_subtype eq "Type1_GC") {
			$windows{$event_chr}{$i}{'simple_NCO'}{'count'}++;
			$windows{$event_chr}{$i}{'Type1_GC'}{'count'}++;
		    }

		    if ($event_subtype eq "Type2_GC") {
			$windows{$event_chr}{$i}{'Type2_GC'}{'count'}++;
		    }
		    if ($event_subtype eq "Type3_GC") {
			$windows{$event_chr}{$i}{'Type3_GC'}{'count'}++;
		    }
		    if ($event_subtype eq "Type4_GC") {
			$windows{$event_chr}{$i}{'Type4_GC'}{'count'}++;
		    }
		    if ($event_subtype eq "Type5_GC") {
			$windows{$event_chr}{$i}{'Type5_GC'}{'count'}++;
		    }

		}
	    }
	}
    }
}


my $output;
if ($normalize eq "yes") {
    $output = "$prefix.recombination_frequency.W${window_size}S${step_size}.txt";
} else {
    $output = "$prefix.recombination_count.W${window_size}S${step_size}.txt";
}
my $output_fh = write_file($output);

my @event_type_subtype_output = qw(all_CO all_GC single_CO double_CO single_CO_associated_GC double_CO_associated_GC simple_NCO double_NCO Type1_CO Type2_CO Type2_GC Type3_GC Type4_GC Type5_GC);

foreach my $chr (@refseq) {
    foreach my $i (sort {$a <=> $b} keys %{$windows{$chr}}) {
	foreach my $subtype (@event_type_subtype_output) {
	    if (not exists $windows{$chr}{$i}{$subtype}{'count'}) {
		$windows{$chr}{$i}{$subtype}{'count'} = 0;
	    }
	    $windows{$chr}{$i}{$subtype}{'freq'} = $windows{$chr}{$i}{$subtype}{'count'}/$tetrad_count;
	    if ($normalize eq "yes") {
		print $output_fh "$chr\t$windows{$chr}{$i}{'midpoint'}\t$windows{$chr}{$i}{$subtype}{'freq'}\t$subtype\n";
	    } else {
		print $output_fh "$chr\t$windows{$chr}{$i}{'midpoint'}\t$windows{$chr}{$i}{$subtype}{'count'}\t$subtype\n";
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
