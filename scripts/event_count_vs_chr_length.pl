#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: event_counts_vs_chr_length.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.10.11
#  description: summarize recombination event counts versus chromosome length
#  example: perl event_count_vs_chr_length.pl -i event_summary.txt -r genome.fa -p prefix
##############################################################

my ($input, $refseq, $prefix);
$input = "Batch_S288C-SK1.S288C.q30.pooled_events.txt";
$refseq = "S288C.genome.raw.relabel.fa";
$prefix = "Batch_S288C-SK1.S288C.q30.pooled_events";

GetOptions('input|i:s' => \$input,
	   'refseq|r:s' => \$refseq,
	   'prefix|p:s' => \$prefix);

my $input_fh = read_file($input);
my %events = parse_events_file($input_fh);

my $refseq_fh = read_file($refseq);
my %refseq = ();
my @refseq = ();
parse_fasta_file($refseq_fh, \%refseq, \@refseq);

my $output = "$prefix.event_count_vs_chr_length.summary.txt";
my $output_fh = write_file($output);
print $output_fh "chr\tchr_length\tevent_type\tsample_size\tevent_count_mean\tevent_count_stderr\n";

my %chr_length = ();
foreach my $chr (@refseq) {
    $chr_length{$chr} = length $refseq{$chr};
}

my %events_by_chr = ();
foreach my $tetrad_id (sort keys %events) {
    foreach my $event_id (sort keys %{$events{$tetrad_id}}) {
	my $event_chr = $events{$tetrad_id}{$event_id}{'chr'};
	my $event_type = $events{$tetrad_id}{$event_id}{'type'};
	my $event_subtype = $events{$tetrad_id}{$event_id}{'subtype'};
	my $event_midpoint = $events{$tetrad_id}{$event_id}{'adjusted_pos_midpoint'};
	if (not exists $events_by_chr{$event_chr}{$event_type}{$tetrad_id}) {
	    $events_by_chr{$event_chr}{$event_type}{$tetrad_id} = 1;
	} else {
	    $events_by_chr{$event_chr}{$event_type}{$tetrad_id}++;
	}
	if ($event_subtype =~ /(Type1_GC|Type6_GC|Type17_GC)/) {
	    if (not exists $events_by_chr{$event_chr}{'NCO'}{$tetrad_id}) {
		$events_by_chr{$event_chr}{'NCO'}{$tetrad_id} = 1;
	    } else {
		$events_by_chr{$event_chr}{'NCO'}{$tetrad_id}++;
	    }
	}
    }
}


foreach my $chr (@refseq) {
    foreach my $event_type (sort keys %{$events_by_chr{$chr}}) {
	# print "chr=$chr, event_type=$event_type\n";
	my @type_count_chr = ();
	foreach my $tetrad_id (sort keys %events) {
	    if (not exists $events_by_chr{$chr}{$event_type}{$tetrad_id}) {
		$events_by_chr{$chr}{$event_type}{$tetrad_id} = 0;
	    }
	    push @type_count_chr, $events_by_chr{$chr}{$event_type}{$tetrad_id};
	}
	# print "check: @type_count_chr\n";
	my $sample_size = scalar @type_count_chr;
	my $type_count_chr_mean = cal_mean(@type_count_chr);
	my $type_count_chr_stdev = cal_stdev(@type_count_chr);
	my $type_count_chr_sterr = $type_count_chr_stdev/sqrt($sample_size);
	print $output_fh "$chr\t$chr_length{$chr}\t$event_type\t$sample_size\t$type_count_chr_mean\t$type_count_chr_sterr\n";
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
