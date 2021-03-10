#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: filter_vcf_by_window.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2019.09.07
#  description: filter variants (both SNP and INDELs) in VCF files that are within a window of defined size 
#  example: perl filter_vcf_by_window.pl -i input.vcf(.gz) -o output.vcf(.gz) -w 10 
##############################################################

my ($input, $output, $window);
$window = 10;

GetOptions('input|i:s' => \$input,
	   'output|o:s' => \$output,
	   'window|w:i' => \$window);

my $input_fh = read_file($input);
my $output_fh = write_file($output);

my %buffer = ();
while (<$input_fh>) {
    chomp;
    (/^\s*$/) and next;
    if (/^#/) {
	print $output_fh "$_\n";
    } else {
	my @line = split /\t/, $_;
	my $chr = $line[0];
	my $pos = $line[1];
	my $ref_allele = $line[3];
	my $alt_allele = $line[4];
	if (not %buffer) {
	    $buffer{'chr'} = $chr;
	    $buffer{'pos'} = $pos;
	    $buffer{'ref_allele'} = $ref_allele;
	    $buffer{'alt_allele'} = $alt_allele;
	    $buffer{'line'} = $_;
	    $buffer{'flag'} = 0;
	} else {
	    if ($chr ne $buffer{'chr'}) {
		if ($buffer{'flag'} == 0) {
		    print $output_fh "$buffer{'line'}\n";
		}
		$buffer{'chr'} = $chr;
		$buffer{'pos'} = $pos;
		$buffer{'ref_allele'} = $ref_allele;
		$buffer{'alt_allele'} = $alt_allele;
		$buffer{'line'} = $_;
		$buffer{'flag'} = 0;
	    } else {
		if ((abs($pos - $buffer{'pos'}) > $window)) {
		    if ($buffer{'flag'} == 0) {
			print $output_fh "$buffer{'line'}\n";
		    }
		    $buffer{'chr'} = $chr;
		    $buffer{'pos'} = $pos;
		    $buffer{'ref_allele'} = $ref_allele;
		    $buffer{'alt_allele'} = $alt_allele;
		    $buffer{'line'} = $_;
		    $buffer{'flag'} = 0;
		} else {
		    my $current_var_type = SNP_or_INDEL($ref_allele, $alt_allele);
		    my $buffer_var_type = SNP_or_INDEL($buffer{'ref_allele'}, $buffer{'alt_allele'});
		    if (($current_var_type eq "SNP") and ($buffer_var_type eq "SNP")) {
			if (abs($pos - $buffer{'pos'}) > 1) {
			    if ($buffer{'flag'} == 0) {
				print $output_fh "$buffer{'line'}\n";
			    }
			    $buffer{'chr'} = $chr;
			    $buffer{'pos'} = $pos;
			    $buffer{'ref_allele'} = $ref_allele;
			    $buffer{'alt_allele'} = $alt_allele;
			    $buffer{'line'} = $_;
			    $buffer{'flag'} = 0;
			} else {
			    $buffer{'chr'} = $chr;
			    $buffer{'pos'} = $pos;
			    $buffer{'ref_allele'} = $ref_allele;
			    $buffer{'alt_allele'} = $alt_allele;
			    $buffer{'line'} = $_;
			    $buffer{'flag'} = 1;
			}
		    } else {
			# at least one of the current and buffer records is INDEL
			$buffer{'chr'} = $chr;
			$buffer{'pos'} = $pos;
			$buffer{'ref_allele'} = $ref_allele;
			$buffer{'alt_allele'} = $alt_allele;
			$buffer{'line'} = $_;
			$buffer{'flag'} = 1;
		    }
		}
	    } 
	}
    }
}
if ((%buffer) and ($buffer{'flag'} == 0)) {
    print $output_fh "$buffer{'line'}\n";
}


sub read_file {
    my $file = shift @_;
    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "gunzip -c $file |") or die "can't open pipe to $file";
    }
    else {
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

sub SNP_or_INDEL {
    my ($ref_allele, $alt_allele) = @_;
    my $var_type;
    if (($ref_allele =~ /^[ATGCNatgcn]$/) and ($alt_allele =~ /^[ATGCNatgcn]$/)) {
	$var_type = "SNP";
    } else {
	$var_type = "INDEL";
    }
    return $var_type;
}





