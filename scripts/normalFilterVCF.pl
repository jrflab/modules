#!/usr/bin/env perl
# filter tumor based on normal
# usage: normalFilterVCF.pl [tumor.vcf] [normal.vcf]

use strict;

if (@ARGV != 2) {
    print "Usage: normalFilterVCF.pl [tumor.vcf] [normal.vcf]\n" and exit(1);
}

my $tumorVCF = $ARGV[0];
my $normalVCF = $ARGV[1];

my $varPosn = {};
open IN, $normalVCF or die("Unable to open " . $normalVCF . "\n");
while (<IN>) {
	next if /^#/;
	my @F = split /\t/;
	my $chr = $F[0];
	my $posn = $F[1];
    my $alt = $F[3];
	$varPosn->{$chr} = {} unless exists $varPosn->{$chr};
	$varPosn->{$chr}{$posn} = {} unless exists $varPosn->{$chr}{$posn};
	$varPosn->{$chr}{$posn}{$alt} = 1;
}
close IN;

open IN, $tumorVCF or die("Unable to open " . $tumorVCF . "\n");
while (<IN>) {
	print and next if /^#/;
	my @F = split /\t/;
	my $chr = $F[0];
	my $posn = $F[1];
    my $alt = $F[3];
	print unless (exists $varPosn->{$chr}{$posn}{$alt});
}
close IN;
