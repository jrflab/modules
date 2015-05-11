#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Std;
my %opt;
getopts('hb:d:', \%opt);

my $usage = <<ENDL;
perl countNonRefPileup.pl -d [min depth] -b [bin size] < [pileup]
ENDL

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}

HELP_MESSAGE if $opt{h};

my @bases = qw/A T C G/;
my $freqThreshold = 0.001;
my $binSize = $opt{b}? $opt{b} : 0.05;
my $minDepth = $opt{d}? $opt{d} : 10;
my $lastBin = 1 / $binSize;

my %nonRefFreq;
while (<>) {
    chomp;
    my @F = split /\t/;
    next unless $F[4];
    my $depth = length($F[4]);
    next if $depth < $minDepth;
    next if $F[4] =~ /[-+]/;
    $F[4] =~ s/[-+][ATCGNatcgn]+//g;
    for my $base (@bases) {
        my $count = () = $F[4] =~ /$base/i;
        my $freq = $count / $depth;
        $nonRefFreq{$F[2]}{$base}[int($freq / $binSize)]++ if $freq > $freqThreshold;
    }
}

print "Ref\tVar";
for my $bin (0..$lastBin) {
    print "\tBin$bin";
}
print "\n";
for my $refBase (@bases) {
    for my $varBase (@bases) {
        next if $refBase eq $varBase;
        print "$refBase\t$varBase";
        for my $bin (0..$lastBin) {
            print "\t";;
            if (defined $nonRefFreq{$refBase}{$varBase}[$bin]) {
                print $nonRefFreq{$refBase}{$varBase}[$bin];
            } else {
                print "0";
            }
        }
        print "\n";
    }
}
