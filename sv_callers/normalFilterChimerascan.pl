#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Std;
my %opt;
getopts('vhw:', \%opt);

my $usage = <<ENDL;
./normalFilterChimerascan.pl [-v] -w [window size] [normal chimscan results] [tumor chimscan results]
    -w: window size [default is 50kb]
    -v: verbose
ENDL

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}

HELP_MESSAGE if $opt{h};

print STDERR "ERROR: Need normal and tumor chimerascan tables\n" and HELP_MESSAGE if @ARGV != 2;

my $windowSize = 50000;
$windowSize = $opt{w} if $opt{w};

my @header;
my %breakpoints;
open IN, $ARGV[0];
while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /chrom5p/) {
        $line =~ s/#//;
        @header = split /\t/, $line;
        next;
    }

    my @arr = split /\t/, $line;
    my %F = map { $_ => shift @arr } @header;

    $breakpoints{$F{"chrom5p"} . ":" . $F{"end5p"}} = $F{"strand5p"};
    $breakpoints{$F{"chrom3p"} . ":" . $F{"start3p"}} = $F{"strand3p"};
}

open IN, $ARGV[1];
while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /chrom5p/) {
        $line =~ s/#//;
        @header = split /\t/, $line;
        print "$line\n" and next;
    }

    my @arr = split /\t/, $line;
    my %F = map { $_ => shift @arr } @header;

    my $filter = 0;
    while (my ($breakpoint, $strand) = each %breakpoints) {
        my ($chr, $pos) = split /:/, $breakpoint;
        if (($F{"chrom5p"} eq $chr &&
            $F{"strand5p"} eq $strand &&
            abs($F{"end5p"} - $pos) <= $windowSize) ||
                ($F{"chrom3p"} eq $chr &&
                $F{"strand3p"} eq $strand &&
                abs($F{"start3p"} - $pos) <= $windowSize)) {
            $filter++;
            last;
        }
    }
    print STDERR "Filtered " . $F{"chimera_cluster_id"} . "\n" if $filter && $opt{v};
    print "$line\n" unless $filter;
}
