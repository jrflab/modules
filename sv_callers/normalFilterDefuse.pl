#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Std;
my %opt;
getopts('vhw:', \%opt);

my $usage = <<ENDL;
./normalFilterDefuse.pl [-v] -w [window size] [normal defuse results] [tumor defuse results]
    -w: window size [default is 50kb]
    -v: verbose
ENDL

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}

HELP_MESSAGE if $opt{h};

print STDERR "ERROR: Need normal and tumor defuse tables\n" and HELP_MESSAGE if @ARGV != 2;

my $windowSize = 50000;
$windowSize = $opt{w} if $opt{w};

my @header;
my %breakpoints;
open IN, $ARGV[0];
while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^cluster_id/) {
        @header = split /\t/, $line;
        next;
    }

    my @arr = split /\t/, $line;
    my %F = map { $_ => shift @arr } @header;

    $breakpoints{$F{"gene_chromosome1"} . ":" . $F{"genomic_break_pos1"}} = $F{"genomic_strand1"};
    $breakpoints{$F{"gene_chromosome2"} . ":" . $F{"genomic_break_pos2"}} = $F{"genomic_strand2"};
}

open IN, $ARGV[1];
while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^cluster_id/) {
        @header = split /\t/, $line;
        print "$line\n" and next;
    }

    my @arr = split /\t/, $line;
    my %F = map { $_ => shift @arr } @header;

    my $filter = 0;
    while (my ($breakpoint, $strand) = each %breakpoints) {
        my ($chr, $pos) = split /:/, $breakpoint;
        if (($F{"gene_chromosome1"} eq $chr &&
            $F{"genomic_strand1"} eq $strand &&
            abs($F{"genomic_break_pos1"} - $pos) <= $windowSize) ||
                ($F{"gene_chromosome2"} eq $chr &&
                $F{"genomic_strand2"} eq $strand &&
                abs($F{"genomic_break_pos2"} - $pos) <= $windowSize)) {
            $filter++;
            last;
        }
    }
    print STDERR "Filtered " . $F{"cluster_id"} . "\n" if $filter && $opt{v};
    print "$line\n" unless $filter;
}
