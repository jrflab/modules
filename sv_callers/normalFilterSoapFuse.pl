#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Std;
my %opt;
getopts('vhw:', \%opt);

my $usage = <<ENDL;
./normalFilterSoapFuse.pl [-v] -w [window size] [normal soapfuse results] [tumor soapfuse results]
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
    if ($line =~ /^Sample/) {
        @header = split /\t/, $line;
        next;
    }

    my @arr = split /\t/, $line;
    my %F = map { $_ => shift @arr } @header;

    $breakpoints{$F{"up_chr"} . ":" . $F{"up_Genome_pos"}} = $F{"up_strand"};
    $breakpoints{$F{"dw_chr"} . ":" . $F{"dw_Genome_pos"}} = $F{"dw_strand"};
}


open IN, $ARGV[1];
while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^Sample/) {
        @header = split /\t/, $line;
        print "$line\n" and next;
    }

    my @arr = split /\t/, $line;
    my %F = map { $_ => shift @arr } @header;

    my $filter = 0;
    while (my ($breakpoint, $strand) = each %breakpoints) {
        my ($chr, $pos) = split /:/, $breakpoint;
        if (($F{"up_chr"} eq $chr &&
            $F{"up_strand"} eq $strand &&
            abs($F{"up_Genome_pos"} - $pos) <= $windowSize) ||
                ($F{"dw_chr"} eq $chr &&
                $F{"dw_strand"} eq $strand &&
                abs($F{"dw_Genome_pos"} - $pos) <= $windowSize)) {
            $filter++;
            last;
        }
    }
    print STDERR "Filtered " . $F{"up_chr"} . ":" . $F{"up_Genome_pos"} . "|" . $F{"dw_chr"} . ":" . $F{"dw_Genome_pos"} . "\n" if $filter && $opt{v};
    print "$line\n" unless $filter;
}
