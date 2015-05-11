#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Std;
my %opt;
getopts('h', \%opt);

my $usage = <<ENDL;
Usage: perl tableToCloneHDBaf.pl [samples] < table.txt > output.txt
ENDL

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}

HELP_MESSAGE if $opt{h};

my @samples = @ARGV;

my $headerLine = <STDIN>;
chomp $headerLine;
$headerLine =~ s/^#//;
my @header = split/\t/, $headerLine;

while (<STDIN>) {
    my @F = split /\t/;
    my %F = map { $_ => shift @F } @header;
    my $line = "$F{CHROM}\t$F{POS}";
    for my $s (@samples) {
        my $ad = $F{$s . ".AD"};
        my @ad = ($ad ne ".")? split /,/, $ad : qw/0 0/;
        my $dp = $ad[0] + $ad[1];
        $line .= "\t$ad[0]\t$dp";
    }
    print $line . "\n";
}


