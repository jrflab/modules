#!/usr/bin/env perl
# join EFF lines

use strict;
use List::MoreUtils qw(first_index indexes);

my $line = <>;
print $line;
chomp $line;
my @header = split /\t/, $line;
my @effIndexes = indexes { /^ANN\[/ } @header;

my %lines;
while (<>) {
    chomp;
    my @F = split /\t/, $_, -1;
    for my $i (0..$#F) {
        $F[$i] = "." unless $F[$i] =~ /\S/;
    }
    push @{$lines{$F[0]}{$F[1]}}, \@F;
}

foreach my $chrom (sort keys %lines) {
    foreach my $posn (sort keys %{$lines{$chrom}}) {
        my $F = pop @{$lines{$chrom}{$posn}};
        while (my $Fn = pop @{$lines{$chrom}{$posn}}) {
            for my $i (@effIndexes) {
                $F->[$i] .= "|" . $Fn->[$i];
            }
        }
        print join("\t", @{$F}) . "\n";
    }
}
    
