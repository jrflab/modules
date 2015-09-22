#!/usr/bin/env perl
# parse samples file to get sample sets (space delimited, normal last)

use strict;
use warnings;

my %tumorSamples;
my %normalSamples;
for my $s (<>) {
    chomp $s;
    my $id;
    if ($s =~ /(\d+)/) {
        $id = $1;
        if ($s =~ m/N$/) {
            print STDERR "Warning: sample $id ($s) has two normals\n" if (exists $normalSamples{$id});
            $normalSamples{$id} = $s;
        } else {
            push @{$tumorSamples{$id}}, $s;
        }
    }
}

while (my ($id, $normal) = each %normalSamples) {
    next and print STDERR "Warning: no tumor samples for $id ($normal)" unless (exists $tumorSamples{$id});
    print join(" ", @{$tumorSamples{$id}}) . " $normal\n";
}

for my $id (keys %tumorSamples) {
    unless (exists $normalSamples{$id}) {
        print STDERR "Warning: no normal sample for $id (" . join(" ", @{$tumorSamples{$id}}) . ")\n";
    }
}
