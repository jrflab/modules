#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Std;

my %opt;
getopts('t:', \%opt);

my $usage = <<ENDL;
Usage: extractCoordsFromDefuse.pl -t [tissue type] defuse_results
-t [tissue type]: either EPI, HEM, MES or AVG
ENDL

my $tissueType;
unless ($opt{t}) {
    $tissueType = "EPI";
} else {
    $tissueType = $opt{t};
}

sub HELP_MESSAGE {
    print STDERR $usage;
    exit(1);
}


my @header;
while (my $line = <>) {
    chomp $line;

    if ($line =~ /^cluster_id/) {
        @header = split /\t/, $line;
        next;
    }

    my @arr = split /\t/, $line;
    my %F = map { $_ => shift @arr } @header;

    my $upstream = ($F{upstream_gene} eq $F{gene_name1})? 1 : 2;
    my $downstream = ($F{downstream_gene} eq $F{gene_name1})? 1 : 2;

    my $upstreamChr = "chr" . $F{"gene_chromosome" . $upstream };
    my $downstreamChr = "chr" . $F{"gene_chromosome" . $downstream };

    # give first/last nt lost
    my $upstreamPosn = ($F{"gene_strand" . $upstream} eq "+")? $F{"genomic_break_pos" . $upstream } + 1 : $F{"genomic_break_pos" . $upstream } - 1;
    my $downstreamPosn = ($F{"gene_strand" . $downstream} eq "+")? $F{"genomic_break_pos" . $downstream } - 1 : $F{"genomic_break_pos" . $downstream } + 1;

    print join("\t", ($upstreamChr, $upstreamPosn, $downstreamChr, $downstreamPosn, $tissueType)) . "\n";
}
