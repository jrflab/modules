#!/usr/bin/env perl

#use strict;
#use warnings;

#use Getopt::Std;
#my %opt;
#getopts('h', \%opt);

#my $usage <<ENDL;

#sub HELP_MESSAGE {
#print STDERR $usage;
#exit(1);
#}

#HELP_MESSAGE if $opt{h};

while (my $line = <>) {
    chomp $line;

    if ($line =~ /^cluster_id/) {
        @header = split /\t/, $line;
        push @header, "upstream_gene";
        push @header, "downstream_gene";
        #print join ("\t", 0..$#header), "\n";
        print join("\t", @header) . "\n";
        next;
    }

    my @arr = split /\t/, $line;
    my %F = map { $_ => shift @arr } @header;

    my $gene1 = $F{"gene_strand1"};
    my $gene2 = $F{"gene_strand2"};
    my $genomic1 = $F{"genomic_strand1"};
    my $genomic2 = $F{"genomic_strand2"};

    if ($gene1 eq $genomic1 && $gene2 ne $genomic2) {
        $F{"upstream_gene"} = $F{"gene_name1"};
        $F{"downstream_gene"} = $F{"gene_name2"};
    } elsif ($gene1 ne $genomic1 && $gene2 eq $genomic2) {
        $F{"upstream_gene"} = $F{"gene_name2"};
        $F{"downstream_gene"} = $F{"gene_name1"};
    } else {
        $F{"upstream_gene"} = "";
        $F{"downstream_gene"} = "";
    }

    if ($F{'orf'} eq "N") { 
        print join("\t", @F{@header}) . "\n";
        next;
    }


    if ($gene1 eq $genomic1) {
        $p5 = $F{"gene_location1"};
    } else {
        $p3 = $F{"gene_location2"};
    }

    if ($gene2 eq $genomic2) {
        $p5 = $F{"gene_location1"};
    } else {
        $p3 = $F{"gene_location2"};
    }
    if ($p5 eq "utr3p" || $p5 eq "downstream") {
        $F{"orf"} = "Y (UTR mid-fusion)";
    } elsif ($p3 eq "utr5p" || $p3 eq "upstream") {
        $F{"orf"} = "Y (UTR mid-fusion)";
    }
    print join("\t", @F{@header}) . "\n";
}
