#!/usr/bin/env perl
# filters paired vcfs using specified normals to determine somatic calls
# also filters LOH positions

use strict;
use warnings;

use List::Util 'sum';

use Getopt::Std;
my %opt;
getopts('n:hf:', \%opt);

# default normal var threshold
$opt{f} = 0.03 unless $opt{f};

my $usage = <<ENDL;
Usage: somaticFilterVCF.pl -n [normal] -f [normal var threshold] [vcf.file]
ENDL

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}

HELP_MESSAGE if $opt{h};
HELP_MESSAGE unless $opt{n};

my @header;
my @samples;
open IN, $ARGV[0] or die("Unable to open " . $ARGV[0] . "\n");
while (my $line = <IN>) {
    chomp $line;
    print "$line\n" and next if ($line =~ /^##/);
    if ($line =~ /^#CHROM/) {
        print "$line\n";
        @header = split /\t/, $line;
        @samples = @header[9..$#header];
        next;
    }

    my @arr = split /\t/, $line;
    my %F = map { $_ => shift @arr } @header;
    my @format = split /:/, $F{FORMAT};
    my $gtMap = {};
    for my $sample (@samples) {
        my @a = split /:/, $F{$sample};
        $gtMap->{$sample} = { map { $_ => shift(@a) } @format } ;
    }
    my $normalAF = 0;
    if (defined $gtMap->{$opt{n}}->{AD}) {
        my @AD = split /,/, $gtMap->{$opt{n}}->{AD};
        if (@AD == 2 && ($AD[0] + $AD[1]) > 0) {
            $normalAF = $AD[1] / ($AD[0] + $AD[1]);
            #    print "normal af: " . $normalAF . "\n";
        }
    }
    # check for variant genotype that isnt homozygous re
    my $nvarGT = 0; # among non-normals
    for my $sample (@samples) {
        $nvarGT++ if ($sample ne $opt{n} && ($gtMap->{$sample}->{GT} ne "0" || $gtMap->{$sample}->{GT} ne "0/0"));
    }


    # check for homozygous ref among non-normal samples
    my $nrefGT = 0;
    for my $sample (@samples) {
        $nrefGT++ if ($sample ne $opt{n} && ($gtMap->{$sample}->{GT} ne "1/1"));
    }
    my $filter = 0;
    if ($gtMap->{$opt{n}}->{GT} eq "./." || $gtMap->{$opt{n}}->{GT} ne "0/0") {
        # filter unknown normal genotypes or non-homo ref normals
        $filter++;
    } elsif ($nvarGT && $normalAF > $opt{f}) {
        # filter above threshold heterozygous
        $filter++;
    }
    #} elsif (!$nrefGT && $normalAF < 1 - $opt{f}) {
        # loss of heterozygosity, i.e. no ref genotypes among non-normal samples and homozygous ref normal with few variant reads
        #$filter++;
        #}
    if ($filter) {
        $F{"FILTER"} =~ s/:normalAD//;
        $F{"FILTER"} =~ s/normalAD//;
        #print "filt " . $F{"FILTER"} . "\n";
        if ($F{"FILTER"} eq "PASS" || $F{"FILTER"} eq "") {
            $F{"FILTER"} = "normalAD";
        } else {
            $F{"FILTER"} .= ":normalAD";
        }
    }
    print(join("\t", map {$F{$_}} @header) . "\n");
}
