#!/usr/bin/env perl
# prepare soapFuse file structure using samples.txt file
# print out soapfuse samples file

use strict;
use warnings;

use File::Path qw/make_path/;

use Getopt::Std;
my %opt;
getopts('h', \%opt);

my $usage = <<ENDL;
Usage: ./prepareSoapFuse.pl [samples_sets.txt]
ENDL

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}

HELP_MESSAGE if $opt{h};

sub getReadLength {
    my $fqFile = $_[0];
    open IN, "zcat $fqFile | head |" or die "Unable to open $fqFile\n";
    <IN>;
    return length(<IN>);
}


while (my $line = <>) {
    chomp $line;
    my $sample = $line;
    my $fq1 = "fastq/$sample.1.fastq.gz";
    my $fq2 = "fastq/$sample.2.fastq.gz";
    die "Cannot find fastq files ($fq1 and $fq2)" unless (-e $fq1 && -e $fq1);
    my $sampleDir = "soapfuse/$sample/$sample";
    make_path($sampleDir);
    system "ln -f $fq1 $sampleDir/${sample}_1.fastq.gz";
    system "ln -f $fq2 $sampleDir/${sample}_2.fastq.gz";
    my $readLength = &getReadLength($fq1);

    print "$sample\t$sample\t$sample\t$readLength\n";
}



