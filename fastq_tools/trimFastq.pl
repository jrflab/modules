#!/usr/bin/env perl
# trim fastq file

use strict;
use warnings;

use Getopt::Std;

my %opt;
getopts('hs:e:l:', \%opt);

my $usage = <<ENDL;
Usage: filterFastq.pl -l [read length] 
-h: this help message
-l [integer]: max length
-s [integer]: trim x bases from start of read
-e [integer]: trim x bases from end of read
ENDL

sub HELP_MESSAGE {
    print STDERR $usage;
    exit(1);
}

print "Missing read-trim length\n" and HELP_MESSAGE() unless ($opt{l} || $opt{s} || $opt{e});

my $i = 0;
while (<STDIN>) {
    chomp;
    if ($i % 2 == 0) {
        print;
    } else {
        my $ss = $_;
        if ($opt{s}) {
            $ss = substr($ss, $opt{s});
        }
        if ($opt{e}) {
            $ss = substr($ss, 0, length($ss) - $opt{e});
        }
        if ($opt{l}) {
            $ss = substr($ss, 0, $opt{l});
        }
        print $ss;
    }
    print "\n";
    $i++;
}
