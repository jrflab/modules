#!/usr/bin/env perl
# convert varscan tab format to vcf

use strict;
use warnings;

use experimental 'smartmatch';

use Getopt::Std;
my %opt;
getopts('hf:n:t:', \%opt);

my $usage = <<ENDL;
perl varscanToVcf.pl -f [ref fasta] -s [sample name]
ENDL

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}

my %iupac = (
    A => ["A"],
    C => ["C"],
    G => ["G"],
    T => ["T"],
    W => ["A", "T"],
    S => ["C", "G"],
    M => ["A", "C"],
    K => ["G", "T"],
    R => ["A", "G"],
    Y => ["C", "T"],
    B => ["C", "G", "T"],
    D => ["A", "G", "T"],
    H => ["A", "C", "T"],
    V => ["A", "C", "G"]
);


HELP_MESSAGE if $opt{h};

my $now = localtime;
my $sample = $opt{s};

open REF, $opt{f} . ".fai" or die "Unable to open reference index file\n";

sub getRefSeq {
    my ($chrom, $pos) = @_;
    open(SAMTOOLS, "samtools faidx $opt{f} $chrom:$pos-$pos |");
    <SAMTOOLS>;
    my $seq = <SAMTOOLS>;
    chomp $seq;
    close(SAMTOOLS);
    $seq;
}

my @contigs = <REF>;

my $vcfHeader = <<ENDL;
##fileformat=VCFv4.1
##fileDate=$now
##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">
##INFO=<ID=VARPVAL,Number=1,Type=Float,Description="variant p-value">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=A,Type=Integer,Description="allelic depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="sample depth">
##FORMAT=<ID=VAF,Number=1,Type=Float,Description="variant freq">
ENDL
for my $contig (@contigs) {
    my @F = split /\t/, $contig;
    $vcfHeader .= "##contig=<ID=$F[0],length=$F[1]>\n";
}
$vcfHeader .= "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample\n";
print $vcfHeader;

my $format = "GT:DP:AD:VAF";

my $headerLine = <>;
chomp $headerLine;
$headerLine =~ s/^#//;
my @header = split/\t/, $headerLine;

while (<>) {
    chomp;
    my @F = split /\t/;
    my %F = map { $_ => shift @F } @header;
    my $chrom = $F{Chrom};
    my $pos = $F{Position};
    my $ref = $F{Ref};
    my $alt = $F{Alt};

    my $gt;
    my $numDel = 0;
    my @alts;
    # sort alt alleles by length: longest to shortest
    if ($F{Var} =~ /[+\-]/) {
        # indels
        my @gts;

        push @gts, "0" if $F{Var} =~ /\*/ || $F{Var} eq $F{Ref};
        my $i = 0;
        for my $alt (sort { length $b <=> length $a } (split /\//, $F{Var})) {
            $i++;
            push @gts, $i if $F{gt} =~ /\Q$alt\E/;
            if ($alt =~ /^\+/) {
                # insertion
                $alt =~ s/^\+//;
                $alt = $ref . $alt;
            } else {
                # deletion
                $alt =~ s/^\-//;
                $numDel++;
                if ($numDel == 1) {
                    # longest deletion is absorbed into reference
                    $ref .= $alt;
                }
                # trim alt from reference sequence to create new alt
                my $newAlt = $ref;
                $newAlt =~ s/\Q$alt\E$//;
                $alt = $newAlt;
            }
            push @alts, $alt;
        }
        $gt = join "/", @gts;
        $gt = "0/0" if $gt eq "0";
        $gt = "1/1" if $gt eq "1";
    } else {
        # snp
        my $alt = $F{Var};

        my @bases = $iupac{$F{Var}};
        $gt = "0/0" if $ref eq $F{Var};
        $gt = "1/1" if $alt eq $F{Var};
        $gt = "0/1" if $ref ~~ @bases && $alt ~~ @bases;
        $gt = "0/1/2" if scalar(@bases) == 3;

        push @alts, $alt;
    }
    my $alt = join ",", @alts;

    my $R = $F{reads1};
    my $A = $F{reads2};
    my $VAF = $F{var_freq};
    $VAF =~ s/%//;
    $VAF /= 100;
    my $DP = $R + $A;
    my $AD = "$R,$A";

    my $varPval = $F{variant_p_value};

    my $info = "DP=$DP;VARPVAL=$varPval";
    my $sFormat = "$GT:$DP:$AD:$VAF";
    print "$chrom\t$pos\t.\t$ref\t$alt\t.\tPASS\t$info\t$format\t$sFormat\n";
}
