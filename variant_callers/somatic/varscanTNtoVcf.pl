#!/usr/bin/env perl
# convert varscan tab format to vcf

use strict;
use warnings;

#use experimental 'smartmatch';

use Getopt::Std;
my %opt;
getopts('hf:n:t:', \%opt);

my $usage = <<ENDL;
perl varscanToVcf.pl -f [ref fasta] -t [tumor sample name] -n [normal sample name]
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
my $normal = $opt{n};
my $tumor = $opt{t};

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
##INFO=<ID=SS,Number=1,Type=String,Description="somatic status">
##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">
##INFO=<ID=VARPVAL,Number=1,Type=Float,Description="variant p-value">
##INFO=<ID=SOMPVAL,Number=1,Type=Float,Description="somatic p-value">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=A,Type=Integer,Description="allelic depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="sample depth">
##FORMAT=<ID=VAF,Number=1,Type=Float,Description="variant freq">
ENDL
##PEDIGREE=<Derived=$tumor,Original=$normal>
for my $contig (@contigs) {
    my @F = split /\t/, $contig;
    $vcfHeader .= "##contig=<ID=$F[0],length=$F[1]>\n";
}
$vcfHeader .= "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$tumor\t$normal\n";
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
    my $chrom = $F{chrom};
    my $pos = $F{position};
    my $ref = $F{ref};

    # skip non-somatic/germline
    next unless $F{somatic_status} eq "Somatic" | $F{somatic_status} eq "Germline";

    my $normalGT;
    my $tumorGT;
    my $numDel = 0;
    my @alts;
    # sort alt alleles by length: longest to shortest
    if ($F{var} =~ /[+\-]/) {
        # indels
        my @normalGTs;
        my @tumorGTs;

        push @normalGTs, "0" if $F{normal_gt} =~ /\*/ || $F{normal_gt} eq $F{ref};
        push @tumorGTs, "0" if $F{tumor_gt} =~ /\*/ || $F{tumor_gt} eq $F{ref};
        my $i = 0;
        for my $alt (sort { length $b <=> length $a } (split /\//, $F{var})) {
            $i++;
            push @normalGTs, $i if $F{normal_gt} =~ /\Q$alt\E/;
            push @tumorGTs, $i if $F{tumor_gt} =~ /\Q$alt\E/;
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
        $normalGT = join "/", @normalGTs;
        $normalGT = "0/0" if $normalGT eq "0";
        $normalGT = "1/1" if $normalGT eq "1";
        $tumorGT = join "/", @tumorGTs;
        $tumorGT = "0/0" if $tumorGT eq "0";
        $tumorGT = "1/1" if $tumorGT eq "1";
    } else {
        # snp
        my $alt = $F{var};

        my @normalBases = $iupac{$F{normal_gt}};
        $normalGT = "0/0" if $ref eq $F{normal_gt};
        $normalGT = "1/1" if $alt eq $F{normal_gt};
        $normalGT = "0/1" if $ref ~~ @normalBases && $alt ~~ @normalBases;
        $normalGT = "0/1/2" if scalar(@normalBases) == 3;

        my @tumorBases = $iupac{$F{tumor_gt}};
        $tumorGT = "0/0" if $ref eq $F{tumor_gt};
        $tumorGT = "1/1" if $alt eq $F{tumor_gt};
        $tumorGT = "0/1" if $ref ~~ @tumorBases && $alt ~~ @tumorBases;
        $tumorGT = "0/1/2" if scalar(@tumorBases) == 3;

        push @alts, $alt;
    }
    my $alt = join ",", @alts;

    my $normalR = $F{normal_reads1};
    my $normalA = $F{normal_reads2};
    my $normalVAF = $F{normal_var_freq};
    $normalVAF =~ s/%//;
    $normalVAF /= 100;
    my $normalDP = $normalR + $normalA;
    my $normalAD = "$normalR,$normalA";

    my $tumorR = $F{tumor_reads1};
    my $tumorA = $F{tumor_reads2};
    my $tumorVAF = $F{tumor_var_freq};
    $tumorVAF =~ s/%//;
    $tumorVAF /= 100;
    my $tumorDP = $tumorR + $tumorA;
    my $tumorAD = "$tumorR,$tumorA";

    my $dp = $normalR + $normalA + $tumorR + $tumorR;
    my $ss = $F{somatic_status};
    my $varPval = $F{variant_p_value};
    my $somPval = $F{somatic_p_value};

    my $info = "DP=$dp;SS=$ss;SOMPVAL=$somPval;VARPVAL=$varPval";
    my $normalFormat = "$normalGT:$normalDP:$normalAD:$normalVAF";
    my $tumorFormat = "$tumorGT:$tumorDP:$tumorAD:$tumorVAF";
    print "$chrom\t$pos\t.\t$ref\t$alt\t.\tPASS\t$info\t$format\t$tumorFormat\t$normalFormat\n";
}
