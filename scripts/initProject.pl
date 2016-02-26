#!/usr/bin/env perl

use Cwd;

# create repo
my $repoName = getcwd;
$repoName =~ s:.*/projects/::;
$repoName =~ s:.*/data/::;
$repoName =~ s:.*kinglab/::;
$repoName =~ s:/:_:g;

system "bb create -c --owner jrflab --protocol ssh $repoName";
system "git init";
system "git remote add origin git\@bitbucket.org:jrflab/$repoName.git";


my $CONFIG = <<ENDL;
---
ref: b37

aligner: bwamem #tophat hisat bwa bowtie tmap
bam_chr1_base_recal: true
bam_dup_type: markdup
bam_no_filter: false
bam_no_recal: false
bam_no_realn: false
bam_no_sort: false
bam_fix_rg: false
bam_phred64: false
bam_reprocess: false

vcf_post_ann_filter_expression: ExAC_AF > 0.05

# targets_file: intervals.bed
exome: true

# gatk options
gatk_hard_filter_snps: true
gatk_pool_snp_recal: false

qsub_priority: -800
...
ENDL

my $MAKEFILE = <<ENDL;
include modules/Makefile
ENDL

unless (-d "modules") {
    system "git clone git\@github.com:jrflab/modules.git -b master";
}

unless (-e "Makefile") {
    open OUT, ">Makefile";
    print OUT $MAKEFILE;
}
close OUT;
unless (-e "config.inc") {
    open OUT, ">config.yaml";
    print OUT $CONFIG;
}
close OUT;
system "git add Makefile";
system "git add config.yaml";
system "git commit -m 'makefile, config.yaml'";
system "git push --set-upstream origin master";

