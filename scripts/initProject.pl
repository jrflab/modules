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


my $MAKEFILE = <<ENDL;
export REF = hg19
DUP_TYPE = markdup

#TARGETS_FILE = intervals.bed
EXOME = true

# gatk options
HARD_FILTER_SNPS = true

QSUB_PRIORITY = -800

include modules/Makefile
ENDL

unless (-d "modules") {
    system "git clone git@github.com:raylim/modules.git -b stable";
}

unless (-e "Makefile") {
    open OUT, ">Makefile";
    print OUT $MAKEFILE;
}
system "git add Makefile";
system "git commit -m 'makefile'";
system "git push --set-upstream origin master";

