#!/usr/bin/env perl

use Cwd;
use File::Copy;

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
include modules/Makefile
ENDL

unless (-d "modules") {
    system "git clone --recursive git\@github.com:jrflab/modules.git -b master";
}

unless (-e "Makefile") {
    open OUT, ">Makefile";
    print OUT $MAKEFILE;
}
close OUT;
unless (-e "project_config.yaml") {
    copy("modules/default_project_config.yaml", "project_config.yaml") or die "Unable to create project_config.yaml: $!";
}
unless (-e "summary_config.yaml") {
    copy("modules/default_summary_config.yaml", "summary_config.yaml") or die "Unable to create summary_config.yaml: $!";
}
unless (-e "sample_attr.yaml") {
    copy("modules/default_sample_attr.yaml", "sample_attr.yaml") or die "Unable to create sample_attr.yaml: $!";
}
close OUT;
system "git add Makefile";
system "git add project_config.yaml";
system "git add summary_config.yaml";
system "git commit -m 'makefile, summary_config.yaml, project_config.yaml'";
system "git push --set-upstream origin master";

