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
    system "git clone git\@github.com:jrflab/modules.git -b master";
}

unless (-e "Makefile") {
    open OUT, ">Makefile";
    print OUT $MAKEFILE;
}
close OUT;
unless (-e "project_config.yaml") {
    copy("modules/default_project_config.yaml", "project_config.yaml") or die "Unable to create project_config.yaml: $!";
}
close OUT;
system "git add Makefile";
system "git add project_config.yaml";
system "git commit -m 'makefile, project_config.yaml'";
system "git push --set-upstream origin master";

