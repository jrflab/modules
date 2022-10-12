#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use File::Copy;

my $MAKEFILE = <<ENDL;
include modules/Makefile
ENDL

unless (-e "Makefile") {
    open OUT, ">Makefile";
    print OUT $MAKEFILE;
}
close OUT;

unless (-e "project_config.yaml") {
    copy("modules/default_yaml/project_config.yaml", "project_config.yaml") or die "Unable to create project_config.yaml: $!";
}

unless (-e "summary_config.yaml") {
    copy("modules/default_yaml/summary_config.yaml", "summary_config.yaml") or die "Unable to create summary_config.yaml: $!";
}
