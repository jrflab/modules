#!/usr/bin/env perl
# qsub wrapper script

use strict;
use warnings;

use Schedule::DRMAAc qw/ :all /;
use File::Temp();
use Cwd qw/ realpath /;
use Cwd;
use Getopt::Std;

my %opt;
getopts('hco:', \%opt);

my $usage = <<ENDL;
Usage: perl qsub.pl -h -- [qsub args]
    -o [file]: check file for non-zero size
    -c: check file across nodes for same file size
ENDL

sub check_file {
    my ($cwd, $file) = @_;
    my @nodes = qw/e01 e02 e03 e04 e05 e06/;
    my $fileSize = `stat -c\%s $cwd/$file`;
    chomp $fileSize;
    #print "checking $cwd/$file on nodes ($fileSize)\n";
    for my $node (@nodes) {
        my $nodeFileSize = `ssh $node stat -c\%s $cwd/$file`;
        chomp $nodeFileSize;
        if ($nodeFileSize eq "" || $fileSize != $nodeFileSize) {
            #print "$node: file size does not match: $fileSize != $nodeFileSize\n";
            return 0;
        }
    }
    #print "$cwd/$file: all file sizes match\n";
    return 1;
}

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}

HELP_MESSAGE if $opt{h};

my $scriptFile = File::Temp->new(TEMPLATE => 'tempXXXXX', DIR => '/home/limr/share/tmp', SUFFIX => '.sge');

my $args = join " ", @ARGV;
while (<STDIN>) {
    print $scriptFile $_;
}
close $scriptFile;

my ($error, $diagnosis) = drmaa_init(undef);
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, my $jt, $diagnosis) = drmaa_allocate_job_template();
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, $diagnosis) = drmaa_set_attribute($jt, $DRMAA_REMOTE_COMMAND, $scriptFile->filename);
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, $diagnosis) = drmaa_set_attribute($jt, $DRMAA_NATIVE_SPECIFICATION, $args);
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, $diagnosis) = drmaa_set_attribute($jt, $DRMAA_WD, getcwd());
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, my $jobid, $diagnosis) = drmaa_run_job($jt);
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, $diagnosis) = drmaa_delete_job_template($jt);
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

sub signalHandler {
    my ($error, $diagnosis) = drmaa_control($jobid, $DRMAA_CONTROL_TERMINATE);
    die drmaa_strerror($error) . "\n" . $diagnosis if $error;
    die "Received interrupt: terminating job\n";
}

$SIG{INT} = \&signalHandler;
$SIG{TERM} = \&signalHandler;

# loop to give a chance to receive sigint/sigterms
my $stat;
do {
    ($error, my $jobidOut, $stat, $diagnosis) = drmaa_wait($jobid, 10);
} until ($error != $DRMAA_ERRNO_EXIT_TIMEOUT);

# pull all exit-related codes
($error, my $exitStatus, $diagnosis) = drmaa_wexitstatus($stat);
die drmaa_strerror($error) . "\n" . $diagnosis if $error;
($error, my $aborted, $diagnosis) = drmaa_wifaborted( $stat );
die drmaa_strerror($error) . "\n" . $diagnosis if $error;
($error, my $signaled, $diagnosis ) = drmaa_wifsignaled( $stat );
die drmaa_strerror($error) . "\n" . $diagnosis if $error;
($error, my $coreDumped, $diagnosis ) = drmaa_wcoredump( $stat );
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, $diagnosis) = drmaa_exit();
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

my $fileStatus = 0;
if ($opt{c} && $opt{o} && (!-e $opt{o} || !-s $opt{o})) {
    sleep 60; # wait for file system to update
    system("rm $opt{o}");
    #print "File not removed\n" if (-e $opt{o});
    $fileStatus = 99;
    print STDERR "$opt{o}: file is size 0\n";
}

sleep 20; # wait for file sync
if ($opt{o} && $opt{o} ne "NULL") {
    my $i = 0;
    while (!check_file(getcwd(), $opt{o})) {
        if ($i++ > 30) {
            $fileStatus = 77;
            print STDERR "file sizes do not match across nodes\n";
        }
        sleep 10;
    }
}

exit $exitStatus + $aborted + $signaled + $coreDumped + $fileStatus;
