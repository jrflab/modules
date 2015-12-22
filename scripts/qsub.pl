#!/usr/bin/env perl
# qsub wrapper script

use strict;
use warnings;

use Schedule::DRMAAc qw/ :all /;
use File::Temp();
use Cwd qw/ realpath /;
use Cwd;
use Getopt::Std;
use Scalar::Util qw(looks_like_number);

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
        my $nodeFileSize = `ssh -x $node stat -c\%s $cwd/$file`;
        chomp $nodeFileSize;
        if ($nodeFileSize =~ /^stat/ || $nodeFileSize =~ /^ssh/ || !looks_like_number($nodeFileSize)) {
            print "$nodeFileSize\n";
            return 0;
        } elsif ($fileSize != $nodeFileSize) {
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

my $exitCodes = $exitStatus + $aborted + $signaled + $coreDumped;
my $fileStatus = 0;
if ($exitCodes == 0) {
    sleep 20; # wait for file sync
    if ($opt{o} && $opt{o} ne "NULL") {
        my $i = 0;
        while (!check_file(getcwd(), $opt{o})) {
            if ($i++ > 30) {
                if (!-e $opt{o}) {
                    $fileStatus = 66;
                    print STDERR "ERROR $opt{o}: file does not exist\n";
                } else {
                    $fileStatus = 77;
                    print STDERR "ERROR $opt{o}: file sizes do not match across nodes\n";
                    system("rm $opt{o}");
                }
                last;
            }
            sleep 20;
        }
    }
}
$exitCodes += $fileStatus;

# check for zero-size output file and remove it
if ($exitCodes == 0) {
    my $i = 0;
    while ($opt{c} && $opt{o} && $opt{o} ne "NULL" && !-s $opt{o}) {
        if ($i++ > 30) { 
            $fileStatus = 99;
            print STDERR "ERROR $opt{o}: file is size 0\n";
            system("rm $opt{o}");
            last;
        }
        sleep 20;
    }

}
$exitCodes += $fileStatus;

exit $exitCodes;
