#!/usr/bin/env perl
# qsub daemon 

use strict;
use warnings;

use Proc::Daemon;
use threads;
use threads::shared;
use Schedule::DRMAAc qw/ :all /;
use File::Temp();
#use IO::Socket qw/ AF_UNIX SOCK_STREAM SOMAXCONN/;
use IO::Socket::INET;
use Path::Class qw/ file /;
use Cwd;

Proc::Daemon::Init;

sub done { exit(0); }
$SIG{TERM} = \&done;
$SIG{INT} = \&done;

# process job
sub job_thread {
    my ($client, $jobid) = @_;
    my $stat;
    my $error;
    my $diagnosis;
    do {
        # delete jobs of disconnected clients
        unless ($client->connected) {
            print "Client disconnected, terminating job\n";
            my ($error, $diagnosis) = drmaa_control($jobid, $DRMAA_CONTROL_TERMINATE);
            print "Error: " .drmaa_strerror($error) . " : " . $diagnosis . "\n" if $error;
        }
        ($error, my $jobidOut, $stat, my $rusage, $diagnosis) = drmaa_wait($jobid, 20);
        print $client "Keep Alive\n";
    } until ($error != $DRMAA_ERRNO_EXIT_TIMEOUT);
    print $client "Error: " .drmaa_strerror($error) . " : " . $diagnosis . "\n" if $error;

    # tell client job is complete
    # pull all exit-related codes
    ($error, my $exitStatus, $diagnosis) = drmaa_wexitstatus($stat);
    print $client "Error: " .drmaa_strerror($error) . " : " . $diagnosis . "\n" if $error;
    ($error, my $aborted, $diagnosis) = drmaa_wifaborted( $stat );
    print $client "Error: " . drmaa_strerror($error) . " : " . $diagnosis . "\n" if $error;
    ($error, my $signaled, $diagnosis ) = drmaa_wifsignaled( $stat );
    print $client "Error: " . drmaa_strerror($error) . " : " . $diagnosis . "\n" if $error;
    ($error, my $coreDumped, $diagnosis ) = drmaa_wcoredump( $stat );
    print $client "Error: " . drmaa_strerror($error) . " : " . $diagnosis . "\n" if $error;

    print $client "Code: " . ($exitStatus + $aborted + $signaled + $coreDumped);
};

#my $socketPath = file($ARGV[0]);
#my $server = new IO::Socket->new (
#Domain => AF_UNIX,
#Type => SOCK_STREAM,
#Local => $socketPath,
#Listen => 3000,
#) or die "Unable to create server socket: $!\n";
#eval 'END { unlink $socketPath } 1' or die $@;

my $server = IO::Socket::INET->new(
    LocalHost => 'localhost',
    LocalPort => '34383',
    Proto => 'tcp',
    Listen => 1000,
    Reuse => 1,
);

my ($error, $diagnosis) = drmaa_init(undef);
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

use threads ('yield',
    'stack_size' => 64*4096,
    'exit' => 'threads_only',
    'stringify');

while (my $client = $server->accept()) {
    print "Client Connected\n";
    my $clientArgs = <$client>;
    chomp $clientArgs;
    print "Client args: $clientArgs\n";
    my $clientCwd = <$client>;
    chomp $clientCwd;
    print "Client wd: $clientCwd\n";
    my $scriptFile = <$client>;
    chomp $scriptFile;
    print "Client script file: $scriptFile\n";

    ($error, my $jt, $diagnosis) = drmaa_allocate_job_template();
    print $client "Error: " . drmaa_strerror($error) . " : " . $diagnosis . "\n" and next if $error;
    
    ($error, $diagnosis) = drmaa_set_attribute($jt, $DRMAA_REMOTE_COMMAND, "$scriptFile");
    print $client "Error: " . drmaa_strerror($error) . " : " . $diagnosis . "\n" and next if $error;

    ($error, $diagnosis) = drmaa_set_attribute($jt, $DRMAA_NATIVE_SPECIFICATION, $clientArgs);
    print $client "Error: " . drmaa_strerror($error) . " : " . $diagnosis . "\n" and next if $error;

    ($error, $diagnosis) = drmaa_set_attribute($jt, $DRMAA_WD, $clientCwd . "/");
    print $client "Error: " . drmaa_strerror($error) . " : " . $diagnosis . "\n" and next if $error;

    ($error, my $jobid, $diagnosis) = drmaa_run_job($jt);
    print $client "Error: " . drmaa_strerror($error) . " : " . $diagnosis . "\n" and next if $error;

    ($error, $diagnosis) = drmaa_delete_job_template($jt);
    print $client "Error: " . drmaa_strerror($error) . " : " . $diagnosis . "\n" and next if $error;

    unlink $scriptFile;

    async(\&job_thread, $client, $jobid)->detach;
}

($error, $diagnosis) = drmaa_exit();
die drmaa_strerror($error) . "\n" . $diagnosis if $error;
