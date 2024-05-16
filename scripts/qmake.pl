#!/usr/bin/env perl
# wrapper script for qmake to remove newlines

use strict;
use warnings;
use Cwd;

my $cwd = getcwd;
my $err_slack = "pipeline_error";
my $fin_slack = "pipeline_finished";

my %slack_map = (
    brownd7 => "W013UH0HWUF",
    selenicp => "W0142HA5LNA",
    dacruzpa => "W01BT68MSSD",
    parejaf => "W01BLNUF7J8",
    zhuy1 => "W013UH382P9",
    peix => "W0147TPN3E1",
    issabhas => "U01V8R1RKQU",
    xiaoy => "U01C8MPBSH5",
    giacomf1 => "U06SW7W6D44"
);

sub HELP_MESSAGE {
    print "Usage: qmake.pl -n [name] -m -r [numAttempts]\n";
    print "-m: e-mail notifications\n";
    print "-s: slack notifications\n";
    print "-c: slack channel\n";
    print "-r: number of attempts (default: 1)\n";
    print "-l: parent log dir (default: log)\n";
    print "-n: job name\n";
    exit(1);
}

sub slack {
    my ($slack_channel, $slack_message) = @_;
    my $slack_url = "";
    if ($slack_channel eq "pipeline_error") {
    	$slack_url = $ENV{SLACK_URL_ERR};
    } elsif ($slack_channel eq "pipeline_finished") {
    	$slack_url = $ENV{SLACK_URL_FIN};
    }
    system "curl -X POST -H 'Content-type: application/json' --data '{\"text\":\"$slack_message\"}' $slack_url &> /dev/null";
}


use File::Basename;
use File::Glob ':glob';
use File::Path;

use Getopt::Std;

my %opt;

getopts('n:smr:l:c:', \%opt);

my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);
my $slackname = $slack_map{$username} || $username;
my $project_name = $cwd;
$project_name =~ s:.*/projects/::;
$project_name =~ s:.*/data/::;
$project_name =~ s:.*kinglab/::;
$project_name =~ s:/:_:g;
my $attempts = 1;
my $name = "qmake";
my $logparent = "log";
$attempts = $opt{r} if defined $opt{r};
$name = $opt{n} if defined $opt{n};
$logparent = $opt{l} if defined $opt{l};

my $qmake = shift @ARGV; 

my $args = join " ", @ARGV;

# makefile processing
=pod
my $orig_args = $args;
$args =~ s;-f (\S+);"-f " . dirname($1) . "/." . basename($1) . ".tmp";e;
my $optf = $1;
my @makefiles;
if (defined $optf) {
    push @makefiles, $optf;
} else {
    if ($args =~ /--/) {
        $args .= " -f .Makefile.tmp";
    } else {
        $args .= "-- -f .Makefile.tmp";
    }
    push @makefiles, "Makefile";
}
do {
    my $makefile = glob(shift(@makefiles));
    
    open IN, "<$makefile" or die "Unable to open $makefile\n";
    my $tmpfile = glob(dirname($makefile) . "/." . basename($makefile) . ".tmp");
    open OUT, ">$tmpfile" or die "Unable to open $tmpfile\n";
    while (<IN>) {
        s/\\\n$//;
        if (!/^include \S+\.tmp/ && s;^include (\S+);"include " . dirname($1) . "/." . basename($1) . ".tmp";e) {
            push @makefiles, $1;
        }
        print OUT $_;
    }
} until (scalar @makefiles == 0);
=cut

my $n = 0;
my $retcode;
do {
    my $logdir = "$logparent/$name";
    my $logfile = "$logdir.log";
    my $i = 0;
    while (-e $logdir || -e $logfile) {
        $logdir = "log/$name.$i";
        $logfile = "$logdir.log";
        $i++;
    }
    mkpath $logdir;
    my $pid = fork;
    if ($pid == 0) {
        #print "$qmake $args &> $logfile\n";
        exec "$qmake $args LOGDIR=$logdir &> $logfile";
    } else {
        my $mail_msg = "Command: $qmake $args\n";
        $mail_msg .=  "Attempt #: " . ($n + 1) . " of $attempts\n";
        $mail_msg .=  "Hostname: " . $ENV{HOSTNAME}. "\n";
        $mail_msg .=  "PID: $pid\n";
        $mail_msg .=  "Dir: $cwd\n";
        $mail_msg .=  "Log dir: $cwd/$logdir\n";
        $mail_msg .=  "Log file: $cwd/$logfile\n";

        if ($opt{m} && ($n == 0 || $n == 1 || $n + 1 == $attempts)) {
            my $mail_subject = "$name: job started ($cwd)";
            $mail_subject .= " Attempt " . ($n + 1) if $n > 0; 
            #open(MAIL, "| mail -s '$mail_subject' $start_email_addrs");
            #print MAIL "$mail_msg";
            #close MAIL;
        }
        waitpid(-1, 0);
        $retcode = $? >> 8; # shift bits to get the real return code
        if ($opt{m} && ($retcode == 0 || $n == 0 || $n == 1 || $n + 1 == $attempts)) {
            #my $addrs = ($retcode > 0)? $err_email_addrs : $fin_email_addrs;
            my $mail_subject = "[$retcode] $name: job finished ($cwd)";
            if ($n + 1 == $attempts) {
                $mail_subject = "**FINAL** $mail_subject";
            }
            $mail_subject .= " Attempt " . ($n + 1) if $n > 0; 
            #open(MAIL, "| mail -s '$mail_subject' $addrs");
            #print MAIL "Return code: $retcode\n";
            #print MAIL "$mail_msg";
            #close MAIL;
        }

        my $pipeline_channel_msg = "<\@${slackname}|cal> $project_name :";
        if ($opt{s} && ($retcode == 0 || $n == 0 || $n + 1 == $attempts)) {
            if ($retcode == 0) {
                # op success
                my $slack_msg = "*COMPLETE* $name :the_horns:";
                &slack($fin_slack, "$pipeline_channel_msg $slack_msg");
                &slack($opt{c}, $slack_msg) if $opt{c};
            } else {
                # op failure
                my $slack_msg = "*FAILURE* $cwd/$logfile";
                if ($n + 1 == $attempts) {
                    # final attempt
                    $slack_msg = ":-1: $slack_msg";
                    &slack($opt{c}, $slack_msg) if $opt{c};
                }
                &slack($err_slack, "$pipeline_channel_msg $slack_msg");
                # wait a bit before retrying to allow cleanup
                sleep 30;
            }
        }
    }
} while ($retcode && ++$n < $attempts);
exit($retcode);
