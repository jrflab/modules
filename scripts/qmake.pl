#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;

my $cwd = getcwd;
my $err_slack = "pipeline_error";
my $fin_slack = "pipeline_finished";

my %slack_map = (
    limr => "raylim",
    debruiji => "debruiji",
    brownd7 => "brownd7",
    lees19 => "lees19",
    ferrandl => "ferrandl",
    dacruzpa => "dacruzpa"
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
    my $slack_url = "\$'https://jrflab.slack.com/services/hooks/slackbot?token=2TWPiY9Hu4EUteoECqCEfYAZ&channel=%23$slack_channel'";
    system "curl --data ' $slack_message' $slack_url &> /dev/null";
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
        }
        
        if ($username eq "") {
            my $pipeline_channel_msg = "\@${slackname} $project_name :";
            if ($opt{s} && ($retcode == 0 || $n == 0 || $n + 1 == $attempts)) {
             if ($retcode == 0) {
                 # op success
                 my $slack_msg = "*FAILURE* $cwd/$logfile";
                 $slack_msg = "$slack_msg :troll:";
                 &slack($opt{c}, $slack_msg) if $opt{c};
             } else {
                 # op failure
                 my $slack_msg = "*FAILURE* $cwd/$logfile";
                 if ($n + 1 == $attempts) {
                     # final attempt
                     $slack_msg = "$slack_msg :troll:";
                     &slack($opt{c}, $slack_msg) if $opt{c};
                 }
                 &slack($err_slack, "$pipeline_channel_msg $slack_msg");
                 # wait a bit before retrying to allow cleanup
                 sleep 30;
            }
          }
        } else {
            my $pipeline_channel_msg = "\@${slackname} $project_name :";
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
                     $slack_msg = "$slack_msg :troll:";
                     &slack($opt{c}, $slack_msg) if $opt{c};
                 }
                 &slack($err_slack, "$pipeline_channel_msg $slack_msg");
                 # wait a bit before retrying to allow cleanup
                 sleep 30;
            }
          }
        }
    }
} while ($retcode && ++$n < $attempts);
exit($retcode);
