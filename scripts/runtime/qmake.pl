#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;

my $cwd = getcwd;
my $err_slack = "pipeline_error";
my $fin_slack = "pipeline_finished";
my $slack_url = "";

my %slack_map = (
    brownd7 => "W013UH0HWUF",
    parejaf => "W01BLNUF7J8"
);

sub HELP_MESSAGE {
    print "Usage: qmake.pl -n [name] -r [numAttempts]\n";
    print "-n: job name\n";
    print "-r: number of attempts (default: 1)\n";
    print "-s: slack notifications\n";
    print "-c: slack channel\n";
    print "-l: parent log dir (default: log)\n";
    exit(1);
}

sub slack {
    my ($slack_channel, $slack_message) = @_;
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
$project_name =~ s:.*/juno/work/bergerm1/Innovation/brownd7/home/::;
$project_name =~ s:.*/lila/data/reis-filho/data/brownd7/home/::;
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
        exec "$qmake $args LOGDIR=$logdir &> $logfile";
    } else {
        waitpid(-1, 0);
        $retcode = $? >> 8;
        my $pipeline_channel_msg = "<\@${slackname}|cal> $project_name :";
        if ($opt{s} && ($retcode == 0 || $n == 0 || $n + 1 == $attempts)) {
            if ($retcode == 0) {
                my $slack_msg = "*COMPLETE* $name :ok_hand:";
                &slack($fin_slack, "$pipeline_channel_msg $slack_msg");
                &slack($opt{c}, $slack_msg) if $opt{c};
            } else {
                my $slack_msg = "*FAILURE* $cwd/$logfile";
                if ($n + 1 == $attempts) {
                    # final attempt
                    $slack_msg = ":-1: $slack_msg";
                    &slack($opt{c}, $slack_msg) if $opt{c};
                }
                &slack($err_slack, "$pipeline_channel_msg $slack_msg");
                sleep 30;
            }
        }
    }
} while ($retcode && ++$n < $attempts);
exit($retcode);
