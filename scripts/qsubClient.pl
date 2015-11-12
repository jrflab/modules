#!/usr/bin/env perl

use strict;
use warnings;

use Cwd qw/ realpath /;
#use IO::Socket qw/ AF_UNIX SOCK_STREAM /;
use IO::Socket::INET;
use IO::Socket::Timeout;
use Errno qw(ETIMEDOUT EWOULDBLOCK);
use Path::Class qw/ file /;
use File::Temp();
use Cwd;
use Storable qw/ nfreeze /;

use Getopt::Std;
my %opt;
getopts('hco:', \%opt);

my $usage = <<ENDL;
Usage: perl qsubClient.pl -h -- [qsub args]
    -c : check file for non-zero size
    -o [file]: output file (checked for uniform size among all nodes)
ENDL

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}


HELP_MESSAGE if $opt{h};

my $args = join " ", @ARGV;
my $socketPath = file($opt{s});
#my $client = IO::Socket->new(
#Domain => AF_UNIX,
#Type => SOCK_STREAM,
#Peer => $socketPath,
#Timeout => 30,
#) or die("Can't connect to server socket: $!\n");
#

my $maxRetry = 10;
my $client;
my $i = 0;
my $port = 34999;
while (!$client && $i < $maxRetry) {
    $client = IO::Socket::INET->new(
        PeerHost => 'localhost',
        PeerPort => "$port",
        Proto => 'tcp',
    );
    unless ($client) {
        print "Can't connect to server socket: $!\n";
        sleep 10;
        print "Retrying... (attempt $i)\n";
        $i++;
    }
}
die "Can't connect to server\n" unless ($client);
eval 'END { close $client } 1' or die $@;

IO::Socket::Timeout->enable_timeouts_on($client);
$client->read_timeout(200);
$client->write_timeout(200);

#print "Connected to server $socketPath\n";
#print "Sending server args: $args\n";
print $client $args . "\n";
#print "Sending server cwd: " . getcwd() . "\n";
print $client getcwd() . "\n";

my $scriptFile = File::Temp->new(TEMPLATE => 'tempXXXXX', DIR => '/home/limr/share/tmp', SUFFIX => '.sge', UNLINK => 0);
chmod 0644, $scriptFile;
while (my $line = <STDIN>) {
    print $scriptFile $line;
}
close $scriptFile;
#print "Sending server script: " . $scriptFile->filename . "\n";
print $client $scriptFile->filename . "\n";

if ($opt{o}) {
    print $client $opt{o} . "\n";
} else {
    print $client "NULL\n";
}

my $exitCode = -1;
while (<$client>) {
    if (! $_ && (0+$! == ETIMEDOUT || 0+$! == EWOULDBLOCK)) {
        print STDERR "Connection timeout\n";
        last;
    } elsif (/^Error:/) {
        print STDERR;
        $exitCode = -1;
        last;
    } elsif (/^Code:/) {
        ($exitCode) = $_ =~ /Code: (\d+)/;
        last;
    }
}

if ($opt{c} && $opt{o} && (!-e $opt{o} || !-s $opt{o})) {
    sleep 60; # wait for file system to update
    system("rm $opt{o}");
    #print "File not removed\n" if (-e $opt{o});
    die "$opt{o}: file is size 0";
}



exit $exitCode;
