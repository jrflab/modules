#!/bin/sh
# gzip old log files
LOCK=~/.gzip_lock
if ! mkdir $LOCK 2> /dev/null; then
    echo "log gzip script is already running"
    exit 1
fi
find /ifs/e63data/reis-filho/data /ifs/e63data/reis-filho/projects/ -name '*.log' -mtime +5 -exec gzip {} \;
rmdir $LOCK
