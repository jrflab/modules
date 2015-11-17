#!/bin/bash

LOCK=~/.backup_lock
LOGFILE=~/.backup.log
if ! mkdir $LOCK 2> /dev/null; then
    echo "backup script is already running"
    exit 1
fi

TMP=`mktemp`;
TOPDIR=/ifs/e63data/reis-filho
if mountpoint -q "/mount/limr/zedshared/"; then
    while [ 1 ]; do
        echo "searching for files in $TOPDIR"
        cd $TOPDIR
        'ls' data/*/*/bam/*.bam* projects/*/bam/*.bam* data/*/wgs*/fastq/*.fastq.gz | \
            rsync --verbose --stats --recursive -a --files-from=- --log-file=$LOGFILE --prune-empty-dirs ./ /mount/limr/zedshared
        if [ "$?" = "0" ]; then
            echo "rsync complete"
            exit
        else
            echo "rsync failure, retrying in 1 minute..."
            sleep 60
        fi
    done
fi

rmdir $LOCK
