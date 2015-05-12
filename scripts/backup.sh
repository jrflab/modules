#!/bin/bash

LOCK=~/.backup_lock
LOGFILE=~/.backup.log
if ! mkdir $LOCK 2> /dev/null; then
    echo "backup script is already running"
    exit 1
fi

TMP=`mktemp`;
if mountpoint -q "/mount/limr/zedshared/"; then
    echo "searching for files in /ifs/e63data/reis-filho/"
    cd /ifs/e63data/reis-filho/ && \
    find data projects -type d \
        \( -name bam -o -name tables -o -name alltables -o -name vcf \) \
        ! -path "*/log/*" ! -path "*/tmap/*" ! -path "*/gatk/*" ! -path "*/hydra/*" ! -path "*/bwa/*" ! -path "*/tophat/*" \
        ! -path "*/varscan/*" ! -path "*/mutect/*" ! -path "*/scalpel/*" ! -path "*/som_sniper/*" ! -path "*/rawdata/*" \
        ! -path "*/unprocessed_bam/*" ! -path "*/defuse/*" ! -path "*/chimscan/*" -print0 > ${TMP}

    while [ 1 ]; do
        cd /ifs/e63data/reis-filho/ && \
            rsync --verbose --stats --recursive -a -0 --files-from=${TMP} --log-file=$LOGFILE --prune-empty-dirs ./ /mount/limr/zedshared
        if [ "$?" = "0" ]; then
            echo "rsync complete"
            exit
        else
            echo "rsync failure, retrying in 1 minute..."
            sleep 60
        fi
    done
    rm ${TMP}
fi

rmdir $LOCK
