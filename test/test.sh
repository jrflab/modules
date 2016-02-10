#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Use gnu corutils gstat for OSX
# http://apple.stackexchange.com/questions/69223
# http://stackoverflow.com/questions/3466166
if [ "$(uname)" == "Darwin" ]; then
    STAT=gstat
else
    STAT=stat
fi

compare_files() {
    diff -q <(${STAT} -c %s $1) <(${STAT} -c %s $2) > /dev/null \
        && echo "success, files the same: ${1} ${2}" \
        || (echo "failed, files differ: ${1} ${2}" && exit 1)
}

compare_bams() {
    (
        samtools quickcheck $1 && \
        samtools quickcheck $2 && \
        diff -q <(samtools view $1) <(samtools view $2) > /dev/null && \
        echo "success, files the same: ${1} ${2}" \
    ) || (echo "failed, files differ: ${1} ${2}" && exit 1)
}

echo "TEST TARGETED SEQUENCING test_data/tseq"
cd ${DIR}/test_data/tseq
echo "==bwamem=="
echo "run"
make bwamem && echo "success" || (echo "failed" && exit 1)
echo "check files"
compare_bams bam/test1T.bam output/bam/test1T.bam

exit 0
