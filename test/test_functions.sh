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
        diff -q <(samtools view $1) <(samtools view $2) > /dev/null && \
        echo "success, files the same: ${1} ${2}" \
    ) || (echo "failed, files differ: ${1} ${2}" && exit 1)
}
