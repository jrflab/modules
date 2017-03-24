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

diff_files() {
    diff -q $1 $2 > /dev/null \
        && echo "success, files the same: ${1} ${2}" \
        || (echo "failed, files differ: ${1} ${2}" && exit 1)
}

compare_vcf() {
    python modules/vcf_tools/compare_vcf.py --ignore_info -q $1 $2 &> /dev/null \
        && echo "success, files the same: ${1} ${2}" \
        || (echo "failed, files differ: ${1} ${2}" && exit 1)
}

upload_image() {
    local image

    # convert pdf to png for uploading to clbin
    if [[ $1 =~ .*pdf$ ]]
    then
        convert -density 150 $1 -quality 90 ${1/pdf/png}
        image=${1/pdf/png}
    else
        image=$1
    fi

    # upload image
    curl -F "clbin=@$image" https://clbin.com
}

compare_images() {
    # first argument will be uploaded if different, so 2nd argument should be
    # the test image
    diff -q <(${STAT} -c %s $1) <(${STAT} -c %s $2) > /dev/null \
        && echo "success, images the same: ${1} ${2}" \
        || (echo "failed, images differ: ${1} ${2}" &&
            upload_image $1 &&
            exit 1)
}

compare_bams() {
    (
        diff -q <(samtools view $1) <(samtools view $2) > /dev/null && \
        echo "success, files the same: ${1} ${2}" \
    ) || (echo "failed, files differ: ${1} ${2}" && exit 1)
}

compare_cncf() {
    python modules/copy_number/compare_facets_cncf.py $1 $2
}

