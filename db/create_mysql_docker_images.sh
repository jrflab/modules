docker run --name fathmm-db \
    -v ~/downloads/fathmm:/docker-entrypoint-initdb.d \
    -e MYSQL_ROOT_PASSWORD=fathmm \
    -e MYSQL_DATABASE=fathmm \
    -e MYSQL_USER=fathmm \
    -e MYSQL_PASSWORD=fathmm \
    -p 9990:3306 \
    -v ~/share/usr/fathmm-db:/var/lib/mysql \
    -d mysql:5.7.14
docker run --name chasm-db \
    -v ~/downloads/chasm:/docker-entrypoint-initdb.d \
    -e MYSQL_ROOT_PASSWORD=chasm \
    -e MYSQL_DATABASE=CHASM \
    -e MYSQL_USER=chasm \
    -e MYSQL_PASSWORD=chasm \
    -p 9991:3306 \
    -v ~/share/usr/chasm3-db:/var/lib/mysql \
    -d mysql:5.7.14
docker run --name ensembl-hs-core-85-37-db \
    -v ~/downloads/homo_sapiens_core_85_37:/docker-entrypoint-initdb.d \
    -e MYSQL_ROOT_PASSWORD=embl \
    -e MYSQL_DATABASE=homo_sapiens_core_85_37 \
    -e MYSQL_USER=embl \
    -e MYSQL_PASSWORD=embl \
    -e MYSQL_ALLOW_EMPTY_PASSWORD=yes \
    -p 9992:3306 \
    -v ~/share/usr/ensembl-hs-core-85-37-db:/var/lib/mysql \
    -d mysql:5.5
