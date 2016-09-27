#!/usr/bin/env sh
docker run -d -v /cbio/ski/reis-filho/home/limr/share/usr/fathmm-db:/var/lib/mysql -p 9990:3306 limr/fathmm-db
docker run -d -v /cbio/ski/reis-filho/home/limr/share/usr/chasm3-db:/var/lib/mysql -p 9991:3306 limr/chasm-db
docker run -d -v /cbio/ski/reis-filho/home/limr/share/usr/ensembl-hs-core-85-37-db/:/var/lib/mysql -p 9992:3306 limr/ensembl-hs-core-85-37-db
