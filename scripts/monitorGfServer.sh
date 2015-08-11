#!/bin/bash
UP=$(pgrep -u limr gfServer | wc -l);
if [ "$UP" -ne 1 ];
then
    echo "gfServer is down.";
    /home/limr/share/usr/bin/gfServer start localhost 88878 -stepSize=5 -log=/home/limr/.blatserver.log /home/limr/share/reference/GATK_bundle/2.3/human_g1k_v37.2bit
else
    echo "All is well.";
fi
