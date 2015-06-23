#!/bin/sh
# gzip old log files
find /ifs/e63data/reis-filho/data /ifs/e63data/reis-filho/projects/ -name '*.log' -mtime +5 -exec gzip {} \;
