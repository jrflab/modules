#!/bin/bash
UP=$(pgrep -u limr mysqld | wc -l);
if [ "$UP" -ne 1 ];
then
    echo "MySQL is down.";
    mysqld --defaults-file=/home/limr/share/usr/mysql/my.cnf
else
    echo "All is well.";
fi
