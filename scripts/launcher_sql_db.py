#!/usr/bin/env python
""" launches the mysql server on a docker node and returns the server ip
"""

import argparse
import MySQLdb
import yaml
from qsub_pbs import Job
import sys
import time

parser = argparse.ArgumentParser(description='launch mysql server on docker node')
parser.add_argument('server_yaml')

args = parser.parse_args()

server_info = yaml.load(open(args.server_yaml, 'r'))
host = server_info['host']
db = server_info['db']

con = None
for attempt in range(2):
    try:
        con = MySQLdb.connect(host=server_info['host'],
                              user=server_info['user'],
                              passwd=server_info['password'],
                              port=server_info['port'],
                              db=server_info['db'])
        break
    except:
        print(("Failed to connect to {} mysql server. Running docker".format(server_info['db'])))
        docker_cmd = "docker run -d -v {}:/var/lib/mysql -p {}:3306 {}".format(server_info['data_dir'],
                                                                               server_info['port'],
                                                                               server_info['docker_repo'])
        print((docker_cmd + "\n"))
        #job = Job(docker_cmd, '-I -l nodes=1:docker -l host={}'.format(server_info['host']))
        #job.run_job()
        #job.wait()
        #time.sleep(90) # wait for mysqld to start

if not con:
    print("Failed to initialize mysql server")
    sys.exit(1)
