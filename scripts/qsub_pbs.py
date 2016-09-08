#!/usr/bin/env python
""" submit jobs to a PBS/torque cluster and wait for them to finish
"""

import argparse
import os
import random
import tempfile
import sys
import signal
import subprocess
import time
import re

DEVNULL = open(os.devnull, 'wb')

class Job:
    _kill_now = False

    def __init__(self, job_script, qsub_args, out_file=None,
                 poll_interval=10):
        self.poll_interval = poll_interval
        self.job_script = job_script
        self.qsub_args = qsub_args
        self.out_file = out_file

        signal.signal(signal.SIGINT, self.__exit_gracefully)
        signal.signal(signal.SIGTERM, self.__exit_gracefully)

    def run_job(self):
        """ qsub job
        """
        self.job_script_file = tempfile.NamedTemporaryFile(mode='w',
                                                           suffix='.sh', delete=False)
        self.job_script_file.write(self.job_script)
        self.job_script_file.close()
        os.chmod(self.job_script_file.name, 0555)

        cmd = "qsub"
        if self.qsub_args is not None:
            cmd += " " + self.qsub_args
        cmd += " " + self.job_script_file.name
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        proc.wait()
        if proc.returncode != 0:
            raise Exception('unable to qsub job: {}'.format(" ".join(cmd)))
        self.job_id = proc.stdout.read().rstrip()
        self.update_qstat()

    def update_qstat(self):
        """ update the job state using qstat
        """
        self.qstat = {}
        proc = subprocess.Popen(['qstat', '-f', self.job_id], stdout=subprocess.PIPE)
        proc.wait()
        job_id_line = proc.stdout.readline().strip()
        if proc.returncode != 0 or self.job_id not in job_id_line:
            raise Exception('unable to qstat job id: {}'.format(self.job_id))
        self.qstat = self._parse_qstat(proc.stdout)

    def _parse_qstat(self, f):
        qstat = {}
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            mo = re.search(r'([^=]+) = (.+)', line)
            if not mo:
                continue
            k = mo.group(1)
            v = mo.group(2)
            while line.endswith(','):
                line = f.readline()
                if not line:
                    break
                line = line.strip()
                v += line
            if v.isdigit():
                qstat[k] = int(v)
            else:
                qstat[k] = v
        return qstat

    def wait(self):
        """ wait for job to finish and return qstat exit_status
        """
        while True:
            self.update_qstat()
            if self.qstat['job_state'] == "C":
                break
            if self.qstat['job_state'] == "H" or self.qstat['job_state'] == "S":
                raise Exception('job halted')
            time.sleep(random.randint(self.poll_interval, self.poll_interval + 60))
            if self._kill_now:
                sys.stderr.write("registering job for deletion\n")
                subprocess.Popen("qdel {}".format(self.job_id), shell=True)
                break
        return self.qstat['exit_status']

    def check_file(self, servers, max_connection_retry=10):
        """ checks the file size of file f on remote servers and returns True
        if it matches local file size
        """
        time.sleep(10)
        if self.out_file is None:
            return True
        for attempt in range(max_connection_retry):
            try:
                local_file_size = subprocess.check_output('stat -c%s "{}"'.format(os.path.abspath(self.out_file)),
                                                        shell=True, stderr=DEVNULL).rstrip()
                break
            except:
                sys.stderr.write("Unable to stat {} locally\n".format(os.path.abspath(self.out_file)))
                time.sleep(10)
            else:
                sys.stderr.write("max retries for local file size check\n")
                return False
        for server in servers:
            for attempt in range(max_connection_retry):
                try:
                    remote_file_size = subprocess.\
                        check_output('ssh {} stat -c%s "{}"'.
                                     format(server, os.path.abspath(self.out_file)),
                                     shell=True, stderr=DEVNULL).rstrip()
                    if remote_file_size != local_file_size:
                        raise ValueError(
                            "{}: remote file size != local file size: {} != {}"
                            .format(server,
                                    remote_file_size,
                                    local_file_size))
                    break
                except ValueError as e:
                    sys.stderr.write(e + "\n")
                    time.sleep(10)
                except:
                    sys.stderr.write("{}: failed remote file size check\n"
                                     .format(server))
                    time.sleep(10)
                else:
                    sys.stderr.write("max connection retries for "
                                     "remote file size check\n")
                    return False
        return True

    def __exit_gracefully(self, signum, frame):
        sys.stderr.write("received interrupt\n")
        self._kill_now = True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='qsub_pbs.py',
                                     description='Submit job to PBS cluster.')
    parser.add_argument('-o', '--out_file',
                        nargs='?', help='file to check for non-zero size')
    parser.add_argument('-i', '--in_file',
                        nargs='?', default = sys.stdin, help='input script file')
    parser.add_argument('--poll_interval', nargs='?',
                        help='qstat polling interval', default=10)
    parser.add_argument('-c', '--check', default=False,
                        action='store_true',
                        help='check this file for non-zero size')
    parser.add_argument('-s', '--servers', nargs='*',
                        default=['gg02', 'gg05', 'gg06'],
                        help='check these servers for '
                        'non-zero output file size')
    parser.add_argument('--args', nargs=argparse.REMAINDER)
    args = parser.parse_args()

    if args.args is not None:
        qsub_args = " ".join(args.args)
    else:
        qsub_args = None

    job = Job(args.in_file.read(), qsub_args, args.out_file, args.poll_interval)
    job.run_job()
    exit_status = job.wait()
    if args.check and exit_status == 0:
        if not job.check_file(args.servers):
            exit_status += 66
    sys.exit(exit_status)
