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
import datetime
import re

DEVNULL = open(os.devnull, 'wb')


def time2secs(s):
    x = time.strptime(s, '%H:%M:%S')
    return datetime.timedelta(hours=x.tm_hour,
                              minutes=x.tm_min,
                              seconds=x.tm_sec).total_seconds()


def secs2time(s):
    m, s = divmod(s, 60)
    h, m = divmod(m, 60)
    return "%02d:%02d:%02d" % (h, m, s)


def human2bytes(s):
    """
    >>> human2bytes('1M')
    1048576
    >>> human2bytes('1G')
    1073741824
    """
    symbols = ('B', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y')
    num, letter = re.match(r"(\d+)([A-Za-z])", s).group(1, 2)
    letter = letter.upper()
    assert num.isdigit() and letter in symbols
    num = float(num)
    prefix = {symbols[0]: 1}
    for i, s in enumerate(symbols[1:]):
        prefix[s] = 1 << (i+1)*10
    return int(num * prefix[letter])


def bytes2human(n, format="%(value)i%(symbol)s"):
    """
    >>> bytes2human(10000)
    '9K'
    >>> bytes2human(100001221)
    '95M'
    """
    symbols = ('B', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y')
    prefix = {}
    for i, s in enumerate(symbols[1:]):
        prefix[s] = 1 << (i+1)*10
    for symbol in reversed(symbols[1:]):
        if n >= prefix[symbol]:
            value = float(n) / prefix[symbol]
            return format % locals()
    return format % dict(symbol=symbols[0], value=n)


class Job:
    _kill_now = False

    def __init__(self, job_script, qsub_args='', out_file=None,
                 poll_interval=10):
        self.poll_interval = poll_interval
        self.job_script = job_script
        if qsub_args is not None:
            self.qsub_args = qsub_args
        else:
            self.qsub_args = ''
        self.out_file = out_file

        signal.signal(signal.SIGINT, self.__exit_gracefully)
        signal.signal(signal.SIGTERM, self.__exit_gracefully)

    def __repr__(self):
        rep = "job_script: " + self.job_script.rstrip() + "\n"
        rep += "qsub_args: " + self.qsub_args.rstrip() + "\n"
        if self.out_file is not None:
            rep += "out_file: " + self.out_file.rstrip() + "\n"
        return rep

    def run_job(self):
        """ qsub job
        """
        self.job_script_file = tempfile.NamedTemporaryFile(mode='w',
                                                           suffix='.sh',
                                                           delete=False)
        self.job_script_file.write(self.job_script)
        self.job_script_file.close()
        os.chmod(self.job_script_file.name, 0555)

        cmd = "qsub"
        cmd += " " + self.qsub_args
        cmd += " " + self.job_script_file.name
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        proc.wait()
        if proc.returncode != 0:
            raise Exception('unable to qsub job: {}'.format(cmd))
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
            time.sleep(random.randint(self.poll_interval, self.poll_interval + 20))
            if self._kill_now:
                sys.stderr.write("registering job for deletion\n")
                subprocess.Popen("qdel {}".format(self.job_id), shell=True)
                self.update_qstat()
                break
        if 'exit_status' in self.qstat:
            return self.qstat['exit_status']
        else:
            return 99

    def hit_mem_limit(self):
        """ True if mem limit was hit
        """
        try:
            mem_used = human2bytes(self.qstat['resources_used.mem'])
            mem_alloc = human2bytes(self.qstat['Resource_List.mem'])
            return mem_used + human2bytes('100M') > mem_alloc
        except:
            return False

    def hit_walltime_limit(self):
        """ True if walltime limit hit
        """
        try:
            walltime_used = time2secs(self.qstat['resources_used.walltime'])
            walltime_alloc = time2secs(self.qstat['Resource_List.walltime'])
            return walltime_used + 5 > walltime_alloc
        except:
            return False

    def restart(self, resource_multiplier):
        """ restart with higher reqs if req'd
        """
        mem = human2bytes(self.qstat['Resource_List.mem'])
        walltime = self.qstat['Resource_List.walltime']
        if self.hit_walltime_limit():
            walltime_secs = time2secs(walltime) * resource_multiplier
            walltime = secs2time(walltime_secs)
            sys.stderr.write("increased walltime to {}\n".format(walltime))
        if self.hit_mem_limit():
            mem = int(mem * resource_multiplier)
            sys.stderr.write("increased mem to {}\n".format(bytes2human(mem)))

        if self.qsub_args is not None and self.qsub_args.find('-l mem') != -1:
            self.qsub_args = re.sub(r'-l mem=(\S+)', '-l mem={}'.format(mem), self.qsub_args)
        else:
            self.qsub_args += ' -l mem={}'.format(mem)

        if self.qsub_args is not None and self.qsub_args.find('-l walltime') != -1:
            self.qsub_args = re.sub(r'-l walltime=(\S+)', '-l walltime={}'.format(walltime), self.qsub_args)
        else:
            self.qsub_args += ' -l walltime={}'.format(walltime)

        self.run_job()

    def _local_check_file(self, max_connection_retry):
        for attempt in range(max_connection_retry):
            try:
                local_file_size = int(subprocess.check_output('stat -c%s "{}"'.format(os.path.abspath(self.out_file)),
                                                            shell=True, stderr=DEVNULL))
                return local_file_size
            except:
                sys.stderr.write("Unable to stat {} locally\n".format(os.path.abspath(self.out_file)))
                time.sleep(10)
            else:
                sys.stderr.write("max retries for local file size check\n")
                return None
            if self._kill_now:
                return None

    def _remote_check_file(self, local_file_size, server, max_connection_retry):
        for attempt in range(max_connection_retry):
            try:
                remote_file_size = int(subprocess.\
                                       check_output('ssh {} stat -c%s "{}"'.
                                                    format(server, os.path.abspath(self.out_file)),
                                                    shell=True, stderr=DEVNULL))
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
            if self._kill_now:
                return False
        return True

    def check_file(self, servers, max_connection_retry=10):
        """ checks the file size of file f on remote servers and returns True
        if it matches local file size and also local file size > 0
        """
        time.sleep(10)
        if self.out_file is None:
            return True
        local_file_size = self._local_check_file(max_connection_retry)
        if local_file_size is None:
            return False
        for server in servers:
            if not self._remote_check_file(local_file_size, server,
                                           max_connection_retry):
                return False
        return local_file_size > 0

    def __exit_gracefully(self, signum, frame):
        sys.stderr.write("received interrupt\n")
        self._kill_now = True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='qsub_pbs.py',
                                     description='Submit job to PBS cluster.')
    parser.add_argument('-o', '--out_file',
                        nargs='?', help='file to check for non-zero size')
    parser.add_argument('-i', '--in_file',
                        nargs='?', default=sys.stdin, help='input script file')
    parser.add_argument('--poll_interval', nargs='?',
                        help='qstat polling interval', default=10)
    parser.add_argument('-c', '--check', default=False,
                        action='store_true',
                        help='check this file for non-zero size')
    parser.add_argument('-s', '--servers', nargs='*',
                        default=['gg02', 'gg05', 'gg06'],
                        help='check these servers for '
                        'non-zero output file size')
    parser.add_argument('--resource_multiplier', default=1.3, type=float,
                        help='job restart resource multiplier (if limit reached)')
    parser.add_argument('--max_restarts', default=10, type=int,
                        help='maximum number of job restarts')
    parser.add_argument('--args', nargs=argparse.REMAINDER)
    args = parser.parse_args()

    if args.args is not None:
        qsub_args = " ".join(args.args)
    else:
        qsub_args = None

    job = Job(args.in_file.read(), qsub_args, args.out_file,
              args.poll_interval)
    job.run_job()
    exit_status = job.wait()
    if args.check and exit_status == 0:
        if not job.check_file(args.servers):
            exit_status += 66

    i = 0
    while i < args.max_restarts and exit_status != 0 and not job._kill_now and \
            (job.hit_walltime_limit() or job.hit_mem_limit()):
        sys.stderr.write("job hit resource limit\n")
        sys.stderr.write("restarting job\n")
        job.restart(args.resource_multiplier)
        sys.stderr.write(str(job))
        exit_status = job.wait()
        i += 1

    sys.exit(exit_status)
