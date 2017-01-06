#!/usr/bin/env python
""" submit jobs to a cluster with drmaa support and wait for them to finish
"""

import drmaa
import argparse
import os
import random
import tempfile
import sys
import signal
import subprocess
import time

DEVNULL = open(os.devnull, 'wb')


class Job:
    _kill_now = False

    def __init__(self, session, job_script, qsub_args, out_file=None,
                 poll_interval=60):
        self.session = session
        self.poll_interval = poll_interval
        self.job_script = job_script
        self.out_file = out_file
        self.job_script_file = tempfile.NamedTemporaryFile(mode='w',
                                                           suffix='.sh', delete=False)
        self.job_script_file.write(job_script)
        self.job_script_file.close()
        os.chmod(self.job_script_file.name, 0555)
        jt = session.createJobTemplate()
        jt.remoteCommand = self.job_script_file.name
        jt.workingDirectory = os.getcwd()
        jt.joinFiles = True
        jt.nativeSpecification = qsub_args
        self.id = session.runJob(jt)
        session.deleteJobTemplate(jt)

        signal.signal(signal.SIGINT, self.__exit_gracefully)
        signal.signal(signal.SIGTERM, self.__exit_gracefully)

    def wait(self):
        retval = None
        while True:
            try:
                retval = self.session.wait(self.id, self.poll_interval +
                                           random.randrange(0, 60))
                if retval.hasExited or retval.hasSignal \
                        or retval.wasAborted or retval.hasCoreDump:
                    break
            except drmaa.ExitTimeoutException:
                pass
            if self._kill_now:
                print "terminating job"
                self.session.control(self.id, drmaa.JobControlAction.TERMINATE)
                print "job terminated"
                break
        return retval

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
                    sys.stderr.write(str(e) + "\n")
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
        print "received interrupt"
        self._kill_now = True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='qsub.py',
                                     description='Submit job to cluster'
                                     ' using drmaa.'
                                     'Read the script from STDIN')
    parser.add_argument('-o', '--out_file',
                        nargs='?', help='file to check for non-zero size')
    parser.add_argument('-i', '--in_file',
                        nargs='?', default=sys.stdin, help='input script file')
    parser.add_argument('--poll_interval', nargs='?',
                        help='job polling interval', default=60)
    parser.add_argument('-c', '--check', default=False,
                        action='store_true',
                        help='check this file for non-zero size')
    parser.add_argument('-s', '--servers',
                        default=['gg02', 'gg05', 'gg06'], nargs='*',
                        help='check these servers for '
                        'non-zero output file size')
    parser.add_argument('--args', nargs=argparse.REMAINDER)
    args = parser.parse_args()

    s = drmaa.Session()
    s.initialize()
    qsub_args = " ".join(args.args)

    job = Job(session=s, job_script=args.in_file.read(), qsub_args=qsub_args, out_file=args.out_file,
              poll_interval=args.poll_interval)
    retval = job.wait()
    exit_status = 0
    if retval is not None:
        exit_status = retval.exitStatus
        if args.check and retval.hasExited and exit_status == 0:
            if not job.check_file(args.servers):
                exit_status += 66
    else:
        exit_status = 1
    s.exit()
    sys.exit(exit_status)
