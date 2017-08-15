""" DRMAA jobs
separate from other cluster jobs because of the required drmaa library
"""
import drmaa
import os
import atexit
import tempfile
import random
from job import ClusterJob


class DRMAAJob(ClusterJob):

    def __init__(self, job_script, qsub_args, out_file=None,
                 poll_interval=60, remote_check_servers=None):
        ClusterJob.__init__(self, out_file, remote_check_servers)
        self.session = drmaa.Session()
        self.session.initialize()
        self.poll_interval = poll_interval
        self.job_script = job_script
        self.qsub_args = qsub_args

        atexit.register(self._cleanup)

    def run_job(self):
        job_script_file = tempfile.NamedTemporaryFile(mode='w',
                                                      suffix='.sh', delete=False)
        job_script_file.write(self.job_script)
        job_script_file.close()
        os.chmod(job_script_file.name, 0o555)
        jt = self.session.createJobTemplate()
        jt.remoteCommand = job_script_file.name
        jt.workingDirectory = os.getcwd()
        jt.joinFiles = True
        jt.nativeSpecification = self.qsub_args
        self.job_id = self.session.runJob(jt)
        self.session.deleteJobTemplate(jt)

    def wait(self):
        self.retval = None
        while True:
            try:
                self.retval = self.session.wait(self.job_id, self.poll_interval +
                                                random.randrange(0, 60))
                if self.retval.hasExited or self.retval.hasSignal \
                        or self.retval.wasAborted or self.retval.hasCoreDump:
                    break
            except drmaa.ExitTimeoutException:
                pass
            if self._kill_now:
                print("terminating job")
                self.session.control(self.job_id, drmaa.JobControlAction.TERMINATE)
                print("job terminated")
                break
        return self.retval

    def is_finished(self):
        return self.retval is not None and self.retval.hasExited and self.retval.exitStatus == 0

    def _cleanup(self):
        self.session.exit()
