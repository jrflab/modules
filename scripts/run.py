#!/usr/bin/env python
""" runs a job either locally or via the cluster
"""

import argparse
import sys
import errno
import os
import job
import math
import re


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='run.py',
                                     description='run jobs')
    parser.add_argument('-i', '--input', default=sys.stdin, type=argparse.FileType('r'), help='job script')
    parser.add_argument('-v', '--env', default=None, help='anaconda env')
    parser.add_argument('--default_env', default=None, help='default anaconda env')
    parser.add_argument('-s', '--soft_memory', default='1G', help='soft memory limit per core')
    parser.add_argument('-m', '--hard_memory', default='1G', help='hard memory limit per core')
    parser.add_argument('-n', '--num_cores', default=1, type=int, help='number of cores')
    parser.add_argument('--sge_parallel_env', default='smp', help='SGE parallel env')
    parser.add_argument('-w', '--walltime', default='24:00:00', help='wall time')
    parser.add_argument('-c', '--check', default=False, action='store_true', help='check for non-zero file size')
    parser.add_argument('-d', '--docker', default=False, action='store_true', help='request docker support')
    parser.add_argument('-I', '--internet', default=False, action='store_true', help='request internet access')
    parser.add_argument('-g', '--cluster_engine', default='sge', help='cluster engine (sge, lfs, or pbs supported)')
    parser.add_argument('-l', '--local', default=False, action='store_true', help='run job locally')
    parser.add_argument('-o', '--out_file', default=None, help='output file to check')
    parser.add_argument('-p', '--project_name', default=None, help='project name')
    parser.add_argument('-N', '--job_name', default=None, help='job name')
    parser.add_argument('-e', '--log_file', default=None, help='log file')
    parser.add_argument('-S', '--shell', default='/bin/bash', help='shell')
    parser.add_argument('--servers', default=None, nargs='*',
                        help='use these servers for checking non-zero output file size')
    args = parser.parse_args()

    cluster_engine = args.cluster_engine.lower()
    total_hard_mem = args.num_cores * job.human2bytes(args.hard_memory)
    total_soft_mem = args.num_cores * job.human2bytes(args.soft_memory)

    env = None
    if args.env is not None:
        env = args.env
    elif args.default_env is not None:
        env = args.default_env

    if args.log_file is not None:
        log_dir = os.path.dirname(args.log_file)
        if not os.path.exists(log_dir):
            try:
                os.makedirs(log_dir)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise

    if args.out_file is not None:
        out_dir = os.path.dirname(args.out_file)
        if not os.path.exists(out_dir):
            try:
                os.makedirs(out_dir)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise

    job_script = args.input.read()
    if env is not None:
        job_script = "tries=0; " \
            "until [[ $tries -gt 10 ]] || source {env}/bin/activate {env}; do " \
            "tries=`expr $tries + 1`; jitter=`expr $RANDOM % 10 + 1`; " \
            "sleep `expr $jitter \* $tries`; done; {script}".format(env=env, script=job_script)

    job_script = "set -o pipefail; umask 002; {script}".format(script=job_script)
    if args.shell is not None:
        job_script = "#!{shell}\n{script}".format(shell=args.shell, script=job_script)

    job_name = None
    if args.job_name is not None and args.project_name is not None:
        job_name = "{}_{}".format(args.project_name, args.job_name)
    elif args.job_name is not None:
        job_name = args.job_name
    elif args.project_name is not None and args.out_file is not None:
        job_name = "{}_{}".format(args.project_name, os.path.basename(args.out_file))

    if args.local or (args.internet and args.cluster_engine != 'lsf'):
        my_job = job.LocalJob(job_script=job_script, out_file=args.out_file, log_file=args.log_file, shell=args.shell)
    elif cluster_engine == 'sge':
        qsub_args = "-V -wd {pwd} -now n -notify -b n -S {shell}".format(pwd=os.getcwd(), shell=args.shell)
        if job_name is not None:
            qsub_args += " -N {}".format(job_name)
        if args.log_file is not None:
            qsub_args += " -j y -o {}".format(args.log_file)
        if args.num_cores > 1:
            qsub_args += " -pe {penv} {num_cores}".format(penv=args.sge_parallel_env,
                                                          num_cores=args.num_cores)
        qsub_args += " -l virtual_free={soft_mem},h_vmem={hard_mem}".format(soft_mem=args.soft_memory,
                                                                            hard_mem=args.hard_memory)
        import drmaa_job
        my_job = drmaa_job.DRMAAJob(job_script=job_script, qsub_args=qsub_args, out_file=args.out_file,
                                    remote_check_servers=args.servers)
    elif cluster_engine == 'pbs':
        qsub_args = "-d {pwd} -S {shell}".format(pwd=os.getcwd(), shell=args.shell)
        if job_name is not None:
            qsub_args += " -N {}".format(job_name)
        if args.log_file is not None:
            qsub_args += " -j oe -o {}".format(args.log_file)
        qsub_args += " -l nodes=1:ppn={num_cores}".format(num_cores=args.num_cores)
        if args.docker:
            qsub_args += ":docker"
        qsub_args += " -l walltime={walltime} -l mem={total_mem}".format(walltime=args.walltime,
                                                                         total_mem=total_hard_mem)
        my_job = job.PBSJob(job_script=job_script, qsub_args=qsub_args, out_file=args.out_file,
                            remote_check_servers=args.servers)
    elif cluster_engine == 'lsf':
        qsub_args = "-cwd {pwd} -L {shell}".format(pwd=os.getcwd(), shell=args.shell)
        if job_name is not None:
            qsub_args += " -J {}".format(job_name)
        if args.log_file is not None:
            qsub_args += " -o {}".format(args.log_file)
        if args.internet:
            qsub_args += ' -R "select[internet]"'
        if args.num_cores > 1:
            qsub_args += ' -n {num_cores} -R "span[hosts=1]"'.format(num_cores=args.num_cores)
        hard_mem_gb = int(math.ceil(job.human2bytes(args.hard_memory) / 1000000000.0))
        soft_mem_gb = int(math.ceil(job.human2bytes(args.soft_memory) / 1000000000.0))
        walltime = re.sub(r'(\d+):(\d+):(\d+)', r'\g<1>:\g<2>', args.walltime)
        qsub_args += ' -M {hard_mem_gb} -R "rusage[mem={soft_mem_gb}]" -W {walltime}'.format(
            hard_mem_gb=hard_mem_gb, soft_mem_gb=soft_mem_gb, walltime=walltime)
        my_job = job.LSFJob(job_script=job_script, qsub_args=qsub_args, out_file=args.out_file,
                            remote_check_servers=None)  # checking filesize on remote servers not supported
    else:
        raise ValueError("invalid cluster engine. supported engines are sge, pbs, and lsf")

    my_job.run_job()
    my_job.wait()

    exit_status = 0
    if my_job.is_finished():
        if args.check and not my_job.check_file():
            exit_status += 66
    else:
        exit_status = 1
    sys.exit(exit_status)
