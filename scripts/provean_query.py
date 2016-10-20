import os
import drmaa
import time
import collections
import qsub_pbs
import requests
import re
import tempfile
import sys
import pandas as pd
from subprocess import Popen


def three_to_one_amino_acid_code(x):
    # 3 > 1 AA code mapping
    d = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
         'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N',
         'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W',
         'Ala': 'A', 'Val': 'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M'}
    # inverse to get 1 > 3 AA code mapping
    for aa in d:
        x = re.sub(aa, d[aa], x)
    return x


class ProveanQuery:
    def __init__(self, ensembl_id, hgvsps, num_provean_threads, provean_script,
                 rest_server='http://grch37.rest.ensembl.org'):
        self.rest_server = rest_server
        self.ensembl_id = ensembl_id
        self.hgvsps = hgvsps
        self.num_threads = num_provean_threads
        self.provean_script = provean_script
        self.fasta = self._get_fasta()
        self.fasta_file = tempfile.NamedTemporaryFile(mode='w', delete=False)
        self.fasta_file.write(self.fasta)
        self.fasta_file.close()

        self.var_file = tempfile.NamedTemporaryFile(mode='w', delete=False)
        for var in hgvsps:
            self.var_file.write(var + "\n")
        self.var_file.close()
        self.out_file = tempfile.mktemp()
        self._cmd = '{} --num_threads {} -q {} -v {} |' \
            'sed "1,/PROVEAN scores/d" > {}'.format(provean_script, num_provean_threads,
                                                    self.fasta_file.name, self.var_file.name,
                                                    self.out_file)

    def _get_fasta(self, max_retry=30):
        """ get fasta sequence from ensembl
        """
        ext = '/sequence/id/{}?'.format(self.ensembl_id)
        for attempt in range(max_retry):
            try:
                sys.stderr.write("querying {}\n".format(self.rest_server + ext))
                r = requests.get(self.rest_server + ext, headers={"Content-Type": "text/x-fasta"})
                if not r.ok:
                    r.raise_for_status()
                return r.text
            except:
                sys.stderr.write("attempt {} failed...\n".format(attempt))
                time.sleep(10)
            else:
                raise Exception('max attempts to retrieve fasta')

    def run_cluster(self, session=None, cluster_mode='SGE', mem_per_thread='1.5G'):
        import qsub
        sys.stderr.write('running provean on the cluster...\n')
        job = None
        if cluster_mode.lower() == 'pbs':
            mem = qsub_pbs.human2bytes(mem_per_thread) * self.num_threads
            job = qsub_pbs.Job(self._cmd, '-l nodes=1:ppn={} -l walltime=12:00:00 '
                               '-l mem={}'.format(self.num_threads, mem))
            job.run_job()
        elif cluster_mode.lower() == 'sge':
            assert session is not None
            job = qsub.Job(session=session, job_script=self._cmd,
                           qsub_args='-pe smp {} -l h_vmem={}'.format(self.num_threads, mem_per_thread))
        else:
            raise Exception("Invalid cluster mode\n")
        return job

    def run_local(self):
        sys.stderr.write('running provean locally...\n')
        process = Popen(self._cmd, shell=True)
        process.wait()
        if process.returncode != 0:
            raise Exception('non-zero exit code; provean query failed: {} {}\n'.format(self.ensembl_id,
                                                                                       ",".join(self.hgvsps)))
            return None

    def process(self, max_retry=10):
        for attempt in range(max_retry):
            try:
                assert os.stat(self._out_file).st_size > 0
                df = pd.read_table(self._out_file)
                assert len(df) > 0
                df.columns = ['hgvsp', 'score']
                return min(df['score'])
            except:
                sys.stderr.write('file is 0\n')
                time.sleep(10)
            else:
                raise Exception('max attempts; provean query failed: {} {}\n'.format(
                    self.ensembl_id, ','.join(self.hgvsps)))


class ProveanQueryManager:
    def __init__(self, records, provean_script, num_provean_threads, mem_per_thread='2G',
                 rest_server='http://grch37.rest.ensembl.org', cluster_mode='none'):
        self.records = records
        self.provean_script = provean_script
        self.num_threads = num_provean_threads
        self.mem_per_thread = mem_per_thread
        self.rest_server = rest_server
        self.cluster_mode = cluster_mode

    def run_queries(self):
        provean_queries = self._create_queries()
        session = None
        if self.cluster_mode.lower() == 'sge':
            session = drmaa.Session()
            session.initialize()
        jobs = []
        for record, queries in provean_queries.iteritems():
            for query in queries:
                if self.cluster_mode.lower() == 'none':
                    query.run_local()
                else:
                    jobs.append(query.run_cluster(session=session, cluster_mode=self.cluster_mode,
                                                    mem_per_thread=self.mem_per_thread))
        for job in jobs:
            job.wait()
        for record, queries in provean_queries.iteritems():
            scores = []
            for query in queries:
                scores.append(query.process())
            record.INFO['provean_protein_id'] = query.ensembl_id
            record.INFO['provean_score'] = min(scores)
            if min(scores) < -2.5:
                record.INFO['provean_pred'] = 'Deleterious'
            else:
                record.INFO['provean_pred'] = 'Neutral'
        if session is not None:
            session.exit()

    def _create_queries(self):
        provean_queries = collections.defaultdict(list)
        for record in self.records:
            hgvsp_ensps = self._get_hgvsp(record)
            # get ensembl id -> HGVSp
            ensp_hgvsp = collections.defaultdict(list)
            for hgvsp_ensp in hgvsp_ensps:
                ensembl_id = hgvsp_ensp.split(":p.")[0]
                hgvsp = three_to_one_amino_acid_code(hgvsp_ensp.split(":p.")[1])
                # ignore frameshifts
                if 'fs' in hgvsp:
                    continue
                ensp_hgvsp[ensembl_id].append(hgvsp)
            for ensp, hgvsps in ensp_hgvsp.iteritems():
                provean_queries[record].append(ProveanQuery(ensembl_id=ensp, hgvsps=hgvsps,
                                                            num_provean_threads=self.num_threads,
                                                            provean_script=self.provean_script))
        return provean_queries

    def _get_hgvsp(self, record, max_retry=30):
        """ get hgvsp from ensembl
        """
        ext = '/vep/human/region/{}:{}-{}/{}/?hgvs=1'.format(record.CHROM, record.POS,
                                                             record.POS + len(record.REF) - 1, record.ALT[0])
        for attempt in range(max_retry):
            try:
                sys.stderr.write("querying {}\n".format(self.rest_server + ext))
                r = requests.get(self.rest_server + ext, headers={"Content-Type": "application/json"})
                if not r.ok:
                    r.raise_for_status()
                decoded = r.json()
                hgvsp_decoded = filter(lambda x: 'hgvsp' in x, decoded[0]['transcript_consequences'])
                hgvsp = map(lambda x: x['hgvsp'], hgvsp_decoded)
                return hgvsp
            except:
                sys.stderr.write("attempt {} failed...\n".format(attempt))
                time.sleep(10)
            else:
                raise Exception('max attempts to retrieve hgvsp from ensembl')
