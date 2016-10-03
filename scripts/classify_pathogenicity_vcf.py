#!/usr/bin/env python

""" classify pathogenicity of vcf records, querying provean as necessary
"""

import os
import argparse
import vcf
import urllib2
from mechanize import Browser
from bs4 import BeautifulSoup
import pandas as pd
import re
import time
import sys
import requests
import tempfile
from subprocess import Popen
if 'SGE_ROOT' in os.environ:
    import drmaa
    import qsub
import qsub_pbs
import collections


class ProveanQuery:
    def __init__(self, ensembl_id, hgvsps, provean_script, num_provean_threads,
                 rest_server='http://grch37.rest.ensembl.org'):
        self.ensembl_id = ensembl_id
        self.fasta = self._get_fasta(ensembl_id, rest_server=rest_server)

        self.hgvsps = hgvsps
        self.provean_script = provean_script
        self.num_threads = num_provean_threads

        self._fasta_file = tempfile.NamedTemporaryFile(mode='w', delete=False)
        self._fasta_file.write(self.fasta)
        self._fasta_file.close()

        self._var_file = tempfile.NamedTemporaryFile(mode='w', delete=False)
        for var in hgvsps:
            self._var_file.write(var + "\n")
        self._var_file.close()
        self._out_file = tempfile.mktemp()
        self._cmd = '{} --num_threads {} -q {} -v {} |' \
            'sed "1,/PROVEAN scores/d" > {}'.format(provean_script, self.num_threads,
                                                    self._fasta_file.name, self._var_file.name,
                                                    self._out_file)

    def _get_fasta(self, ensembl_id, max_retry=30, rest_server='http://grch37.rest.ensembl.org'):
        """ get fasta sequence from ensembl
        """
        ext = '/sequence/id/{}?'.format(ensembl_id)
        for attempt in range(max_retry):
            try:
                sys.stderr.write("querying {}\n".format(rest_server + ext))
                r = requests.get(rest_server + ext, headers={"Content-Type": "text/x-fasta"})
                if not r.ok:
                    r.raise_for_status()
                return r.text
            except:
                sys.stderr.write("attempt {} failed...\n".format(attempt))
                time.sleep(10)
            else:
                sys.stderr.write('max attempts\n')
                return None

    def run_cluster(self, session=None, cluster_mode='SGE', mem_per_thread='1.5G'):
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
            sys.stderr.write('non-zero exit code; provean query failed: {} {}\n'.format(self.ensembl_id,
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
                sys.stderr.write('max attempts\n')
                sys.stderr.write('provean query failed: {} {}\n'.format(self.ensembl_id, ','.join(self.hgvsps)))
                return None


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


def add_provean_info(records, max_retry=3, rest_server='http://grch37.rest.ensembl.org',
                     provean_script='provean.sh', cluster_method='sge',
                     mem_per_thread='1.5G', num_provean_threads=4):
    """ add provean results using remote server or locally if necessary
    """
    query = ""
    for record in records:
        query += "{},{},{},{}\n".format(record.CHROM, record.POS, record.REF, record.ALT[0])

    # get the job URL
    job_url = None
    for attempt in range(max_retry):
        try:
            br = Browser()
            sys.stderr.write("Querying Provean:\n{}".format(query))
            br.open('http://provean.jcvi.org/genome_submit_2.php?species=human')
            br.form = list(br.forms())[1]  # select the chrpos form
            control = br.form.find_control("CHR")
            control.value = query
            br.submit()
            job_url = br.geturl()
            sys.stderr.write("job url: {}\n".format(job_url))
            if 'genome_prg_2' in job_url:
                job_url = None
                break
            if 'jobid' not in job_url:
                raise Exception("jobid not in job url")
            break
        except:
            sys.stderr.write("query attempt {} failed...\n".format(attempt))
            time.sleep(10)
        else:
            sys.stderr.write('max query attempts\n')
            break
    # parse job result page
    if job_url is not None:
        for attempt in range(max_retry):
            try:
                page = urllib2.urlopen(job_url).read()
                soup = BeautifulSoup(page, 'html.parser')
                link = soup.find('a', href=re.compile('one\.tsv'))
                url = 'http://provean.jcvi.org/' + link.get('href')
                df = pd.read_table(url)
                break
            except:
                sys.stderr.write("attempt {} failed...\n".format(attempt))
                time.sleep(10)
            else:
                sys.stderr.write('max attempts\n')
                return None
        for idx, record in enumerate(records):
            if df.ix[idx, 'PROTEIN_ID'] != 'record not found':
                record.INFO['provean_protein_id'] = df.ix[idx, 'PROTEIN_ID']
                record.INFO['provean_pred'] = df.ix[idx, 'PREDICTION (cutoff=-2.5)']
                record.INFO['provean_score'] = df.ix[idx, 'SCORE']
            else:
                record.INFO['provean_protein_id'] = '.'
                record.INFO['provean_pred'] = '.'
                record.INFO['provean_score'] = '.'
    else:
        add_provean_info_local(records, rest_server=rest_server, provean_script=provean_script,
                               cluster_mode=cluster_method, mem_per_thread=mem_per_thread,
                               num_provean_threads=num_provean_threads)

    return records


def get_hgvsp(record, max_retry=30, rest_server='http://grch37.rest.ensembl.org'):
    """ get hgvsp from ensembl
    """
    ext = '/vep/human/region/{}:{}-{}/{}/?hgvs=1'.format(record.CHROM, record.POS,
                                                         record.POS + len(record.REF) - 1, record.ALT[0])
    for attempt in range(max_retry):
        try:
            sys.stderr.write("querying {}\n".format(rest_server + ext))
            r = requests.get(rest_server + ext, headers={"Content-Type": "application/json"})
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
            sys.stderr.write('max attempts\n')
            return None


def add_provean_info_local(records, rest_server='http://grch37.rest.ensembl.org',
                           provean_script='provean.sh', cluster_mode='SGE',
                           num_provean_threads=4, mem_per_thread='1.5G', max_retry=30):
    """ run provean locally
    """
    session = None
    if cluster_mode.lower() == 'sge':
        session = drmaa.Session()
        session.initialize()
    provean_queries = collections.defaultdict(list)
    for record in records:
        hgvsp_ensps = get_hgvsp(record, rest_server=rest_server)
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
                                                        num_provean_threads=num_provean_threads,
                                                        provean_script=provean_script))
    jobs = []
    for record, queries in provean_queries.iteritems():
        for query in queries:
            if cluster_mode.lower() == 'none':
                query.run_local()
            else:
                jobs.append(query.run_cluster(session=session, cluster_mode=cluster_mode,
                                              mem_per_thread=mem_per_thread))
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


def is_fs_splice_stop(record):
    ann_effect = get_ann_effect(record)
    return any([c in ef for ef in ann_effect for c in ["frameshift", "splice_donor", "splice_acceptor", "stop_gained"]])


def is_missense(record):
    ann_effect = get_ann_effect(record)
    return any(["missense_variant" in ef for ef in ann_effect])


def get_fs_splice_stop_pathogenicity(record):
    if (is_loh(record) or is_hap_insuf(record)) and is_cancer_gene(record):
        return "pathogenic"
    elif is_loh(record) or is_hap_insuf(record) or is_cancer_gene(record):
        return "potentially_pathogenic"
    else:
        return "passenger"


def is_mt_pathogenic(record):
    pathogenic = False
    if 'dbNSFP_MutationTaster_pred' in record.INFO:
        pathogenic = 'D' in record.INFO['dbNSFP_MutationTaster_pred'] or \
            'A' in record.INFO['dbNSFP_MutationTaster_pred']
    elif 'MutationTaster_pred' in record.INFO and not all(x is None for x in record.INFO['MutationTaster_pred']):
        pathogenic = any(['disease' in info for info in record.INFO['MutationTaster_pred']])
    return pathogenic


def classify_missense_pathogenicity(record):
    if not is_mt_pathogenic(record) and not is_chasm_pathogenic(record):
        return "passenger"
    else:
        if is_fathmm_pathogenic(record) or is_chasm_pathogenic(record):
            return "pathogenic" if is_cancer_gene(record) else "potentially_pathogenic"
        else:
            return "passenger"


def is_chasm_pathogenic(record):
    cp = filter(lambda x: x.endswith('chasm_score'), record.INFO.keys())
    chasm_scores = [min(record.INFO[x]) for x in cp]
    return any([x <= 0.3 for x in chasm_scores])


def is_fathmm_pathogenic(record):
    return "CANCER" in record.INFO['fathmm_pred'] if 'fathmm_pred' in record.INFO else False


def is_cancer_gene(record):
    return 'lawrence' in record.INFO or 'kandoth' in record.INFO or 'cancer_gene_census' in record.INFO


def is_hap_insuf(record):
    return 'hap_insuf' in record.INFO


def get_ann_effect(record):
    return [x.split('|')[1] for x in record.INFO['ANN']]


def is_loh(record):
    return record.INFO["facetsLCN_EM"] == 0 if 'facetsLCN_EM' in record.INFO else False


def get_missense_pathogenicity(record):
    if is_mt_pathogenic(record) or is_chasm_pathogenic(record):
        if is_fathmm_pathogenic(record) or is_chasm_pathogenic(record):
            return "pathogenic" if is_cancer_gene(record) else "potentially_pathogenic"
        else:
            return "passenger"
    else:
        return "passenger"


def is_provean_pathogenic(record):
    return 'provean_pred' in record.INFO and record.INFO['provean_pred'] == 'Deleterious'


def get_provean_pathogenicity(record):
    if not is_mt_pathogenic(record) and not is_provean_pathogenic(record):
        return "passenger"
    else:
        if (is_loh(record) or is_hap_insuf(record)) and is_cancer_gene(record):
            return "pathogenic"
        elif is_loh(record) or is_hap_insuf(record) or is_cancer_gene(record):
            return "potentially_pathogenic"
        else:
            return "passenger"


def classify_pathogenicity(record):
    """ classify pathogenicity for a vcf record
    """
    if is_fs_splice_stop(record):
        record.INFO["pathogenicity"] = get_fs_splice_stop_pathogenicity(record)
    elif is_missense(record):
        record.INFO["pathogenicity"] = get_missense_pathogenicity(record)
    elif is_inframe(record):
        record.INFO["pathogenicity"] = get_provean_pathogenicity(record)
    else:
        record.INFO["pathogenicity"] = None


def is_inframe(record):
    ann_effect = get_ann_effect(record)
    return any(["inframe" in ef for ef in ann_effect])


def should_run_provean(record):
    return not is_fs_splice_stop(record) and not is_missense(record) and is_inframe(record) and not is_mt_pathogenic(record)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='classify_pathogenicity_vcf.py',
                                     description='Add pathogenicity to vcf file')
    parser.add_argument('vcf_infile')
    parser.add_argument('--mem_per_thread', nargs='?', default='1.5G', help='memory per provean thread')
    parser.add_argument('--provean_script', nargs='?', default='provean.sh', help='provean script')
    parser.add_argument('--cluster_mode', nargs='?', default='SGE', help='cluster mode')
    parser.add_argument('--qsub_queue', nargs='?', default='jrf.q,all.q', help='qsub queue')
    parser.add_argument('--num_provean_threads', nargs='?', default=4, type=int, help='number of provean threads')
    parser.add_argument('--run_local', action='store_true', default=False, help='run provean locally')
    args = parser.parse_args()

    vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))

    assert "hap_insuf" in vcf_reader.infos
    assert "ANN" in vcf_reader.infos
    assert "kandoth" in vcf_reader.infos
    assert "lawrence" in vcf_reader.infos
    assert "facetsLCN_EM" in vcf_reader.infos

    # add necessary info headers
    vcf_reader.infos['pathogenicity'] = vcf.parser._Info(id='pathogenicity', num=-1, type='String',
                                                         desc="Classification of pathogenicity")
    vcf_reader.infos['provean_protein_id'] = vcf.parser._Info(id='provean_protein_id', num=-1, type='String',
                                                              desc="provean protein id (run if necessary)")
    vcf_reader.infos['provean_pred'] = vcf.parser._Info(id='provean_pred', num=-1, type='String',
                                                        desc="provean prediction (run if necessary)")
    vcf_reader.infos['provean_score'] = vcf.parser._Info(id='provean_score', num=-1, type='Float',
                                                         desc="provean score (run if necessary)")

    vcf_writer = vcf.Writer(sys.stdout, vcf_reader)

    records = list()
    provean_records = list()
    for record in vcf_reader:
        if should_run_provean(record):
            provean_records.append(record)
        records.append(record)
    if len(provean_records) > 0:
        if args.run_local:
            add_provean_info_local(provean_records, mem_per_thread=args.mem_per_thread,
                                   provean_script=args.provean_script,
                                   cluster_mode=args.cluster_mode,
                                   num_provean_threads=args.num_provean_threads)
        else:
            add_provean_info(provean_records, mem_per_thread=args.mem_per_thread,
                             provean_script=args.provean_script,
                             num_provean_threads=args.num_provean_threads)

    for record in records:
        classify_pathogenicity(record)
        vcf_writer.write_record(record)
    vcf_writer.close()
