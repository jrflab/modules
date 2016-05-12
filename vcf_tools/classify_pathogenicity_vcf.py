#!/usr/bin/env python

""" classify pathogenicity of vcf records, querying provean as necessary
"""

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


def add_provean_info(records, max_retry=30, rest_server='http://grch37.rest.ensembl.org',
                     provean_script='provean.sh', qsub_script='perl modules/qsub.pl', qsub_queue='jrf.q,all.q',
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
            sys.stderr.write("Querying Provean: {}\n".format(query))
            br.open('http://provean.jcvi.org/genome_submit_2.php?species=human')
            br.form = list(br.forms())[1]  # select the chrpos form
            control = br.form.find_control("CHR")
            control.value = query
            br.submit()
            job_url = br.geturl()
            sys.stderr.write("job url: {}\n".format(job_url))
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
        for record in records:
            add_provean_info_local(record, rest_server=rest_server, provean_script=provean_script,
                                   qsub_script=qsub_script, qsub_queue=qsub_queue, mem_per_thread=mem_per_thread,
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


def get_fasta(ensembl_id, max_retry=30, rest_server='http://grch37.rest.ensembl.org'):
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


def add_provean_info_local(records, rest_server='http://grch37.rest.ensembl.org',
                           provean_script='provean.sh', qsub_script='perl modules/scripts/qsub.pl',
                           qsub_queue='jrf.q,all.q', num_provean_threads=4, mem_per_thread='1.5G'):
    """ run provean locally
    """
    provean_queries = []
    for record in records:
        provean_queries.extend(run_provean_query(record))
    process_provean_queries(provean_queries)


def run_provean_query(record, rest_server='http://grch37.rest.ensembl.org',
                      provean_script='provean.sh', qsub_script='perl modules/scripts/qsub.pl',
                      qsub_queue='jrf.q,all.q', num_provean_threads=4, mem_per_thread='1.5G'):
    sys.stderr.write('running provean locally...\n')
    hgvsp_ensps = get_hgvsp(record, rest_server=rest_server)
    ensp_hgvsp = {}
    for hgvsp_ensp in hgvsp_ensps:
        ensembl_id = hgvsp_ensp.split(":p.")[0]
        hgvsp = three_to_one_amino_acid_code(hgvsp_ensp.split(":p.")[1])
        # ignore frameshifts
        if 'fs' in hgvsp:
            continue
        if ensembl_id not in ensp_hgvsp:
            ensp_hgvsp[ensembl_id] = []
        ensp_hgvsp[ensembl_id].append(hgvsp)
    provean_queries = []
    for ensp in ensp_hgvsp:
        fasta = get_fasta(ensp, rest_server=rest_server)
        tmp_fasta = tempfile.NamedTemporaryFile(mode='w', delete=False)
        tmp_fasta.write(fasta)
        tmp_var = tempfile.NamedTemporaryFile(mode='w', delete=False)
        for var in ensp_hgvsp[ensp]:
            tmp_var.write(var + "\n")
        tmp_fasta.close()
        tmp_var.close()
        tmp_out = tempfile.mktemp()
        cmd = 'echo "{} --num_threads {} -q {} -v {} |' \
            'sed \\"1,/PROVEAN scores/d\\" > {}"'.format(provean_script, num_provean_threads,
                                                         tmp_fasta.name, tmp_var.name, tmp_out)
        if qsub_script is not None:
            cmd += '| {} -- -j y -q {} -V -b n -o err.test -pe smp {} -l h_vmem={}'.format(qsub_script, qsub_queue,
                                                                                           num_provean_threads,
                                                                                           mem_per_thread)
        else:
            cmd += "| bash"
        sys.stderr.write('running: {}\n'.format(cmd))
        process = Popen(cmd, shell=True)
        query = {'record': record,
                 'out': tmp_out,
                 'process': process,
                 'ensp': ensp,
                 'hgvsp': ensp_hgvsp[ensp],
                 'files': [tmp_fasta, tmp_var, tmp_out]}
        provean_queries.append(query)
    return provean_queries


def process_provean_queries(provean_queries):
    for query in provean_queries:
        query['process'].wait()
        if query['process'].returncode == 0:
            df = pd.read_table(query['out'], comment='#')
            df.columns = ['hgvsp', 'score']
            if 'provean_protein_id' not in record.INFO:
                record.INFO['provean_protein_id'] = []
                record.INFO['provean_score'] = []
            record.INFO['provean_protein_id'].append(query['ensp'])
            record.INFO['provean_score'].append(df['score'])
        else:
            sys.stderr.write('provean query failed: {} {}\n'.format(query['ensp'], query['hgvsp']))
    if 'provean_score' in record.INFO:
        record.INFO['provean_pred'] = ['Deleterious' if x < -2.5 else 'Neutral' for x in record.INFO['provean_score']]


def is_fs_splice_stop(record):
    ann_effect = get_ann_effect(record)
    return any([c in ef for ef in ann_effect for c in ["frameshift", "splice_donor", "splice_acceptor", "stop_gained"]])


def is_missense(record):
    ann_effect = get_ann_effect(record)
    return any(["missense_variant" in ef for ef in ann_effect])


def get_fs_splice_stop_pathogenicity(record):
    if (is_loh(record) or is_hap_insuf(record)) and is_cancer_gene(record):
        record.INFO["pathogenicity"] = "pathogenic"
    elif is_loh(record) or is_hap_insuf(record) or is_cancer_gene(record):
        record.INFO["pathogenicity"] = "potentially_pathogenic"
    else:
        record.INFO["pathogenicity"] = "passenger"


def is_mt_pathogenic(record):
    pathogenic = False
    if 'MutationTaster_pred' in record.INFO:
        pathogenic = 'disease' in record.INFO['MutationTaster_pred']
    elif 'dbNSFP_MutationTaster_pred' in record.INFO:
        pathogenic = record.INFO['dbNSFP_MutationTaster_pred'] == 'D' or \
            record.INFO['dbNSFP_MutationTaster_pred'] == 'A'
    return pathogenic


def classify_missense_pathogenicity(record):
    if ~is_mt_pathogenic(record) and ~is_chasm_pathogenic(record):
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
    return record.INFO['fathmm_pred'] == "CANCER" if 'fathmm_pred' in record.INFO else False


def is_cancer_gene(record):
    return 'lawrence' in record.INFO or 'kandoth' in record.INFO


def is_hap_insuf(record):
    return 'hap_insuf' in record.INFO


def get_ann_effect(record):
    return [x.split('|')[1] for x in record.INFO['ANN']]


def is_loh(record):
    return record.INFO["facetsLCN_EM"] == 0 if 'facetsLCN_EM' in record.INFO else False


def get_missense_pathogenicity(record):
    if ~is_mt_pathogenic(record) and ~is_chasm_pathogenic(record):
        return "passenger"
    else:
        if is_fathmm_pathogenic(record) or is_chasm_pathogenic(record):
            return "pathogenic" if is_cancer_gene else "potentially_pathogenic"
        else:
            return "passenger"


def is_provean_pathogenic(record):
    return 'provean_pred' in record.INFO and any([x == 'Deleterious' for x in record.INFO['provean_pred']])


def get_provean_pathogenicity(record):
    if ~is_mt_pathogenic(record) and ~is_provean_pathogenic(record):
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
    elif is_provean_record(record):
        record.INFO["pathogenicity"] = get_provean_pathogenicity(record)


def is_inframe(record):
    ann_effect = get_ann_effect(record)
    return any(["inframe" in ef for ef in ann_effect])


def is_provean_record(record):
    return ~is_fs_splice_stop(record) and ~is_missense(record) and is_inframe(record) and ~is_mt_pathogenic(record)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='classify_pathogenicity_vcf.py',
                                     description='Add pathogenicity to vcf file')
    parser.add_argument('vcf_infile')
    parser.add_argument('--mem_per_thread', nargs='?', default='1.5G', help='memory per provean thread')
    parser.add_argument('--provean_script', nargs='?', default='provean.sh', help='provean script')
    parser.add_argument('--qsub_script', nargs='?', default='perl modules/scripts/qsub.pl', help='qsub script')
    parser.add_argument('--qsub_queue', nargs='?', default='jrf.q,all.q', help='qsub queue')
    parser.add_argument('--num_provean_threads', nargs='?', default=4, type=int, help='number of provean threads')
    args = parser.parse_args()

    vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))

    assert "hap_insuf" in vcf_reader.infos
    assert "ANN" in vcf_reader.infos
    assert "kandoth" in vcf_reader.infos
    assert "lawrence" in vcf_reader.infos
    assert "facetsLCN_EM" in vcf_reader.infos
    assert any(["MutationTaster_pred" in x for x in vcf_reader.infos])

    # add necessary info headers
    vcf_reader.infos['pathogenicity'] = vcf.parser._Info(id='pathogenicity', num=-1, type='String',
                                                         desc="Classification of pathogenicity",
                                                         source=None, version=None)
    vcf_reader.infos['provean_protein_id'] = vcf.parser._Info(id='provean_protein_id', num=-1, type='String',
                                                              desc="provean protein id (run if necessary)",
                                                              source=None, version=None)
    vcf_reader.infos['provean_pred'] = vcf.parser._Info(id='provean_pred', num=-1, type='String',
                                                        desc="provean prediction (run if necessary)",
                                                        source=None, version=None)
    vcf_reader.infos['provean_score'] = vcf.parser._Info(id='provean_score', num=-1, type='Float',
                                                         desc="provean score (run if necessary)",
                                                         source=None, version=None)

    vcf_writer = vcf.Writer(sys.stdout, vcf_reader)

    records = list()
    provean_records = list()
    for record in vcf_reader:
        if is_provean_record(record):
            provean_records.append(record)
        records.append(record)
    add_provean_info(provean_records, mem_per_thread=args.mem_per_thread, provean_script=args.provean_script,
                     qsub_script=args.qsub_script, num_provean_threads=args.num_provean_threads)

    for record in records:
        classify_pathogenicity(record)
        vcf_writer.write_record(record)
    vcf_writer.close()
