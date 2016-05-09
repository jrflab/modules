#!/usr/bin/env python

import argparse
import vcf
import urllib2
from mechanize import Browser
from bs4 import BeautifulSoup
import pandas as pd
import re
import time
import sys

parser = argparse.ArgumentParser(prog='classify_pathogenicity_vcf.py',
                                 description='Add pathogenicity to vcf file')
parser.add_argument('vcf_infile')

args = parser.parse_args()

vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))

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

assert "hap_insuf" in vcf_reader.infos
assert "ANN" in vcf_reader.infos
assert "kandoth" in vcf_reader.infos
assert "lawrence" in vcf_reader.infos
assert "facetsLCN_EM" in vcf_reader.infos
assert any(["MutationTaster_pred" in x for x in vcf_reader.infos])

vcf_writer = vcf.Writer(sys.stdout, vcf_reader)


def query_provean(record, max_retry):
    query = str(record.CHROM) + "," + str(record.POS) + "," + str(record.REF) + "," + str(record.ALT[0])
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
            time.sleep(30)
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
                return {'protein_id': list(df['PROTEIN_ID']),
                        'score': list(df['SCORE']),
                        'pred': list(df['PREDICTION (cutoff=-2.5)'])}
            except:
                sys.stderr.write("attempt {} failed...\n".format(attempt))
                time.sleep(30)
            else:
                sys.stderr.write('max attempts\n')
                break
    return None


def classify_pathogenicity(record):
    is_loh = record.INFO["facetsLCN_EM"] == 0 if 'facetsLCN_EM' in record.INFO else False
    ann_effect = [x.split('|')[1] for x in record.INFO['ANN']]
    hap_insuf = 'hap_insuf' in record.INFO
    cp = filter(lambda x: x.endswith('chasm_score'), record.INFO.keys())
    chasm_scores = [min(record.INFO[x]) for x in cp]
    is_chasm_pathogenic = any([x <= 0.3 for x in chasm_scores])
    is_fathmm_pathogenic = record.INFO['fathmm_pred'] == "CANCER" if 'fathmm_pred' in record.INFO else False
    is_cancer_gene = 'lawrence' in record.INFO or 'kandoth' in record.INFO
    is_mt_pathogenic = False
    if 'MutationTaster_pred' in record.INFO:
        is_mt_pathogenic = 'disease' in record.INFO['MutationTaster_pred']
    elif 'dbNSFP_MutationTaster_pred' in record.INFO:
        is_mt_pathogenic = record.INFO['dbNSFP_MutationTaster_pred'] == 'D' or \
            record.INFO['dbNSFP_MutationTaster_pred'] == 'A'

    if any([c in ef for ef in ann_effect for c in ["frameshift", "splice_donor", "splice_acceptor", "stop_gained"]]):
        if (is_loh or hap_insuf) and is_cancer_gene:
            record.INFO["pathogenicity"] = "pathogenic"
        elif is_loh or hap_insuf or is_cancer_gene:
            record.INFO["pathogenicity"] = "potentially_pathogenic"
        else:
            record.INFO["pathogenicity"] = "passenger"
    elif any(["missense_variant" in ef for ef in ann_effect]):
        if ~is_mt_pathogenic and ~is_chasm_pathogenic:
            record.INFO["pathogenicity"] = "passenger"
        else:
            if is_fathmm_pathogenic or is_chasm_pathogenic:
                record.INFO["pathogenicity"] = "pathogenic" if is_cancer_gene else "potentially_pathogenic"
            else:
                record.INFO["pathogenicity"] = "passenger"
    elif any(["inframe" in ef for ef in ann_effect]):
        if ~is_mt_pathogenic:
            provean_res = query_provean(record, 500)
            if provean_res is None:
                record.INFO["pathogenicity"] = "unknown"
            else:
                is_provean_pathogenic = any([x == 'Deleterious' for x in provean_res['pred']])
                record.INFO['provean_protein_id'] = provean_res['protein_id']
                record.INFO['provean_pred'] = provean_res['pred']
                record.INFO['provean_score'] = provean_res['score']
        if ~is_mt_pathogenic and provean_res is not None and ~is_provean_pathogenic:
            record.INFO["pathogenicity"] = "passenger"
        else:
            if (is_loh or hap_insuf) and is_cancer_gene:
                record.INFO["pathogenicity"] = "pathogenic"
            elif is_loh or hap_insuf or is_cancer_gene:
                record.INFO["pathogenicity"] = "potentially_pathogenic"
            else:
                record.INFO["pathogenicity"] = "passenger"

if __name__ == "__main__":
    for record in vcf_reader:
        classify_pathogenicity(record)
        vcf_writer.write_record(record)
    vcf_writer.close()
