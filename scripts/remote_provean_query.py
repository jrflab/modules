import urllib2
from mechanize import Browser
from bs4 import BeautifulSoup
import sys
import time
import re
import pandas as pd

class TooManyQueriesException(Exception):
    pass

class RemoteProveanQuery:
    def __init__(self, records, rest_server='http://grch37.rest.ensembl.org'):
        self.records = records
        self.rest_server = rest_server
        self.query = ""
        for record in self.records:
            self.query += "{},{},{},{}\n".format(record.CHROM, record.POS, record.REF, record.ALT[0])

    def _get_job_url(self, max_retry=10):
        for attempt in range(max_retry):
            try:
                br = Browser()
                sys.stderr.write("Querying Provean:\n{}".format(self.query))
                br.open('http://provean.jcvi.org/genome_submit_2.php?species=human')
                br.form = list(br.forms())[1]  # select the chrpos form
                control = br.form.find_control("CHR")
                control.value = self.query
                br.submit()
                job_url = br.geturl()
                sys.stderr.write("job url: {}\n".format(job_url))
                if 'genome_prg_2' in job_url:
                    raise TooManyQueriesException("too many provean queries in 24 hr")
                if 'jobid' not in job_url:
                    raise ValueError("jobid not in job url")
                return job_url
            except TooManyQueriesException:
                raise
            except Exception:
                sys.stderr.write("query attempt {} failed...\n".format(attempt))
                time.sleep(10)
            else:
                raise ValueError('max query attempts to provean')

    def _parse_job_results(self, job_url, max_retry=10):
        for attempt in range(max_retry):
            try:
                page = urllib2.urlopen(job_url).read()
                soup = BeautifulSoup(page, 'html.parser')
                link = soup.find('a', href=re.compile('one\.tsv'))
                url = 'http://provean.jcvi.org/' + link.get('href')
                df = pd.read_table(url)
                return df
            except:
                sys.stderr.write("attempt {} failed...\n".format(attempt))
                time.sleep(10)
            else:
                raise ValueError('max query attempts to provean')

    def run_query(self, max_retry=10):
        job_url = self._get_job_url()
        df = self._parse_job_results(job_url)
        for idx, record in enumerate(self.records):
            if df.ix[idx, 'PROTEIN_ID'] != 'record not found':
                #record.INFO['provean_protein_id'] = df.ix[idx, 'PROTEIN_ID']
                record.INFO['provean_pred'] = df.ix[idx, 'PREDICTION (cutoff=-2.5)']
                record.INFO['provean_score'] = df.ix[idx, 'SCORE']
            else:
                #record.INFO['provean_protein_id'] = '.'
                record.INFO['provean_pred'] = 'none'
                record.INFO['provean_score'] = '.'
