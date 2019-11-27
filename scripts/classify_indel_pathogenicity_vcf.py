#!/usr/bin/env python

""" classify pathogenicity of vcf records, querying provean as necessary
"""

import vcf
import argparse
import sys
import remote_provean_query as rpq
import mutation_taster_query as mtq
import classify_pathogenicity_vcf as cp


def query_mutation_taster(record):
    try:
        prediction, score = mtq.query(record.CHROM, record.POS, record.REF, record.ALT[0])
        record.INFO['MT_pred'] = prediction
        record.INFO['MT_score'] = score
    except:
        record.INFO['MT_pred'] = 'none'



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
    parser.add_argument('--no_remote', action='store_true', default=True, help='no remote queries: can only call potentially pathogenic')
    parser.add_argument('--no_mt_provean', action='store_true', default=False, help='do not run mutation taster / provean')
    args = parser.parse_args()

    vcf_reader = vcf.Reader(open(args.vcf_infile, 'r'))

    assert "hap_insuf" in vcf_reader.infos
    assert "ANN" in vcf_reader.infos
    assert "kandoth" in vcf_reader.infos
    assert "lawrence" in vcf_reader.infos
    assert "facetsLOH" in vcf_reader.infos or "LOH" in vcf_reader.infos

    # add necessary info headers
    vcf_reader.infos['pathogenicity'] = vcf.parser._Info(id='pathogenicity', num=1, type='String',
                                                         desc="Classification of pathogenicity",
                                                         source=None, version=None)
    vcf_reader.infos['provean_protein_id'] = vcf.parser._Info(id='provean_protein_id', num=1, type='String',
                                                              desc="provean protein id (run if necessary)",
                                                              source=None, version=None)
    vcf_reader.infos['provean_pred'] = vcf.parser._Info(id='provean_pred', num=1, type='String',
                                                        desc="provean prediction (run if necessary)",
                                                        source=None, version=None)
    vcf_reader.infos['provean_score'] = vcf.parser._Info(id='provean_score', num=1, type='Float',
                                                         desc="provean score (run if necessary)",
                                                         source=None, version=None)
    vcf_reader.infos['MT_pred'] = vcf.parser._Info(id='MT_pred', num=1, type='String',
                                                        desc="mutation taster prediction (run if necessary)",
                                                        source=None, version=None)
    vcf_reader.infos['MT_score'] = vcf.parser._Info(id='MT_score', num=1, type='Float',
                                                         desc="mutation taster score (run if necessary)",
                                                         source=None, version=None)

    vcf_writer = vcf.Writer(sys.stdout, vcf_reader)

    records = list()
    query_records = list()
    for record in vcf_reader:
        if len(record.FILTER) == 0 and record.is_indel and cp.requires_mt_provean(record):
            query_records.append(record)
        records.append(record)
    if len(query_records) > 0 and not args.no_remote:
        # run mutation taster
        for record in records:
            if 'MT_pred' not in record.INFO:
                query_mutation_taster(record)
        if args.run_local:
            import provean_query
            query_manager = provean_query.ProveanQueryManager(query_records, args.provean_script,
                                                              args.num_provean_threads, args.mem_per_thread,
                                                              cluster_mode=args.cluster_mode)
            query_manager.run_queries()
        else:
            query = rpq.RemoteProveanQuery(query_records)
            query.run_query()

    for record in records:
        if len(record.FILTER) == 0 and record.is_indel:
            cp.classify_pathogenicity(record, no_remote=args.no_remote)
        vcf_writer.write_record(record)
    vcf_writer.close()
