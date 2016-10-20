#!/usr/bin/env python

""" classify pathogenicity of vcf records, querying provean as necessary
"""

import vcf
import argparse
import sys
import remote_provean_query as rpq
import mutation_taster_query as mtq


def is_fs_splice_stop(record):
    ann_effect = get_ann_effect(record)
    return any([c in ef for ef in ann_effect for c in ["frameshift", "splice_donor", "splice_acceptor", "stop_gained"]])


def get_fs_splice_stop_pathogenicity(record):
    if (is_loh(record) or is_hap_insuf(record)) or is_cancer_gene(record):
        return "likely_pathogenic"
    else:
        return "passenger"


def query_mutation_taster(record):
    prediction, score = mtq.query(record.CHROM, record.POS, record.REF, record.ALT[0])
    record.INFO['MutationTaster_pred'] = [prediction]
    record.INFO['MutationTaster_score'] = [score]


def is_mt_pathogenic(record):
    return any(['disease' in info for info in record.INFO['MutationTaster_pred']])


def is_chasm_pathogenic(record):
    cp = filter(lambda x: x.endswith('chasm_score'), record.INFO.keys())
    chasm_scores = [min(record.INFO[x]) for x in cp]
    return any([x <= 0.3 for x in chasm_scores])


def is_cancer_gene(record):
    return 'lawrence' in record.INFO or 'kandoth' in record.INFO or 'cancer_gene_census' in record.INFO


def is_hap_insuf(record):
    return 'hap_insuf' in record.INFO


def get_ann_effect(record):
    return [x.split('|')[1] for x in record.INFO['ANN']]


def is_loh(record):
    return record.INFO["facetsLCN_EM"] == 0 if 'facetsLCN_EM' in record.INFO else False


def is_provean_pathogenic(record):
    return 'provean_pred' in record.INFO and record.INFO['provean_pred'] == 'Deleterious'


def get_inframe_pathogenicity(record):
    if (is_loh(record) or is_hap_insuf(record) or is_cancer_gene(record)) and \
            is_mt_pathogenic(record) and is_provean_pathogenic(record):
        return 'likely_pathogenic'
    else:
        return 'passenger'


def classify_pathogenicity(record):
    """ classify pathogenicity for a vcf record
    """
    if is_fs_splice_stop(record):
        record.INFO["pathogenicity"] = get_fs_splice_stop_pathogenicity(record)
    elif is_inframe(record):
        record.INFO["pathogenicity"] = get_inframe_pathogenicity(record)
    else:
        record.INFO["pathogenicity"] = None


def is_inframe(record):
    ann_effect = get_ann_effect(record)
    return any(["inframe" in ef for ef in ann_effect])


def requires_mt_provean(record):
    return is_inframe(record) and (is_loh(record) or is_hap_insuf(record) or is_cancer_gene(record))


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
    query_records = list()
    for record in vcf_reader:
        if requires_mt_provean(record):
            query_records.append(record)
        records.append(record)
    if len(query_records) > 0:
        # run mutation taster
        for record in records:
            if 'MutationTaster_pred' not in record.INFO:
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
        classify_pathogenicity(record)
        vcf_writer.write_record(record)
    vcf_writer.close()
