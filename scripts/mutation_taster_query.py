import pandas as pd
import sys


_prediction_priority = ['disease_causing_automatic',
                        'disease_causing',
                        'polymorphism_automatic',
                        'polymorphism']


def query(chrom, pos, ref, alt):
    pred = 'none'
    score = None
    url = "http://www.mutationtaster.org/cgi-bin/MutationTaster/MT_ChrPos.cgi" \
        "?chromosome={}&position={}&ref={}&alt={}".format(chrom, pos, ref, alt)
    sys.stderr.write("Querying: {}\n".format(url))
    dfs = pd.read_html(url)
    if (len(dfs) < 2):
        raise Exception('query failed')
    else:
        # summary is the second dataframe
        summary = dfs[1][1:]
        summary.columns = dfs[1].iloc[0]
        if 'prediction' in summary.columns:
            pred_score = {}
            for i, row in summary.iterrows():
                if row['prediction'] not in pred_score:
                    pred_score[row['prediction']] = []
                    pred_score[row['prediction']].append(float(row['probability']))
                    for p in _prediction_priority:
                        if p in pred_score:
                            pred = p
                            score = max(pred_score[p])
                            break
    return (pred, score)
