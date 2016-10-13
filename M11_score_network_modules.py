__author__ = 'james'

# coding: utf-8
import argparse
import csv
import os

from collections import defaultdict

import numpy as np

from S00_common_io import read_input


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--cluster_edges', nargs='+')
    parser.add_argument('--clusters', nargs='+')
    parser.add_argument('--ref_clusters', nargs='*')
    parser.add_argument('--output')

    parser.add_argument('--strict', action='store_true')

    args = parser.parse_args()

    h_set = set()
    cluster_scores = []

    cluster_files = {os.path.basename(fn):fn for fn in args.clusters}
    edge_files = {os.path.basename(fn):fn for fn in args.cluster_edges}
    if args.ref_clusters:
        ref_clusters = {os.path.basename(fn):fn for fn in args.ref_clusters}
    else:
        ref_clusters = False

    for fn in sorted(set(cluster_files) & set(edge_files) & set(ref_clusters)):
        with open(edge_files[fn]) as f:
            rdr = csv.reader(f, delimiter='\t')
            edges = [(row[0], row[1], float(row[2]), row[3].split('-')) for row in rdr]

        if ref_clusters:
            with open(ref_clusters[fn]) as f:
                rdr = csv.reader(f, delimiter='\t')
                h = rdr.next()
                flips = [r[0].split('_') for r in rdr]
                flips = {(d,g):int(flip[1:-1]) for d,g,flip in flips}
        else:
            flips = defaultdict(lambda: 1)

        with open(cluster_files[fn]) as f:
            rdr = csv.reader(f, delimiter='\t')
            h = rdr.next()[1:]
            h_set.add(tuple(h))

            lbls,data = zip(*[(row[0], [float(v) if v != 'NA' else -9999 for v in row[1:]])
                              for row in rdr])
            dts = {lbl.split('_', 1)[0] for lbl in lbls}

            data_i = {dt:{lbl.split('_', 2)[1]:i for i,lbl in enumerate(lbls) if lbl.startswith(dt)}
                      for dt in dts}
            data = np.ma.masked_values(data, -9999)

            score = np.zeros(len(h))

            for e1,e2,w,dt in edges:
                for d in dt:
                    for e in (e1,e2):
                        if d in data_i and e in data_i[d]:
                            if args.strict:
                                score = score + (w * flips[d,e] * data[data_i[d][e], :])
                            else:
                                score += (w * flips[d,e] * data[data_i[d][e], :]).filled(0)

        if args.strict:
            cluster_scores.append((os.path.splitext(fn)[0],
                                   {h[i]:('NA' if (np.ma.is_masked(score) and score.mask[i]) else ('%f' % score[i]))
                                    for i in range(score.shape[0])}))
        else:
            cluster_scores.append((os.path.splitext(fn)[0],
                                   {h[i]:('%f' % score[i]) for i in range(score.shape[0])}))

    h_u = sorted(reduce(set.union, h_set, set()))

    with open(args.output, 'w') as OUT:
        print >> OUT, 'cluster\t%s' % '\t'.join(h_u)
        print >> OUT, '\n'.join('{}\t{}'.format(cn, '\t'.join(d.get(c, 'NA') for c in h_u))
                                for cn,d in cluster_scores)
