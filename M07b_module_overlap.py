__author__ = 'james'


import argparse
import csv
import itertools
import re

from collections import defaultdict

import networkx
import numpy as np

from S00_common_io import read_infomap_clusters


def mp_network_overlap(input):
    scp,n,c,graph,cluster = input

    cluster_size = len(cluster)
    cluster_graph = graph.subgraph(cluster)

    cluster_overlap = cluster_graph.size()

    return scp, n, c, cluster_size, cluster_overlap


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--clusters', nargs='+')
    parser.add_argument('--output')

    parser.add_argument('--labels')

    parser.add_argument('--network')
    parser.add_argument('--cutoff', type=float, default=2.0)

    parser.add_argument('--random', action='store_true')
    parser.add_argument('--random_seed', type=int, default=0)

    parser.add_argument('--i', type=int, default=0)
    parser.add_argument('--proc', type=int, default=1)

    args = parser.parse_args()

    with open(args.labels) as f:
        genes = [line.strip() for line in f]

    if args.random:
        np.random.seed(args.random_seed)

    with open(args.network) as f:
        G = networkx.Graph()
        G.add_edges_from((row[0], row[1])
                         for row in csv.reader(f, delimiter='\t')
                         if (not args.cutoff) or float(row[2]) >= args.cutoff)

    node_set = set(G.nodes()).intersection(genes)
    G = G.subgraph(node_set)

    cluster_overlap = defaultdict(dict)
    cluster_size = defaultdict(dict)

    fnre = re.compile(r'^.+output_([0-9]+)/.+_([.0-9]+)(?:_pow([0-9]+))?(?:_relax([.0-9]+))?\.tree$', re.I)
    gre = lambda fn: '\t'.join(g or '1' for g in fnre.match(fn).groups())

    cluster_list = [(gre(fn), read_infomap_clusters(fn, genes, deepest=True, random=args.random))
                    for fn in sorted(args.clusters)[args.i::args.proc]]
    mjobs = [(scp, len(clusters), c, G, clusters[c] & node_set)
             for scp,clusters in cluster_list for c in clusters]

    mjobs = [cgc for cgc in mjobs if len(cgc[4])]

    res = map(mp_network_overlap, mjobs)

    for scp,n,c,c_size,c_overlap in res:
        cluster_size[(scp,n)][c] = c_size
        cluster_overlap[(scp,n)][c] = c_overlap


    with open(args.output, 'w') as OUT:
        print >> OUT, '\t'.join(('seed', 'cutoff', 'power', 'relax',
                                 'n_clusters', 'avg_size', 'mean',
                                 'total', 'theoretical_total'))
        fmt_str = '%s\t%d\t%.3f\t%.3f\t%d\t%d'

        for scp,n in sorted(cluster_size):
            c_overlap = cluster_overlap[(scp,n)].values()
            print >> OUT, fmt_str % (scp, n,
                                     float(sum(cluster_size[(scp,n)].values())) / n,
                                     np.mean(c_overlap),
                                     sum(c_overlap),
                                     sum(v*(v-1)/2 for v in cluster_size[(scp,n)].values()))
