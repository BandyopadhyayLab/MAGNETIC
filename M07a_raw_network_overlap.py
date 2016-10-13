__author__ = 'james'


import argparse
import csv

import networkx
import numpy as np



def get_overlap(G, input_network, genes, random=False):
    if random:
        np.random.shuffle(genes)

    with open(input_network) as f:
        rdr = csv.reader(f, delimiter='\t')
        overlap = sum(G.has_edge(genes[int(r[0])], genes[int(r[1])])
                      for r in rdr)

    return overlap


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--input_network', nargs='*')

    parser.add_argument('--labels')

    parser.add_argument('--network')
    parser.add_argument('--cutoff', type=float, default=2.0)

    parser.add_argument('--random', action='store_true')
    parser.add_argument('--random_seed', type=int, default=0)

    parser.add_argument('--i', type=int)

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

    overlap = get_overlap(G, args.input_network[args.i], genes, random=args.random)

    print '{}\t{:d}\t{:d}'.format(args.input_network[args.i], args.random_seed, overlap)
