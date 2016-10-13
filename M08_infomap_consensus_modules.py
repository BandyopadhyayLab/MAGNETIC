# coding: utf-8
import argparse
import csv
import itertools

import numpy as np
import scipy.spatial.distance as sdist
import scipy.cluster.hierarchy as sch

from S00_common_io import read_infomap_clusters


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--clusters', nargs='+')
    parser.add_argument('--labels')

    parser.add_argument('--min_size', type=int, default=3)
    parser.add_argument('--output')

    args = parser.parse_args()

    with open(args.labels) as f:
        genes = [line.strip() for line in f]


    clusters = []

    for input_file in args.clusters:
        cs = read_infomap_clusters(input_file, genes, deepest=True)
        clusters.append({c:cs[c] for c in cs if len(cs[c]) >= args.min_size})

    clustered_genes = sorted({g for cs in clusters for c in cs for g in cs[c]})

    print len(clustered_genes)

    gene_d = {g:i for i,g in enumerate(clustered_genes)}

    nc2 = lambda n: n * (n - 1) / 2 # n-choose-2 function
    gc2 = nc2(len(gene_d)) # n-choose-2 for number of genes
    ij = lambda i,j: gc2 - nc2(len(gene_d) - i) + (j - i - 1) # go from matrix coords to vector

    gene_m = np.zeros(gc2, dtype=float)

    for cs in clusters:
        for c in cs:
            for g1,g2 in itertools.combinations(sorted(cs[c]), 2):
                gene_m[ij(gene_d[g1], gene_d[g2])] += 1.0

    gene_m /= len(clusters)

    gene_m = sdist.squareform(gene_m)
    gene_m += np.eye(gene_m.shape[0])

    with open(args.output + '.mat', 'w') as OUT:
        gene_m.tofile(OUT)

    with open(args.output + '_genes.txt', 'w') as OUT:
        print >> OUT, '\n'.join(clustered_genes)

    modules = set()
    for i in range(len(clustered_genes)):
        idx = gene_m[i, :] == 1
        m_g = [g for i,g in zip(idx, clustered_genes) if i]

        if len(m_g) >= args.min_size:
            modules.add(tuple(m_g))

    module_list = sorted((sorted(m) for m in modules), key=len, reverse=True)

    with open(args.output + '_modules.txt', 'w') as OUT:
        print >> OUT, '\n'.join('module_{}\t{}'.format(i, '\t'.join(m_g)) for i,m_g in enumerate(module_list))
