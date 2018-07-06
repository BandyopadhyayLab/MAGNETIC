__author__ = 'james'

# script to go through the network file and count up edges in humannet by threshold


import argparse
import gzip
import csv
import itertools
import os

from collections import defaultdict

import numpy as np

def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks"""
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return itertools.izip_longest(fillvalue=fillvalue, *args)


def proc_file(input_file, edges, d0_genes, d1_genes, amax, i, output_dir, n_boot):
    thresholds = np.arange(0, amax, i)

    if amax < 0.0:
        f_e = lambda e,t: (e[...,np.newaxis] <= t).sum(0)
        output_file = os.path.join(output_dir,
                                   os.path.splitext(os.path.basename(input_file))[0] + '_negative_enrichment.txt')
    else:
        f_e = lambda e,t: (e[...,np.newaxis] >= t).sum(0)
        output_file = os.path.join(output_dir,
                                   os.path.splitext(os.path.basename(input_file))[0] + '_enrichment.txt')

    if os.path.exists(output_file):
        return

    n_d0 = sum(len(v) for v in d0_genes.values())
    n_d1 = sum(len(v) for v in d1_genes.values())
    input_data = np.fromfile(input_file, dtype=np.float64).reshape((n_d0, n_d1))

    d0_egenes = {g0 for g0,g1 in edges}
    d1_egenes = {g1 for g0,g1 in edges}

    len_te = 0
    total_counts = np.zeros_like(thresholds)

    for glist in grouper(((g0,g1) for g0,g1 in itertools.product(d0_genes, d1_genes)
                          if g0 != g1 and d0_genes[g0] & d0_egenes and d1_genes[g1] & d1_egenes),
                         1000):
        glist = filter(None, glist)
        if not glist:
            continue

        len_te += len(glist)
        t_ix0,t_ix1 = zip(*[(i0,i1) for g0,g1 in glist for i0,i1
                            in itertools.product(d0_genes[g0], d1_genes[g1])])

        total_counts += f_e(input_data[t_ix0, t_ix1], thresholds)


    # convert to list for bootstrapping
    edgelist0 = np.array(list(edges), dtype=int)

    # background rate: [# of network edges] / [# possible edges]
    bkrd = float(len(edges)) / len_te

    counts_boot = np.zeros((3, thresholds.shape[0]))
    enrich_boot = np.zeros((3, thresholds.shape[0]))

    enrich_mean = np.zeros(thresholds.shape[0])
    enrich_std = np.zeros(thresholds.shape[0])

    assert thresholds.shape[0] % 10 == 0

    das_boot = np.random.randint(0, len(edgelist0), (n_boot, len(edgelist0)))

    for i in range(0, thresholds.shape[0], 10):
        # counts: number of network edges above each threshold
        counts = np.zeros((n_boot, 10), dtype=int)
        enrichment = np.zeros_like(counts, dtype=float)

        # bootstrap n_boot different edge lists
        for n in range(n_boot):
            edgeboot = das_boot[n, :]

            counts[n,:] = f_e(input_data[edgelist0[edgeboot, 0], edgelist0[edgeboot, 1]],
                              thresholds[i:i+10]).astype(float)

        counts_boot[:, i:i+10] = np.vstack(np.percentile(counts, (2.5, 50.0, 97.5), 0))

        tci = total_counts[i:i+10] > 0
        enrichment[:, tci] = (counts[:, tci] / total_counts[i:i+10][tci]) / bkrd

        enrich_boot[:, i:i+10] = np.vstack(np.percentile(enrichment, (2.5, 50.0, 97.5), 0))
        enrich_mean[i:i+10] = enrichment.mean(0)
        enrich_std[i:i+10] = enrichment.std(0)

    with open(output_file, 'w') as OUT:
        print >> OUT, '\t'.join(('Bin', 'Total',
                                 'Hnet_lo', 'HNet_median', 'Hnet_hi',
                                 'Enrichment_mu', 'Enrichment_std',
                                 'Enrichment_lo', 'Enrichment_median', 'Enrichment_hi'))

        fmt_str = '%g\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g'
        for j in range(thresholds.shape[0]):
            print >> OUT, fmt_str % (thresholds[j], total_counts[j],
                                     counts_boot[0,j], counts_boot[1,j], counts_boot[2,j],
                                     enrich_mean[j], enrich_std[j],
                                     enrich_boot[0,j], enrich_boot[1,j], enrich_boot[2,j])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('input', nargs='+')

    parser.add_argument('--output_dir')
    parser.add_argument('--labels', nargs='+')
    parser.add_argument('--network')
    parser.add_argument('--network_threshold', type=float, default=2.0)

    parser.add_argument('--i', type=float)
    parser.add_argument('--max', type=float)
    parser.add_argument('--n', type=int, default=100000)

    parser.add_argument('--j', type=int)

    parser.add_argument('--random_seed', type=int, default=0)

    args = parser.parse_args()

    if args.max < 0.0 and args.i > 0.0:
        raise ValueError("Can't mix positive index with negative range")

    np.random.seed(args.random_seed)

    input_file = args.input[args.j]
    d0,d1 = os.path.basename(input_file).split('_')[0].split('-')
    label_files = {os.path.basename(lf)[:-4].split('_')[2]:lf for lf in args.labels}

    with open(label_files[d0]) as f:
        rdr = csv.reader(f, delimiter='\t')
        rdr.next()
        d0_genes = [row[0] for row in rdr]
        assert d0_genes == sorted(d0_genes)

        d0_genes_d = defaultdict(set)
        for i,g in enumerate(d0_genes):
            d0_genes_d[g.split('_')[0]].add(i)
        # d0_genes = {g:i for i,g in enumerate(sorted(d0_genes))}

    if d0 != d1:
        with open(label_files[d1]) as f:
            rdr = csv.reader(f, delimiter='\t')
            rdr.next()
            d1_genes = [row[0] for row in rdr]
            assert d1_genes == sorted(d1_genes)

            d1_genes_d = defaultdict(set)
            for i,g in enumerate(d1_genes):
                d1_genes_d[g.split('_')[0]].add(i)
            # d1_genes = {g:i for i,g in enumerate(sorted(d1_genes))}
    else:
        # d1_genes = d0_genes
        d1_genes_d = d0_genes_d

    edges = set()

    with open(args.network) as f:
        rdr = csv.reader(f, delimiter='\t')
        for row in rdr:
            if args.network_threshold == 0.0 or float(row[2]) >= args.network_threshold:
                assert row[0] != row[1]
                if row[0] in d0_genes_d and row[1] in d1_genes_d:
                    # edges.add((d0_genes[row[0]], d1_genes[row[1]]))
                    edges.update(itertools.product(d0_genes_d[row[0]], d1_genes_d[row[1]]))
                if row[1] in d0_genes_d and row[0] in d1_genes_d:
                    # edges.add((d0_genes[row[1]], d1_genes[row[0]]))
                    edges.update(itertools.product(d0_genes_d[row[1]], d1_genes_d[row[0]]))

    proc_file(input_file, edges, d0_genes_d, d1_genes_d, # d0_genes, d1_genes
              args.max, args.i, args.output_dir, args.n)
