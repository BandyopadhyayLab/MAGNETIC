__author__ = 'james'

import argparse
import csv
import itertools
import os

import numpy as np
import scipy.stats as st


# calculate for edge score t in input_file 0 ... max?
# t -> median correlation in BCL of all edges > t
def proc_file(tcga_file, bcl_file,
              tcga_d0_genes, tcga_d1_genes,
              bcl_d0_genes, bcl_d1_genes,
              output_file, same_d=False):

    tcga_data = np.fromfile(tcga_file, dtype=np.float32).reshape(
        (len(tcga_d0_genes), len(tcga_d1_genes)))
    bcl_data = np.fromfile(bcl_file, dtype=np.float32).reshape(
        (len(bcl_d0_genes), len(bcl_d1_genes)))

    d0_int = set(tcga_d0_genes) & set(bcl_d0_genes)
    d1_int = set(tcga_d1_genes) & set(bcl_d1_genes)

    tcga_i0 = np.array([(g in d0_int) for g in tcga_d0_genes], dtype=bool)
    tcga_i1 = np.array([(g in d1_int) for g in tcga_d1_genes], dtype=bool)

    bcl_i0 = np.array([(g in d0_int) for g in bcl_d0_genes], dtype=bool)
    bcl_i1 = np.array([(g in d1_int) for g in bcl_d1_genes], dtype=bool)

    tcga_data = tcga_data[np.ix_(tcga_i0, tcga_i1)]
    bcl_data = bcl_data[np.ix_(bcl_i0, bcl_i1)]
    assert tcga_data.shape == bcl_data.shape

    thresholds = np.linspace(-1.0, 1.0, 41)

    if same_d:
        tcga_data[np.eye(tcga_data.shape[0], dtype=bool)] = thresholds[0]

    stats = []
    for i,t in enumerate(thresholds[:-1]):
        idx = (tcga_data > t) & (tcga_data <= thresholds[i+1])
        if idx.sum() == 0:
            continue

        stats.append((idx.sum(), t, thresholds[i+1], np.mean(bcl_data[idx]), np.std(bcl_data[idx]))
                     + tuple(np.percentile(bcl_data[idx], (1.0, 2.5, 5.0, 25.0, 50.0, 75.0, 95.0, 97.5, 99.0))))

    with open(output_file, 'w') as OUT:
        print >> OUT, 'count\tt0\tt1\tmean\tstd\t1\t2.5\t5\t25\t50\t75\t95\t97.5\t99'
        for s in stats:
            print >> OUT, '\t'.join('{:g}'.format(st) for st in s)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # matrix of edge scores
    parser.add_argument('--tcga_matrix', nargs='+')
    # matrix of edge correlations
    parser.add_argument('--bcl_matrix_dir')

    parser.add_argument('--tcga_labels', nargs='+')
    parser.add_argument('--bcl_labels', nargs='+')

    parser.add_argument('--output_dir')

    parser.add_argument('--j', type=int)

    args = parser.parse_args()

    tcga_matrix = args.tcga_matrix[args.j]

    d0,d1 = os.path.basename(tcga_matrix).split('_')[0].split('-')
    if d0 > d1:
        import sys
        sys.exit()

    bcl_matrix = os.path.join(args.bcl_matrix_dir, '{}-{}_matrix.dat'.format(d0,d1))
    assert os.path.exists(bcl_matrix)


    tcga_label_files = {os.path.basename(lf).split('_')[2]:lf for lf in args.tcga_labels}
    bcl_label_files = {os.path.basename(lf).split('_')[2]:lf for lf in args.bcl_labels}

    with open(tcga_label_files[d0]) as f:
        tcga_d0_genes = sorted({line.strip() for line in f})

    with open(bcl_label_files[d0]) as f:
        bcl_d0_genes = sorted({line.strip() for line in f})

    if d0 != d1:
        with open(tcga_label_files[d1]) as f:
            tcga_d1_genes = sorted({line.strip() for line in f})

        with open(bcl_label_files[d1]) as f:
            bcl_d1_genes = sorted({line.strip() for line in f})
    else:
        tcga_d1_genes = tcga_d0_genes
        bcl_d1_genes = bcl_d0_genes

    proc_file(tcga_matrix, bcl_matrix,
              tcga_d0_genes, tcga_d1_genes,
              bcl_d0_genes, bcl_d1_genes,
              os.path.join(args.output_dir, '{}-{}.txt'.format(d0, d1)),
              same_d=(d0 == d1))
