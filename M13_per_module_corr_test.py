__author__ = 'james'

import argparse
import csv
import itertools
import os

import numpy as np


def proc_file(cluster_files, d0, d1, bcl_data,
              tcga_d0_genes, tcga_d1_genes,
              bcl_d0_genes, bcl_d1_genes, output_file):

    d0_int = set(tcga_d0_genes) & set(bcl_d0_genes)
    d1_int = set(tcga_d1_genes) & set(bcl_d1_genes)

    bcl_i0 = np.array([(g in d0_int) for g in bcl_d0_genes], dtype=bool)
    bcl_i1 = np.array([(g in d1_int) for g in bcl_d1_genes], dtype=bool)

    bcl_data = bcl_data[np.ix_(bcl_i0, bcl_i1)]

    d0_gene_d = {g:i for i,g in enumerate(sorted(d0_int))}
    d1_gene_d = {g:i for i,g in enumerate(sorted(d1_int))}

    bkg_u = bcl_data.mean()
    bkg_d = bcl_data.std()

    with open(output_file, 'w') as OUT:
        for cluster_file in cluster_files:
            cluster_name = os.path.basename(cluster_file)[:-4]
            with open(cluster_file) as f:
                rdr = csv.reader(f, delimiter='\t')
                rows = list(rdr)

            g_cluster_edges = [(r[0], r[1]) for r in rows if r[3] == '{}-{}'.format(d0, d1)
                               and r[0] in d0_gene_d and r[1] in d1_gene_d]

            gce_set = set(g_cluster_edges)

            cluster_nodes = sorted({r[i] for r in rows for i in (0,1)})

            cn0 = [g for g in cluster_nodes if g in d0_gene_d]
            cn1 = [g for g in cluster_nodes if g in d1_gene_d]

            if d0 == d1:
                assert cn0 == cn1
                ng_cluster_edges = [(g0, g1) for g0,g1 in itertools.product(cn0, cn1)
                                    if d0_gene_d[g0] < d1_gene_d[g1] and (g0,g1) not in gce_set]
            else:
                ng_cluster_edges = [(g0, g1) for g0,g1 in itertools.product(cn0, cn1)
                                    if g0 != g1 and (g0,g1) not in gce_set]

            if len(g_cluster_edges) >= 1:
                c0,c1 = zip(*[(d0_gene_d[g0], d1_gene_d[g1]) for g0,g1 in g_cluster_edges])

                bcl_s = '{:3f}\t{:3f}'.format(bkg_u, bkg_d)
                bcl_m_s = '\t'.join('{:.4f}'.format(v) for v in bcl_data[c0, c1])

                print >> OUT, '{}\t{}-{}\tTrue\t{:d}\t{}\t{}'.format(
                    cluster_name, d0, d1, len(c0), bcl_s, bcl_m_s)

            if len(ng_cluster_edges) >= 1:
                c0,c1 = zip(*[(d0_gene_d[g0], d1_gene_d[g1]) for g0,g1 in ng_cluster_edges])

                bcl_s = '{:3f}\t{:3f}'.format(bkg_u, bkg_d)
                bcl_m_s = '\t'.join('{:.4f}'.format(v) for v in bcl_data[c0, c1])

                print >> OUT, '{}\t{}-{}\tFalse\t{:d}\t{}\t{}'.format(
                    cluster_name, d0, d1, len(c0), bcl_s, bcl_m_s)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # matrix of edge correlations
    parser.add_argument('--bcl_matrix', nargs='+')

    parser.add_argument('--tcga_labels', nargs='+')
    parser.add_argument('--bcl_labels', nargs='+')
    parser.add_argument('--cluster_edges', nargs='+')

    parser.add_argument('--output_dir')

    parser.add_argument('--j', type=int)

    args = parser.parse_args()

    bcl_matrix = args.bcl_matrix[args.j]
    output_file = os.path.join(
        args.output_dir, '{}_stats.txt'.format(os.path.basename(bcl_matrix)[:-4])
    )

    d0,d1 = os.path.basename(bcl_matrix).split('_')[0].split('-')
    if d0 > d1 or 'mut' in (d0,d1):
        import sys
        sys.exit()

    tcga_label_files = {os.path.basename(lf)[:-4].split('_')[2]:lf for lf in args.tcga_labels}
    bcl_label_files = {os.path.basename(lf)[:-4].split('_')[2]:lf for lf in args.bcl_labels}

    with open(tcga_label_files[d0]) as f:
        rdr = csv.reader(f, delimiter='\t')
        tcga_d0_genes = sorted({r[0] for r in itertools.islice(rdr, 1, None)})

    with open(bcl_label_files[d0]) as f:
        rdr = csv.reader(f, delimiter='\t')
        bcl_d0_genes = sorted({r[0] for r in itertools.islice(rdr, 1, None)})

    if d0 != d1:
        with open(tcga_label_files[d1]) as f:
            rdr = csv.reader(f, delimiter='\t')
            tcga_d1_genes = sorted({r[0] for r in itertools.islice(rdr, 1, None)})

        with open(bcl_label_files[d1]) as f:
            rdr = csv.reader(f, delimiter='\t')
            bcl_d1_genes = sorted({r[0] for r in itertools.islice(rdr, 1, None)})
    else:
        tcga_d1_genes = tcga_d0_genes
        bcl_d1_genes = bcl_d0_genes

    bcl_data = np.fromfile(bcl_matrix, dtype=np.float64).reshape(
        (len(bcl_d0_genes), len(bcl_d1_genes)))

    proc_file(args.cluster_edges, d0, d1, bcl_data,
              tcga_d0_genes, tcga_d1_genes,
              bcl_d0_genes, bcl_d1_genes,
              output_file)
