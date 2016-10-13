# coding: utf-8
import argparse
import csv
import itertools
import os

import numpy as np


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--cluster_files', nargs='+')
    parser.add_argument('--cluster_edges', nargs='+')
    parser.add_argument('--dats', nargs='+')
    parser.add_argument('--output_dir')

    args = parser.parse_args()

    # file_list = zip(sorted(args.clusters), sorted(args.cluster_edges))
    # assert all(os.path.basename(f1) == os.path.basename(f2) for f1,f2 in file_list)

    edge_files = {os.path.basename(fn):fn for fn in args.cluster_edges}

    for cluster_file in args.cluster_files:
        with open(cluster_file) as f:
            rdr = csv.reader(f, delimiter='\t')
            h = rdr.next()[1:]
            rows = list(rdr)

        lbls = [row[0] for row in rows]
        data = np.array([[float(v) if v != 'NA' else -9999. for v in row[1:]]
                         for row in rows])
        ix = (data == -9999).sum(0) == 0
        h = [c for i,c in enumerate(h) if ix[i]]
        data = data[:, ix]
        n_rows = data.shape[0]

        cluster_is = set()

        dateq = lambda f1,f2: (os.path.basename(f1)[:-4] == os.path.basename(f2)[:-4]
                               or (len(os.path.basename(f1).split('_')) == 3
                                   and os.path.basename(f1).rsplit('_', 1)[0] == os.path.basename(f2)[:-4]))

        dat_files = [df for df in args.dats if dateq(df, cluster_file)]

        if not dat_files:
            continue

        for dat_file in dat_files:
            d = np.fromfile(dat_file, dtype=int)
            cluster_is.add(tuple(d.tolist()))


        if len(cluster_is) <= 2:
            flips = np.array(max(cluster_is, key=lambda i: np.abs(sum(i))))
        else:
            import matplotlib.pyplot as plt
            flips = np.vstack([np.array(fs) for fs in cluster_is])

            with open(edge_files[os.path.basename(cluster_file)]) as f:
                rdr = csv.reader(f, delimiter='\t')
                edges = [(row[0], row[1], float(row[2]), row[3].split('-')) for row in rdr]

            dts = {lbl.split('_', 1)[0] for lbl in lbls}
            data_i = {dt:{lbl.split('_', 2)[1]:i for i,lbl in enumerate(lbls) if lbl.startswith(dt)} for dt in dts}

            scores = np.array([np.vstack([w * flips[i,data_i[d][e]] * data[data_i[d][e], :]
                                          for e1,e2,w,dt in edges
                                          for d in dt for e in (e1,e2)
                                          if d in data_i and e in data_i[d]]).sum(0)
                               for i in range(flips.shape[0])])

            scores = scores[:, scores.mean(0).argsort()]

            assert np.all(np.abs(np.corrcoef(scores)) > 0.99)

            print cluster_file
            plt.plot(scores.T)
            plt.show()

            min_flips = np.argmax(np.abs(flips.sum(1)))
            flips = flips[min_flips, :]


        if flips.sum() < 0:
            flips *= -1

        data *= flips[:, None]

        # print (data * flips[:, None]).var(0).sum()

        output_file = os.path.join(args.output_dir, os.path.basename(cluster_file))

        with open(output_file, 'w') as OUT:
            print >> OUT, 'row_id\t{}'.format('\t'.join(h))
            print >> OUT, '\n'.join('{}_({})\t{}'.format(lbls[i], flips[i],
                                                         '\t'.join('%f' % v for v in data[i,:]))
                                    for i in range(n_rows))

