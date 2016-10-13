# coding: utf-8
import argparse
import csv
import itertools
import os

from collections import defaultdict

import numpy as np


def feval(d, i):
    return np.var(d * i[:, None], 0).sum()


def brute_force(data, n_urows, rrows):
    i2v = lambda i: np.array([(-1)**int(j) for j in bin(i)[2:].zfill(n_urows)])
    best_i = min(xrange(2**(n_urows - 1)), key=lambda i: feval(data, i2v(i)[rrows]))
    return i2v(best_i)[rrows]


def mhmc(datac, n_urows, rrows, n_runs):
    temperatures = [300., 100., 30., 10., 3., 1., 0.3, 0.1, 0.03, 0.01]

    # initialize i randomly
    current_i = (-1.0)**np.random.randint(0, 2, n_urows)
    current_score = feval(datac, current_i[rrows])
    best_score = current_score
    best_i = current_i.copy()

    for T in temperatures:
        for i in range(n_urows * n_runs):
            m = np.random.randint(1, n_urows / 2 + 1)
            j = np.random.choice(np.arange(n_urows), m, replace=False)
            current_i[j] *= -1
            new_score = feval(datac, current_i[rrows])
            if new_score <= current_score:
                current_score = new_score
                if new_score < best_score:
                    best_score = new_score
                    best_i[:] = current_i
            else: # accept with MH probability
                if np.random.random() < np.exp((current_score - new_score) / T):
                    current_score = new_score
                else:
                    current_i[j] *= -1

    return best_i[rrows]


def bounded_bfs(datac, n_rows, r):
    idx = np.hstack([[-1.], np.ones(n_rows - 1)])

    ci = np.ones(n_rows)
    cs = feval(datac, ci)

    q = [(cs, set(), ci)]

    n_row_s = frozenset(range(n_rows))

    k = 0

    while k < len(q):
        cs,n,ci = q[k]
        k += 1

        if len(n) > (n_rows / 2):
            continue

        for s in (n_row_s - n):
            ni = ci * np.roll(idx, s)
            ns = feval(datac, ni)

            if ns < cs:
                q.append((ns, n | {s}, ni))


    best_s,best_n,best_i = min(q, key=lambda t: t[0])

    return best_i


def greedy(datac, n_rows):
    current_i = (-1.0)**np.random.randint(0, 2, n_rows)
    current_score = feval(datac, current_i)

    idx = np.hstack([[-1.], np.ones(n_rows - 1)])

    fi = lambda i: feval(current_i * np.roll(idx, i))

    for i in xrange(n_rows * 10):
        j = min(range(n_rows), key=fi)
        current_i[j] *= -1
        new_score = feval(datac, current_i)
        if current_score <= new_score:
            current_i[j] *= -1
            break
        else:
            current_score = new_score

    return current_i


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('clusters', nargs='+')
    parser.add_argument('--output_dir', default=None)
    parser.add_argument('--trials', type=int, default=1)
    parser.add_argument('--runs', type=int, default=100)
    parser.add_argument('--seed', type=int, default=None)

    parser.add_argument('--i', type=int, default=None)

    args = parser.parse_args()

    if args.seed is not None:
        np.random.seed(args.seed)

    if args.i is not None:
        clusters = [sorted(args.clusters,
                           key=lambda fn: int(fn[:-4].split('_')[-1]))[args.i]]
    else:
        clusters = args.clusters


    for cluster_file in clusters:
        with open(cluster_file) as f:
            rdr = csv.reader(f, delimiter='\t')
            h = rdr.next()[1:]
            rows = list(rdr)

        lbls = [row[0] for row in rows]

        urows = dict()
        rrows = []

        for i,row in enumerate(rows):
            urow = tuple(row[1:])
            if urow not in urows:
                urows[urow] = len(urows)
            rrows.append(urows[urow])

        data = np.array([[float(v) if v != 'NA' else -9999. for v in row[1:]]
                         for row in rows])
        rrows = np.array(rrows)
        n_urows = len(urows)

        ix = (data == -9999).sum(0) == 0
        h = [c for i,c in enumerate(h) if ix[i]]
        data = data[:, ix]

        if n_urows < 24:
            results = [brute_force(data, n_urows, rrows)]
        else:
            results = []

            for k in range(args.trials):
                r = (-1)**np.random.randint(0, 2, n_urows)[rrows]

                datac = data * r[:, None]

                # results.append(greedy(n_urows))
                results.append(mhmc(datac, n_urows, rrows, args.runs) * r)
                # results.append(bounded_bfs(datac, n_urows, r) * r)

        best_i = min(results, key=lambda i: feval(data, i))
        print len(set(map(lambda i: feval(data, i), results))), feval(data, best_i)

        output_file = os.path.join(args.output_dir, os.path.basename(cluster_file)[:-4])
        if args.seed is not None:
            output_file += '_{}.dat'.format(args.seed)
        else:
            output_file += '.dat'

        with open(output_file, 'w') as OUT:
            best_i.astype(int).tofile(OUT)
