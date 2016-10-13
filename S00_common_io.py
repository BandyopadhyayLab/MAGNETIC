__author__ = 'james'

import csv
import os

from collections import defaultdict
from multiprocessing import Pool

import numpy as np


# reads a giant struct file in big chunks to speed up IO
def bytes_from_file(f, chunksize, bitesize):
    chunk = ''
    while True:
        chunk += f.read(chunksize)
        if chunk:
            for i in xrange(0, len(chunk) - bitesize + 1, bitesize):
                yield chunk[i:i+bitesize]
            chunk = chunk[i+bitesize:]
        else:
            break


def read_input(input_file):
    labels = []
    data = []

    dtype_starts = dict()
    dtype_ends = dict()

    with open(input_file) as f:
        rdr = csv.reader(f, delimiter='\t')
        rdr.next()

        dtype = None
        dtypes_seen = set()

        for i,row in enumerate(rdr):
            new_dtype,lbl = row[0].split('_', 1)

            labels.append(lbl)
            data.append(row[1:])

            if new_dtype != dtype:
                assert new_dtype not in dtypes_seen

                dtypes_seen.add(new_dtype)

                dtype_starts[new_dtype] = i + 1
                if dtype is not None:
                    dtype_ends[dtype] = i

                dtype = new_dtype

        dtype_ends[dtype] = i + 1

    data = np.array([[float(v) for v in row] for row in data])

    dtypes = sorted(dtypes_seen)

    return (labels, data, dtype_starts, dtype_ends, dtypes)


def read_stats(stats_file):
    edge_means = dict()
    edge_stds = dict()

    mins = []

    with open(stats_file) as f:
        rdr = csv.reader(f, delimiter='\t')

        for row in rdr:
            d1,d2 = row[0],row[1]
            mu = float(row[5])
            edge_means[(d1,d2)] = mu
            edge_means[(d2,d1)] = mu

            dev = float(row[6])
            edge_stds[(d1,d2)] = dev
            edge_stds[(d2,d1)] = dev

            mins.append((float(row[3]) - mu) / dev)

    dist_pad = -np.floor(min(mins))

    return edge_means, edge_stds, dist_pad


# just some reusable code to create a Pool of size proc (or length of the queue, if it's less),
# but only if there is more than one job (in which case, avoiding a Pool is good for debugging)
def adaptive_pool(func, jobs, proc, maxtasksperchild=None, chunksize=1):
    if proc == 1 or (hasattr(jobs, '__len__') and len(jobs) == 1):
        return map(func, jobs)
    else:
        mpool = Pool(processes=min(proc, len(jobs) if hasattr(jobs, '__len__') else proc),
                     maxtasksperchild=maxtasksperchild)
        return mpool.map(func, jobs, chunksize=chunksize)


def read_mcl_clusters(fn, genes, label_fn=None, random=False):
    if label_fn is not None:
        with open(label_fn) as f:
            label_d = dict(line.split() for line in f)
    else:
        label_d = {str(i):i for i in range(len(genes))}

    with open(fn) as f:
        clusters = defaultdict(set)
        line = f.readline()
        while line[0] != '0':
            line = f.readline()
        while line[0] != ')':
            s = line.split()
            if line[0] == ' ':
                clusters[c].update(int(label_d[i]) for i in s if i != '$')
            else:
                c = s[0]
                clusters[c].update(int(label_d[i]) for i in s[1:] if i != '$')
            line = f.readline()

    clusters = {c:clusters[c] for c in clusters if len(clusters[c]) > 2}

    genes_n,genes_g = zip(*[(n,genes[n]) for n in reduce(set.union, clusters.values(), set())])

    if random:
        genes_g = list(genes_g)
        np.random.shuffle(genes_g)

    genes_d = dict(zip(genes_n, genes_g))

    clusters = {c:{genes_d[n] for n in clusters[c]} for c in clusters}

    return clusters


def read_infomap_clusters(input_file, genes, deepest=False, random=False):
    if (not deepest) and random:
        raise NotImplementedError('Ehhhhh')

    with open(input_file) as f:
        f.readline()
        rows = [line.strip().split() for line in f]

    rows2 = [(tuple(row[0].split(':')), int(row[3])) for row in rows]

    n_levels = max(len(r[0]) for r in rows2)

    if n_levels == 1:
        # no clusters here
        return dict()

    clusters = defaultdict(lambda: defaultdict(set))
    deepness = defaultdict(int)

    for r in rows2:
        for i in range(1, n_levels):
            if len(r[0]) > i:
                clusters[i][r[0][:i]].add(r[1])
                deepness[r[1]] = max(deepness[r[1]], i)

    if deepest:
        clusters = {c:{n for n in clusters[i][c] if deepness[n]  == i} for i in clusters for c in clusters[i]}
        clusters = {c:clusters[c] for c in clusters if len(clusters[c]) > 2}

        genes_n,genes_g = zip(*[(n,genes[n]) for n in reduce(set.union, clusters.values(), set())])

        if random:
            genes_g = list(genes_g)
            np.random.shuffle(genes_g)

        genes_d = dict(zip(genes_n, genes_g))

        clusters = {c:{genes_d[n] for n in clusters[c]} for c in clusters}
    else:
        clusters = {i:{c:clusters[i][c] for c in clusters[i] if len(clusters[i][c]) > 2} for i in clusters}

    return clusters


def get_c(n, d='BCL', c=3, fl=True):
    p = {('BCL',3): '/Volumes/webberj/20150723_BCL/i3_data',
         ('TCGA',3): '/Volumes/webberj/20150722_TCGA/i3_data'}

    with open(os.path.join(p[d,c] + ('_flipped_data' if fl else ''),
                           'cluster_{}.txt'.format(n))) as f:
        rdr = csv.reader(f, delimiter='\t')
        h = rdr.next()[1:]
        rows = {row[0]:np.array([(float(v) if v != 'NA' else -9999) for v in row[1:]])
                for row in rdr}

    return h,rows
