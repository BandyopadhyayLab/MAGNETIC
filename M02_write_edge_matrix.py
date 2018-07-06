import argparse
import csv
import itertools
import os

import numpy as np


def read_dtype_file(file_name):
    with open(file_name) as f:
        rdr = csv.reader(f, delimiter='\t')
        h = rdr.next()[1:]
        labels = [row[0] for row in rdr]

    assert list(labels) == sorted(labels)

    return h,labels


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--input', nargs='+')
    parser.add_argument('--labels', nargs='+')
    parser.add_argument('--output')
    parser.add_argument('--i', type=int)
    parser.add_argument('--j', type=int, default=0)

    parser.add_argument('--p', type=float)

    args = parser.parse_args()

    test = args.output.format('blah', 'blah', args.p, args.j)

    np.random.seed(args.j)

    dat_1,dat_2 = list(itertools.combinations_with_replacement(sorted(args.input), 2))[args.i]
    input_1,input_2 = list(itertools.combinations_with_replacement(sorted(args.labels), 2))[args.i]

    d1 = os.path.splitext(os.path.basename(input_1))[0].rsplit('_', 1)[1]
    d2 = os.path.splitext(os.path.basename(input_2))[0].rsplit('_', 1)[1]

    assert d1 == os.path.splitext(os.path.basename(dat_1))[0].rsplit('_', 2)[1]
    assert d2 == os.path.splitext(os.path.basename(dat_2))[0].rsplit('_', 2)[1]
    assert d1 <= d2

    h1,labels1 = read_dtype_file(input_1)
    n_d1 = len(labels1)
    data1 = np.memmap(dat_1, mode='r', dtype=np.float64, shape=(n_d1, len(h1)))

    if d1 != d2:
        h2,labels2 = read_dtype_file(input_2)
        n_d2 = len(labels2)
        data2 = np.memmap(dat_2, mode='r', dtype=np.float64, shape=(n_d2, len(h2)))

        h_int = set(h1) & set(h2)
        h1_idx = np.array([(c in h_int) for c in h1])
        h2_idx = np.array([(c in h_int) for c in h2])

        assert np.array_equal(np.array(h1)[h1_idx], np.array(h2)[h2_idx])

        data1 = data1[:, h1_idx]
        data2 = data2[:, h2_idx]
        assert data1.shape[1] == data2.shape[1]

        if args.p < 1.0:
            n_sub = int(np.round(args.p * data1.shape[1]))
            subsample = np.random.choice(np.arange(data1.shape[1]), replace=False, size=n_sub)
        else:
            subsample = slice(None)

        nz1 = np.any(data1[:, subsample] != 0.0, 1)
        nz2 = np.any(data2[:, subsample] != 0.0, 1)

        cc = np.corrcoef(data1[np.ix_(nz1, subsample)],
                         data2[np.ix_(nz2, subsample)])[:nz1.sum(), nz1.sum():]

        out = np.memmap(args.output.format(d1, d2, args.p, args.j),
                        dtype=np.float64, mode='w+', shape=(n_d1, n_d2))

        out[np.ix_(nz1, nz2)] = cc
        out.flush()
    else:
        if args.p < 1.0:
            n_sub = int(np.round(args.p * data1.shape[1]))
            subsample = np.random.choice(np.arange(data1.shape[1]), replace=False, size=n_sub)
        else:
            subsample = slice(None)

        nz1 = np.any(data1[:, subsample] != 0.0, 1)
        cc = np.corrcoef(data1[np.ix_(nz1, subsample)])

        out = np.memmap(args.output.format(d1, d1, args.p, args.j),
                        dtype=np.float64, mode='w+', shape=(n_d1, n_d1))

        out[np.ix_(nz1, nz1)] = cc
        out.flush()
