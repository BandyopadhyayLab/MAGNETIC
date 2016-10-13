import argparse
import csv
import os
import re

from collections import defaultdict

import numpy as np


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--matrix', nargs='+', help='Correlation matrices')
    parser.add_argument('--enrichment', nargs='+',
                        help='Enrichment files. The filenames matter!')
    parser.add_argument('--i', type=int, default=None)
    parser.add_argument('--output')

    args = parser.parse_args()

    fnre = re.compile(r'^(.+-.+)_matrix_(?:(n)egative_)?enrichment.txt', flags=re.I)
    get_f = lambda fn: fnre.match(os.path.basename(fn)).groups()

    fns = {get_f(fn): [(float(row['Bin']), float(row['Enrichment_median']))
                       for row in csv.DictReader(open(fn), delimiter='\t')]
           for fn in args.enrichment}

    all_evals = dict()
    bins = defaultdict(list)

    for d,n in fns:
        if n is None:
            assert [b for b, e in fns[d, n]] == sorted([b for b, e in fns[d, n]])
        else:
            assert [b for b, e in fns[d, n]] == sorted([b for b, e in fns[d, n]], reverse=True)

        bins[d].extend(b for b, e in fns[d, n])

    for d in bins:
        bins[d].sort()
        all_evals[d] = [fns[d, 'n'][0][1]]

        for b, e in fns[d, 'n'][1:]:
            all_evals[d].append(max(e, all_evals[d][-1]))
        all_evals[d] = all_evals[d][::-1]

        all_evals[d].append(fns[d, None][0][1])
        for b, e in fns[d, None][1:]:
            all_evals[d].append(max(e, all_evals[d][-1]))

    all_bins = {tuple(bins[d]) for d in bins}

    assert len(all_bins) == 1
    all_bins = np.array(all_bins.pop())

    assert all((len(all_evals[d]) == len(all_bins)) for d in all_evals)

    evalues = {d:np.log(all_evals[d]) for d in all_evals}

    if args.i is not None:
        matrix_files = [args.matrix[args.i]]
    else:
        matrix_files = args.matrix


    for input_file in matrix_files:
        d = os.path.basename(input_file).split('_')[0]

        d_evalues = evalues[d]
        zero_val = d_evalues[all_bins == 0.0].mean()

        input_data = np.fromfile(input_file, dtype=np.float32)

        for i in xrange(input_data.size):
            if input_data[i] < 0.0:
                input_data[i] = d_evalues[(all_bins < input_data[i]).sum()]
            elif input_data[i] > 0.0:
                input_data[i] = d_evalues[(all_bins <= input_data[i]).sum() - 1]
            else:
                input_data[i] = zero_val

        with open(args.output.format(d), 'wb') as OUT:
            input_data.tofile(OUT)
