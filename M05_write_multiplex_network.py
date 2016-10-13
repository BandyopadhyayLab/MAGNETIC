# coding: utf-8
import argparse
import gzip
import os
import itertools

import numpy as np

def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks"""
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return itertools.izip_longest(fillvalue=fillvalue, *args)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('input', nargs="+")
    parser.add_argument('--labels', nargs="+")
    parser.add_argument('--output')

    parser.add_argument('--power', type=float, default=1.0)
    parser.add_argument('--cutoff', type=float, default=0.0)

    parser.add_argument('--format', choices=('multiplex', 'monolayer', 'summed'))
    parser.add_argument('--cutoff_type', choices=('independent', 'group'))

    args = parser.parse_args()

    labels = dict()
    for label_file in args.labels:
        d1 = os.path.basename(label_file).split('_')[2]
        with open(label_file) as f:
            labels[d1] = {line.strip():i for i,line in enumerate(f)}

    label_u = sorted(reduce(set.union, labels.values(), set()))

    mmaps = dict()
    for input_file in args.input:
        d1,d2 = os.path.basename(input_file).split('_')[0].split('-')

        mmaps[(d1,d2)] = np.memmap(input_file, dtype=np.float32,
                                   mode='r', shape=(len(labels[d1]), len(labels[d2])))

    dtypes = list(enumerate(sorted(reduce(set.union, mmaps, set()))))

    with gzip.open(args.output, 'w') as OUT:
        if args.format == 'multiplex':
            print >> OUT, '*Vertices {:d}'.format(len(label_u))
            print >> OUT, '\n'.join('{:d} "{}"'.format(i,g) for i,g in enumerate(label_u))
            print >> OUT, '*Multiplex'

        for g1g2 in grouper(itertools.combinations(enumerate(label_u), 2), 10000):
            lines = []
            for (i1,g1),(i2,g2) in itertools.ifilter(None, g1g2):
                pair_edges = []
                for (j1,d1),(j2,d2) in itertools.product(dtypes, dtypes):
                    if g1 in labels[d1] and g2 in labels[d2]:
                        if j1 <= j2:
                            e = mmaps[(d1,d2)][labels[d1][g1], labels[d2][g2]]
                        else:
                            e = mmaps[(d2,d1)][labels[d2][g2], labels[d1][g1]]

                        if ((args.cutoff_type == 'independent' and e > args.cutoff)
                            or (args.cutoff_type == 'group' and e > 0.0)):

                            pair_edges.append((j1, j2, e))

                if pair_edges and (args.cutoff_type != 'group'
                                   or sum(e for j1,j2,e in pair_edges) > args.cutoff):
                    if args.format == 'multiplex':
                        lines.extend('{:d}\t{:d}\t{:d}\t{:d}\t{:f}'.format(j1, i1, j2, i2, e ** args.power)
                                     for j1,j2,e in pair_edges)
                    elif args.format == 'monolayer':
                        lines.extend('{:d}\t{:d}\t{:f}'.format(i1, i2, e ** args.power)
                                     for j1,j2,e in pair_edges)
                    elif args.format == 'summed':
                        lines.append('{:d}\t{:d}\t{:f}'.format(i1, i2, sum(e for j1,j2,e in pair_edges) ** args.power))

            if lines:
                print >> OUT, '\n'.join(lines)
