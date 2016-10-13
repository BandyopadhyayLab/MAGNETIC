__author__ = 'james'

# coding: utf-8
import argparse
import csv
import os

from collections import defaultdict

import numpy as np

from S00_common_io import read_input


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--data', nargs='+')

    parser.add_argument('--module_edges', nargs='+')
    parser.add_argument('--output_dir')

    parser.add_argument('--all_data', action='store_true',
                        help='Get all available data, rather than graph-only')

    args = parser.parse_args()

    h_u = set()
    lines = dict()

    for label_file in args.data:
        dtype = os.path.splitext(label_file)[0].split('_')[-1]
        with open(label_file) as f:
            rdr = csv.DictReader(f, delimiter='\t')
            h_u.update(rdr.fieldnames[1:])
            lines[dtype] = list(rdr)

    h_u = sorted(h_u)

    for edge_file in args.module_edges:
        with open(edge_file) as f:
            rdr = csv.reader(f, delimiter='\t')
            edges = [((row[0], row[1]), row[3].split('-')) for row in rdr]

        if args.all_data:
            module = {(d,g) for gs,ds in edges for g in gs for d in lines}
        else:
            module = {(d,g) for gs,ds in edges for d,g in zip(ds, gs)}

        with open(os.path.join(args.output_dir,
                               os.path.basename(edge_file)), 'w') as OUT:
            print >> OUT, 'row_id\t{}'.format('\t'.join(h_u))
            print >> OUT, '\n'.join('{}_{}\t{}'.format(d, line['row_id'],
                                                       '\t'.join(line.get(c, 'NA') for c in h_u))
                                    for d in sorted(lines) for line in lines[d]
                                    if (d, line['row_id']) in module)
