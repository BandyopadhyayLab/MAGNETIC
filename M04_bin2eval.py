import argparse
import csv
import os

import numpy as np


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--matrix', nargs='+', help='Correlation matrices')
    parser.add_argument('--enrichment_dir', help='Location of enrichment files')
    parser.add_argument('--i', type=int)
    parser.add_argument('--output_dir')

    args = parser.parse_args()

    matrix_file = args.matrix[args.i]
    matrix_base = os.path.splitext(os.path.basename(matrix_file))[0]

    positive_enrichment_file = os.path.join(
        args.enrichment_dir, '{}_enrichment.txt'.format(matrix_base)
    )
    negative_enrichment_file = os.path.join(
        args.enrichment_dir, '{}_negative_enrichment.txt'.format(matrix_base)
    )

    output_file = os.path.join(
        args.output_dir, '{}_evalues.dat'.format(matrix_base)
    )

    with open(positive_enrichment_file) as f:
       pos_enr = [(float(row['Bin']), float(row['Enrichment_median']))
                  for row in csv.DictReader(f, delimiter='\t')]

    with open(negative_enrichment_file) as f:
        neg_enr = [(float(row['Bin']), float(row['Enrichment_median']))
                   for row in csv.DictReader(f, delimiter='\t')]

    assert [b for b,e in pos_enr] == sorted([b for b,e in pos_enr])
    assert [b for b,e in neg_enr] == sorted([b for b,e in neg_enr], reverse=True)

    bins = np.array(sorted([b for enr in (neg_enr, pos_enr) for b,e in enr]))

    all_evals = [max(neg_enr[0][1], 0.5)]

    for b,e in neg_enr[1:]:
        all_evals.append(max(e, all_evals[-1]))
    all_evals = all_evals[::-1]

    all_evals.append(max(pos_enr[0][1], 0.5))
    for b,e in pos_enr[1:]:
        all_evals.append(max(e, all_evals[-1]))

    d_evalues = np.log(all_evals)

    zero_val = d_evalues[bins == 0.0].mean()

    input_data = np.fromfile(matrix_file, dtype=np.float64)

    for i in xrange(input_data.size):
        if input_data[i] < 0.0:
            input_data[i] = d_evalues[(bins < input_data[i]).sum()]
        elif input_data[i] > 0.0:
            input_data[i] = d_evalues[(bins <= input_data[i]).sum() - 1]
        else:
            input_data[i] = zero_val

    with open(output_file, 'wb') as OUT:
        input_data.tofile(OUT)

if __name__ == '__main__':
    main()
