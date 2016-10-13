# coding: utf-8
import argparse
import csv
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument('--drug_data', default='gb-2013-14-10-r110-s1.txt')
parser.add_argument('--output_file')
parser.add_argument('--output_raw_file')

args = parser.parse_args()

with open(args.drug_data, 'rU') as f:
    rdr = csv.DictReader(f, delimiter='\t')
    drugs = rdr.fieldnames[11:]

    gbd = {row['Cell line']:{d:row[d] for d in drugs} for row in rdr}

cls = sorted(gbd)
cl_count = {c:sum(gbd[c][d] == 'NA' for d in drugs) for c in cls}
cls = [c for c in sorted(gbd) if cl_count[c] < 48]

zd = lambda d: (d[d != -9999] - d[d != -9999].mean()) / d[d != -9999].std()

gb_data = {d:np.array([(float(gbd[c][d]) if gbd[c][d] != 'NA' else -9999) for c in cls]) for d in drugs}
gb_zdata = dict()

for d in gb_data:
    dd =  gb_data[d].copy()
    dd[dd != -9999] = zd(dd)
    gb_zdata[d] = dd

with open(args.output_file, 'w') as OUT:
    print >> OUT, 'Drug\t{}'.format('\t'.join(cls))
    print >> OUT, '\n'.join('{}\t{}'.format(d, '\t'.join(str(v) for v in gb_zdata[d])) for d in drugs)

with open(args.output_raw_file, 'w') as OUT:
    print >> OUT, 'Drug\t{}'.format('\t'.join(cls))
    print >> OUT, '\n'.join('{}\t{}'.format(d, '\t'.join(str(v) for v in gb_data[d])) for d in drugs)
