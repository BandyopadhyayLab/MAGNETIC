# coding: utf-8
import argparse
import csv
import itertools
import os

from collections import Counter, defaultdict

import numpy as np

from S00_common_io import read_input


def translate_genes(gene_list, mapped, unmapped, hgnc_file):
    gene_list = {g.upper() for g in gene_list}

    if isinstance(hgnc_file, str):
        with open(hgnc_file) as f:
            rdr = csv.DictReader(f, delimiter='\t')
            hgnc_genes = [row for row in rdr]
    else:
        hgnc_genes = hgnc_file

    hgnc_set = {row['symbol'].upper() for row in hgnc_genes}

    old2new = defaultdict(set)
    syn2new = defaultdict(set)
    for row in hgnc_genes:
        if row['prev_symbol']:
            for o in row['prev_symbol'].split('|'):
                old2new[o.upper()].add(row['symbol'].upper())

        if row['alias_symbol']:
            for s in row['alias_symbol'].split('|'):
                if s.startswith('h'):
                    if s[1].isupper():
                        syn2new[s[1:].upper()].add(row['symbol'].upper())
                    elif s[2] == '-':
                        syn2new[s[2:].upper()].add(row['symbol'].upper())
                syn2new[s.upper()].add(row['symbol'].upper())

    # throw out ambiguous mappings
    old2new = {o:list(old2new[o])[0] for o in old2new
               if len(old2new[o]) == 1}
    syn2new = {s:list(syn2new[s])[0] for s in syn2new
               if len(syn2new[s]) == 1}
    ambiguous = {a for a in set(old2new).intersection(syn2new) if old2new[a] != syn2new[a]}

    translated_genes = dict()
    lost_genes = set()

    # just going to say: current ID trumps old, via-old trumps via-synonym
    # if a symbol could map to two genes, I throw it out
    for g in gene_list:
        if g in hgnc_set:
            translated_genes[g] = g
        elif g in old2new and g not in ambiguous:
            translated_genes[g] = old2new[g]
        elif g in syn2new and g not in ambiguous:
            translated_genes[g] = syn2new[g]
        else:
            lost_genes.add(g)

    print '%d entries mapped to new names' % sum(k != v for k,v in translated_genes.items())
    print '%d entries could not be mapped' % len(lost_genes)

    # print 'old', Counter(len(v) for v in old2new.values())
    # print 'syn', Counter(len(v) for v in syn2new.values())
    # print "'%s'" % max(old2new, key=lambda k: len(old2new[k]))
    # print "'%s'" % max(syn2new, key=lambda k: len(syn2new[k]))

    # print out the new values
    if mapped:
        with open(mapped, 'w') as OUT:
            print >> OUT, '\n'.join(sorted(set(translated_genes.values())))

    if unmapped:
        with open(unmapped, 'w') as OUT:
            print >> OUT, '\n'.join(sorted(lost_genes))

    return translated_genes


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--input', nargs='*')
    parser.add_argument('--i', type=int, default=None)
    parser.add_argument('--output')

    parser.add_argument('--labels')
    parser.add_argument('--hgnc_file')

    parser.add_argument('--mapped', default=None)
    parser.add_argument('--unmapped', default=None)
    # parser.add_argument('--gene_union', default=None)

    args = parser.parse_args()

    genef = lambda g: g.split('_')[0].upper()


    genes,data,dtype_starts,dtype_ends,dtypes = read_input(args.labels)

    print 'total', len(set(genef(lbl) for lbl in genes)),
    gene_translate = translate_genes({genef(lbl) for lbl in genes},
                                     mapped=(args.mapped % 'all') if args.mapped else None,
                                     unmapped=(args.unmapped % 'all') if args.unmapped else None,
                                     hgnc_file=args.hgnc_file)

    g_sets = dict()
    g_translated = dict()

    for d in dtypes:
        g_sets[d] = {genef(lbl) for lbl in genes[(dtype_starts[d] - 1):dtype_ends[d]]}
        print d, len(g_sets[d]),
        g_translated[d] = translate_genes(g_sets[d],
                                          mapped=(args.mapped % d) if args.mapped else None,
                                          unmapped=(args.unmapped % d) if args.unmapped else None,
                                          hgnc_file=args.hgnc_file)


    # these are the *new* gene names
    gene_union = sorted(set(gene_translate.values()))
    gene_union_d = {g:i for i,g in enumerate(gene_union)}
    print 'union of translated genes', len(gene_union)

    d2i = defaultdict(dict)
    for d in dtypes:
        for i in range(dtype_starts[d] - 1, dtype_ends[d]):
            g = genef(genes[i])
            if g in gene_translate:
                # index in dtype matrix -> the index in the unified matrix
                d2i[d][i - dtype_starts[d] + 1] = gene_union_d[gene_translate[g]]


    if args.i is not None:
        input_files = [args.input[args.i]]
    else:
        input_files = args.input

    for input_file in input_files:
        d1,d2 = os.path.basename(input_file).split('_')[0].split('-')

        pairwise = np.zeros((len(gene_union), len(gene_union)), dtype=np.float32)
        pair_counts = np.zeros_like(pairwise, dtype=np.int16)

        n_d1 = dtype_ends[d1] - dtype_starts[d1] + 1
        n_d2 = dtype_ends[d2] - dtype_starts[d2] + 1

        input_data = np.fromfile(input_file, dtype=np.float32).reshape((n_d1,n_d2))
        npairs = 0

        for i1,i2 in np.ndindex(*input_data.shape):
            if i1 in d2i[d1] and i2 in d2i[d2]:
                npairs += 1
                pairwise[d2i[d1][i1],d2i[d2][i2]] += input_data[i1,i2]
                pair_counts[d2i[d1][i1],d2i[d2][i2]] += 1
                if d2i[d1][i1] != d2i[d2][i2]:
                    pairwise[d2i[d2][i2],d2i[d1][i1]] += input_data[i1,i2]
                    pair_counts[d2i[d2][i2],d2i[d1][i1]] += 1


        pairwise[pair_counts > 0] /= pair_counts[pair_counts > 0]

        print ('%s-%s\t%s\t%d\t%g\t%g\t%d\t%d'
               % (d1, d2, pairwise.shape,
                  (pairwise != 0).sum(), pairwise.mean(), np.std(pairwise),
                  (input_data != 0).sum(), npairs))

        with open(args.output % (d1,d2), 'wb') as OUT:
            pairwise.tofile(OUT)
