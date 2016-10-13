# coding: utf-8
import argparse
import csv
import numpy as np

from collections import defaultdict, Counter

from sklearn import preprocessing
from sklearn.neighbors import NearestNeighbors

from S01_hgnc_translate import translate_genes


# clean and normalize (to uppercase) the gene labels
genef = lambda g: g.split('_')[0].upper()

# function to process the tumor IDs into a uniform format
hsh = lambda h,c: '-'.join(h.split(c))[:15] # this seems fragile but works as of 20150427!


# function to process mutation data
def process_muts(rdr, cutoff):
    var_set = {'Missense_Mutation',
               'Nonsense_Mutation',
               'Frame_Shift_Del',
               'In_Frame_Del',
               'Frame_Shift_Ins',
               'Nonstop_Mutation',
               'In_Frame_Ins'}

    muts = defaultdict(set)

    for row in rdr:
        if row['Variant_Classification'] in var_set:
            muts[row['Hugo_Symbol']].add(hsh(row['Tumor_Sample_Barcode'], '-'))

    h = reduce(set.union, muts.values(), set())
    cutoff *= len(h)
    muts = {k:muts[k] for k in muts if len(muts[k]) > cutoff}
    h2 = sorted(reduce(set.union, muts.values(), set()))

    return h2,muts


def get_rows(rdr, i, j, hgnc_file, id_dict=None):
    if id_dict:
        labels,rows = zip(*[(genef(id_dict[row[j]]),row[i:]) for row in rdr if row[j] in id_dict])
    else:
        labels,rows = zip(*[(genef(row[j]),row[i:]) for row in rdr])

    rows = [[float(v) if v not in ('', 'NA') else -9999.0 for v in row]
            for row in rows]

    gene_translate = translate_genes(set(labels), mapped=None, unmapped=None,
                                     hgnc_file=hgnc_file)

    gene_counts = defaultdict(list)

    for lbl,row in zip(labels,rows):
        if lbl in gene_translate:
            gene_counts[gene_translate[lbl]].append(row)

    labels2,rows2 = zip(*[(g,max(gene_counts[g], key=lambda r: sum(v for v in r if v != -9999.)))
                          for g in sorted(gene_counts)])

    return labels2,rows2


def process_rows(labels, hdrs, rows):
    data = np.array(rows)

    print ((data == -9999.0).sum(0) > 0).sum(), 'samples with missing values'
    sample_i = ((data == -9999.0).sum(0) <= (0.2 * len(hdrs)))
    print (len(hdrs) - sample_i.sum()), 'samples with >= 20% missing values'

    print ((data == -9999.0).sum(1) > 0).sum(), 'rows with missing values'
    gene_i = ((data == -9999.0).sum(1) <= (0.2 * len(labels)))
    print (len(labels) - gene_i.sum()), 'rows with >= 20% missing values'

    # filter out rows/columns with > 20% missing values (not many)
    labels2 = [labels[i] for i in np.nonzero(gene_i)[0]]
    hdrs2 = [hdrs[i] for i in np.nonzero(sample_i)[0]]
    data2 = data[np.ix_(gene_i, sample_i)]

    assert data2.shape == (len(labels2), len(hdrs2))

    # robust scaling of the data
    for i in range(data2.shape[0]):
        if np.any(data2[i,:] == -9999):
            data2[i, data2[i,:] != -9999] = preprocessing.robust_scale(data2[i, data2[i,:] != -9999].reshape((1,-1)), 1)
        else:
            data2[i,:] = preprocessing.robust_scale(data2[[i],:].reshape((1,-1)), 1)

    # find genes that are present in all samples
    # a slightly-better approach would adjust this for every missing value, but this will be faster
    gene_js = ((data[:,sample_i] == -9999).sum(1) == 0.0)

    for i in range(data2.shape[0]):
        if np.any(data2[i,:] == -9999):
            # find all samples that have this value
            # note: using original copy, don't want to include imputed stuff
            sample_j = (data[i,sample_i] != -9999)

            data3 = data2[np.ix_(gene_js, sample_j)].T
            # find the nearest neighbor for each missing value
            knn = NearestNeighbors(n_neighbors=20, algorithm='brute').fit(data3)

            for j in np.nonzero(data2[i,:] == -9999)[0]:
                data2[i,j] = np.median(data2[i, knn.kneighbors(data2[gene_js,j].reshape((1,-1)),
                                                               return_distance=False).ravel()])

    assert (data2 == -9999).sum() == 0

    # re-order sample labels to be in alphabetical order
    hdrs2_i = {h:i for i,h in enumerate(hdrs2)}
    h_idx = [hdrs2_i[h] for h in sorted(hdrs2)]

    return labels2,sorted(hdrs2),data2[:,h_idx]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--source', choices=('TCGA', 'BCL', 'metabric'), default='TCGA')
    parser.add_argument('--output', help='Output file name template')

    # input files
    parser.add_argument('--exp', default='sblab/20130712_TCGA/BRCA.exp.547.med.txt')
    parser.add_argument('--cnv', default='sblab/20130712_TCGA/GISTIC2/all_data_by_genes.txt')
    parser.add_argument('--methyl', default='20140110_TCGA/BRCA.methylation.27k.450k.filtered_new.txt')
    parser.add_argument('--rppa', default='sblab/20130712_TCGA/rppaData-403Samp-171Ab-Trimmed.txt')
    parser.add_argument('--muts', default='sblab/20130712_TCGA/genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.3.2.0.somatic.maf')

    parser.add_argument('--rppa_map', default='20140110_TCGA/rppa_id_map.txt',
                        help='Map from antibody names to gene symbols')
    parser.add_argument('--methyl_map', default='2013_GrayLab/synapse_data/LBNL_Gray_BCCL_methylation_annotation_v1.txt',
                        help='Map from methylation array ID to gene symbols (BCL only)')
    parser.add_argument('--mut_cutoff', type=float, default=0.02,
                        help='Minimum mutation frequency to consider')
    parser.add_argument('--hgnc_file', default='HGNC/20150923_protein_coding_gene.txt',
                        help='File of HGNC symbols')

    args = parser.parse_args()

    assert args.output.format('blah', 'blah')

    with open(args.hgnc_file, 'rU') as f:
        rdr = csv.DictReader(f, delimiter='\t')
        hgnc_genes = [row for row in rdr]

    with open(args.rppa_map, 'rU') as f:
        rppa_id = {row[0]:row[1] for row in csv.reader(f, delimiter='\t')}

    with open(args.methyl_map, 'rU') as f:
        methyl_id = {row['Name']:row['Symbol'] for row in csv.DictReader(f, delimiter='\t')}

    # name, file, first data column, gene id map, sample id separator
    if args.source == 'TCGA':
        output_tuples = [('cnv', args.cnv, 3, 0, None, '-'),
                         ('exp', args.exp, 1, 0, None, '-'),
                         ('methyl', args.methyl, 1, 0, None, '-'),
                         ('rppa', args.rppa, 1, 0, rppa_id, '.')]
    elif args.source == 'BCL':
        output_tuples = [('cnv', args.cnv, 5, 4, None, None),
                         ('exp', args.exp, 1, 0, None, None),
                         ('methyl', args.methyl, 1, 0, methyl_id, None),
                         ('rppa', args.rppa, 1, 0, rppa_id, None)]
    elif args.source == 'metabric':
        output_tuples = [('cnv', args.cnv, 1, 0, None, None),
                         ('exp', args.exp, 1, 0, None, None)]

    dtype_h = dict()

    for dtype,fname,i,j,id_map,delim in output_tuples:
        with open(fname, 'rU') as f:
            rdr = csv.reader(f, delimiter='\t')
            if args.source == 'TCGA':
                h = [hsh(c,delim) for c in rdr.next()[i:]]
            else:
                h = rdr.next()[i:]

            labels,rows = get_rows(rdr, i, j, hgnc_genes, id_map)
            print dtype, len(h), 'columns',
            labels,h,data = process_rows(labels, h, rows)
            dtype_h[dtype] = set(h)

            with open(args.output.format(args.source, dtype), 'w') as OUT:
                print >> OUT, "row_id\t%s" % '\t'.join(h)
                for i,g in enumerate(labels):
                    print >> OUT, '%s\t%s' % (g, '\t'.join('%.5f' % data[i,j] for j in range(len(h))))

    if args.source == 'TCGA':
        with open(args.muts) as f:
            f.readline()
            rdr = csv.DictReader(f, delimiter='\t')
            h,muts = process_muts(rdr, args.mut_cutoff)

            dtype_h['mut'] = set(h)

            with open(args.output.format(args.source, 'mut'), 'w') as OUT:
                print >> OUT, 'row_id\t%s' % '\t'.join(h)
                for g in sorted(muts):
                    print >> OUT, '%s\t%s' % (g, '\t'.join('%d' % int(c in muts[g]) for c in h))

    for fn in glob.glob(args.output.format(args.source, '*')):
        with open(fn) as f:
            rdr = csv.reader(f, delimiter='\t')
            rdr.next()
            data = np.array([map(float, r[1:]) for r in rdr])
            with open(fn[:-4] + '.dat', 'w') as OUT:
                data.tofile(OUT)

    print 'dtype sample overlaps:'
    print 'dtype\t%s' % '\t'.join(sorted(dtype_h))
    for d0 in sorted(dtype_h):
        print '%s\t%s' % (d0, '\t'.join('%d' % len(dtype_h[d0] & dtype_h[d1]) for d1 in sorted(dtype_h)))
