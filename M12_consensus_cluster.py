__author__ = 'james'

'''Consensus clustering using scikit-learn. Based on ConsensusClusterPlus
package from Bioconductor.'''

import argparse
import csv
import os

import numpy as np

import sklearn.cluster
from sklearn.cross_validation import _validate_shuffle_split
from sklearn.base import clone
from sklearn.utils import check_random_state
from sklearn.neighbors import kneighbors_graph

# class ConsensusCluster(object):
#     def __init__(self,

def consensus_cluster(
                 X, # data to be clustered [n_samples, n_features]
                 inner_estimator, # clustering object to use for repeated clustering
                 final_estimator, # clustering object to use on the consensus matrix
                 n_iter=100, # number of subsamples
                 sample_ratio=0.8, # proportion of items to sample
                 feature_ratio=1.0, # proportion of features to sample
                 sample_weight=None, # weights for sampling items
                 feature_weight=None, # weights for sampling features
                 random_state=None, # seed=NULL, # random state
                 # n_jobs=1, # number of jobs to run on
                 # verbose=False, # verbosity
                 # distance="pearson", # distance function for clustering (in estimators)
                 # n_clusters, # maxK=3, # maximum cluster number to evaluate (in estimators)
                 # innerLinkage="average", # linkage method for subsampling (in estimator)
                 # finalLinkage="average", # linkage method for consensus matrix (in estimator)
                 # d = None, # data to be clustered (should be passed to fit method)
                 # title="untitled_consensus_cluster", # title for plot (no plots)
                 # ml=NULL, # for plotting a prior result (no plots)
                 # tmyPal=NULL, # for colors of consensus matrix (no plots)
                 # plot=NULL, # not going to plot
                 # writeTable=FALSE, # will return this
                 # corUse="everything" # how to handle missing data...N/A in Python
                 ):

    rng = check_random_state(random_state)
    n_samples, n_features = X.shape

    if sample_ratio < 1.0:
        sample_sub = _validate_shuffle_split(n_samples, sample_ratio, None)[0]
    else:
        sample_sub = n_samples

    if feature_ratio < 1.0:
        feature_sub = _validate_shuffle_split(n_features, feature_ratio, None)[0]
    else:
        feature_sub = n_features

    consensus_matrix = np.zeros((n_samples, n_samples), dtype=float)
    indicator_matrix = np.zeros_like(consensus_matrix)

    for i in xrange(n_iter):
        s_i = rng.choice(n_samples, sample_sub, False, sample_weight)
        s_i.sort()
        f_i = rng.choice(n_features, feature_sub, False, feature_weight)
        f_i.sort()

        s_ix = np.ix_(s_i, s_i)

        indicator_matrix[s_ix] += 1.0
        X_sub = X[np.ix_(s_i, f_i)]

        lbls = inner_estimator.fit(X_sub).labels_
        consensus_matrix[s_ix] += (lbls == lbls[:, None])

    consensus_matrix /= indicator_matrix

    final_estimator.fit(consensus_matrix)

    return consensus_matrix


def get_method(method_name, n_clusters, connectivity=None, full_tree='auto'):
    if connectivity:
        kng = lambda d: kneighbors_graph(d, connectivity, include_self=False)
    else:
        kng = None

    if method_name == 'ward':
        return sklearn.cluster.AgglomerativeClustering(n_clusters=n_clusters,
                                                       linkage=method_name,
                                                       connectivity=kng,
                                                       compute_full_tree=full_tree)
    elif method_name == 'kmeans':
        return sklearn.cluster.KMeans(n_clusters=n_clusters)
    elif method_name == 'mbkmeans':
        return sklearn.cluster.MiniBatchKMeans(n_clusters=n_clusters)
    else:
        raise ValueError("Invalid method name: {}".format(method_name))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--input')
    parser.add_argument('--output_dir')

    parser.add_argument('--n_clusters', type=int)
    parser.add_argument('--n_iter', type=int, default=1000)

    parser.add_argument('--inner_method', default='ward', choices=('ward', 'kmeans', 'mbkmeans'))
    parser.add_argument('--outer_method', default='ward', choices=('ward', 'kmeans', 'mbkmeans'))
    parser.add_argument('--connectivity', type=int, default=None)

    args = parser.parse_args()

    if args.connectivity:
        if args.inner_method != 'ward' and args.outer_method != 'ward':
            raise ValueError("Connectivity only applies to agglomerative clustering methods")

    with open(args.input) as f:
        rdr = csv.reader(f, delimiter='\t')
        h = rdr.next()[1:]
        rows = list(rdr)

    data = np.array(zip(*[map(float, row[1:]) for row in rows]))

    inner_est = get_method(args.inner_method, args.n_clusters, args.connectivity)
    outer_est = get_method(args.outer_method, args.n_clusters, args.connectivity, True)

    consensus_matrix = consensus_cluster(data,
                                         inner_estimator=inner_est,
                                         final_estimator=outer_est,
                                         n_iter=args.n_iter)

    with open(os.path.join(args.output_dir, "consensus_labels_{}.dat".format(args.n_clusters)), 'w') as OUT:
        outer_est.labels_.tofile(OUT)

    with open(os.path.join(args.output_dir, "consensus_tree_{}.dat".format(args.n_clusters)), 'w') as OUT:
        outer_est.children_.tofile(OUT)

    with open(os.path.join(args.output_dir, "consensus_matrix_{}.dat".format(args.n_clusters)), 'w') as OUT:
        consensus_matrix.tofile(OUT)
