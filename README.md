# MAGNETIC
The repository for the MAGNETIC project from the Bandyopadhyay lab.

### Workflow

This file should cover all of the steps to get from the public data files to a list of modules and module-drug correlations. When I refer to the "previous workflow" I mean what I did to get my current modules&mdash;since then I have experimented a little with improving the method but haven't seen enough to change. Still, the improvements don't hurt and are likely better.

I'm naming all the Python scripts as they are used. When something is run on the cluster there's also a shell script used to submit it, I will list those at the end.

1. Get the data (from TCGA, CCLE/CTRP)
2. Process data (convert all of it to HGNC gene IDs, filter and normalize)
  * Two scripts for this: `process_files.py` and `process_files_BCL.py` for TCGA and cell lines respectively.
  * Each row was normalized by z-scoring the values: (x - µ) / σ
  * I discard samples/features with lots of missing values, and use k-nearest-neighbors to infer other missing values from similar samples.
3. Calculate all-by-all correlation
  * I use `write_edge_matrix.py` to do this on the cluster&mdash;it's not that slow but it's a memory hog.
  * Calculates Pearson's &rho; for every pair of data types and every pair of genes.
4. Enrichment scoring
  * `matrix_edge_enrichment_bootstrap.py` (on cluster) takes a matrix of correlation values and bootstraps an enrichment score using a network of interactions.
  * `corr2enrichment.py` takes the output of the above and makes some nice plots.
  * `format_enrichment_file.py` reformats the output into one file to use for the next step.
  * `bin2eval.py` (on cluster) converts the correlation matrices into enrichment matrices.
5. Creating network files
  * `write_multiplex_network.py` takes all of those enrichment matrices and writes them into one big text file, for use with Infomap.
     * We apply a cutoff to this network, based on how well the resulting clusters enrich for known interactions. So I make a bunch of different network files with different edge cutoffs to take into the next step.
     * Currently I'm using a cutoff of 3, I don't know whether that will be same in a new cancer type.
     * I experimented with different ways of representing the network but it seems like a "flat" representation (one node per gene, multiple edges between nodes) works the best.
6. Running Infomap (runs on cluster, uses lots of memory)
  * Infomap (www.mapequation.org) is a network clustering algorithm that use random walks to find tightly connected areas of the network.
  * It takes a network file from the previous step and outputs a hierarchical clustering.
  * Because it has a stochastic component, I run it many times and calculate a consensus clustering once I've chosen the best parameters.
7. Picking a threshold based on enrichment
  * `network_overlap.py` (on cluster) takes a directory of Infomap outputs (not the consensus) as well as a known network (e.g. HumanNet) and calculates the enrichment for known edges in the identified modules.
    * This value is corrected for the inherent enrichment of a network at a given threshold
    * The script takes a --random argument to compute the background distribution
  * `combine_neto_files.py` and `plot_network_overlap.py` are scripts for processing the output, and make a bar plot of enrichment. I picked a threshold where it seemed to level off.
8. Consensus Infomap modules
  * `infomap_consensus_modules.py` makes a .dat file of the consensus matrix, as well as a gene list, and outputs the consensus clusters
9. Calculating module scores
  * First I use `extract_module_data.py` to make a file per module that includes all the data used in that module. For each node, if a given data type had no edges in the graph, it won't be included.
  * To account for the fact that some nodes are anti-correlated with the others, I flip some values (multiply by -1) to minimize the variance across nodes.
    * This is not trivially solvable (I think) so I use the script `flip_module.py` and (for big clusters that require lots of trials) `flip_module_proc.py`. These get run on the cluster.
    * `flip_module.py` just checks all flip possibilities for small modules, which is O(2<sup>n</sup>). The other way is to do simulated annealing which requires many tries but they tend to be similar in the end.
  * A module score is a weighted sum of the node values in the modules. The weights are determined by the edges the node has.
    * `score_network_modules.py` takes the module edges and flipped module data and calculates the scores (this is pretty fast)
10. Merging gene modules
  * A lot of module scores are highly correlated, even though they are not connected by any edges.
  * Our solution to this was just to merge all modules that correlate at a high level&mdash;I calculated all-by-all correlation of the module scores and merged anything that was in the top 0.5% of correlations.
  * After doing this, I have to calculate the module scores again as in step 9.
11. Module analysis 0: Identifying/annotating modules
  * There isn't one right way to do this&mdash;we've used gene-enrichment tools and done correlation/ANOVA against clinical covariates try to figure out what the modules mean.
  * Another likely annotation is amplicons: I check to see whether a module is all or mostly in a single genomic area. The majority of the BRCA modules were amplicons of some kind. This is less common in LUAD though.
12. Module analysis 1: Clustering of TCGA samples
 * `consensus_cluster.py` (run on cluster) takes the matrix of module scores across samples and calculates consensus clustering for a given number of clusters _k_
13. Module analysis 2: Correlation with drug response


Copyright 2016 Regents of University of California
