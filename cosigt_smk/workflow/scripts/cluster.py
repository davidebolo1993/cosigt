#!/usr/bin/python3

import os
import sys
import json
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from scipy.cluster import hierarchy
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score


def read_tsv(tsv):

	'''
	Read tsv with distances between haplotypes from odgi similarity
	'''

	return 1-pd.read_table(tsv).set_index(['group.a', 'group.b'])['jaccard.similarity'].sort_index().unstack()

def agglomerative(mtx,prefix):

	'''
	Perform agglomerative clustering and plot dendrogram
	'''

	sil_score=[]
	result=[]
	d=dict()

	agg=AgglomerativeClustering(distance_threshold=0, n_clusters=None, metric='precomputed', linkage='average')
	cluster=agg.fit(mtx)
	linkage_matrix=LinkageMatrix(agg)

	font = {'size': 8}
	plt.rc('font', **font)
	plt.figure()
	dendrogram = hierarchy.dendrogram(linkage_matrix)
	plt.ylabel('jaccard dissimilarity')
	plt.tick_params(labelbottom = False, bottom = False)
	plt.savefig(os.path.join(prefix, 'dendrogram.jaccard.dissimilarity.pdf'))
	plt.close()

	range_n_clusters = range(2,mtx.shape[0])

	for i,n_clusters in enumerate(range_n_clusters):

		clusterer = AgglomerativeClustering(n_clusters=n_clusters, metric='precomputed', linkage='average')
		cluster_labels = clusterer.fit_predict(mtx)
		silhouette_avg = silhouette_score(mtx, cluster_labels, metric='precomputed')
		sil_score.append((n_clusters,silhouette_avg))

	#get best
	k=sorted(sil_score, key=lambda x: x[1], reverse=True)[0][0]

	#plot best
	plt.figure(figsize=(10,6))
	dendrogram = hierarchy.dendrogram(linkage_matrix, truncate_mode = 'lastp', p = k)
	plt.xlabel('haplotype group')
	plt.ylabel('jaccard dissimilarity')
	plt.xticks(rotation = 45)
	plt.savefig(os.path.join(prefix, 'dendrogram.jaccard.bestcut.pdf'))
	plt.close()
	agg=AgglomerativeClustering(n_clusters=k, metric='precomputed', linkage='average')
	cluster=agg.fit(mtx)
	groups=set(cluster.labels_)

	#result
	for g in groups:
		
		group=list(np.take(mtx.columns.tolist(),np.where(cluster.labels_ == g))[0])
		results.append(group)

	with open(os.path.join(prefix, 'dendrogram.jaccard.bestcut.tsv'), 'w') as outfile:

		for g in result:

			outfile.write(','.join(r) + '\n')	

			for l in g:

				d[l]=[x for x in g if x != l]


	with open(os.path.join(prefix, 'dendrogram.jaccard.bestcut.json'), 'w') as outfile:
		
		json.dump(d, outfile, indent = 4)

def LinkageMatrix(model):

	'''
	Create scipy linkage matrix from sklearn model
	'''

	counts = np.zeros(model.children_.shape[0])
	n_samples = len(model.labels_)

	for i, merge in enumerate(model.children_):

		current_count = 0

		for child_idx in merge:

			if child_idx < n_samples:

				current_count += 1  # leaf node

			else:

				current_count += counts[child_idx - n_samples]

		counts[i] = current_count

	linkage_matrix = np.column_stack([model.children_, model.distances_,counts]).astype(float)

	return linkage_matrix


if __name__ == '__main__':

    mtx=read_tsv(sys.argv[1])
    agglomerative(mtx,sys.argv[2])