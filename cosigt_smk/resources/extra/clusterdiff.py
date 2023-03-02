
#!/usr/bin/python3


import os
import sys
import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from matplotlib import pyplot as plt
from scipy.cluster import hierarchy
import matplotlib.cm as cm
from scipy.cluster import hierarchy
from sklearn.metrics import silhouette_samples, silhouette_score


def read_tsv(tsv):

	'''
	Read tsv with distances between haplotypes
	'''

	return 1-pd.read_table(tsv).set_index(['path.a', 'path.b'])['jaccard'].sort_index().unstack()


def agglomerative(mtx,n):

	'''
	Perform agglomerative clustering and plot dendrogram
	'''

	result=[]

	agg=AgglomerativeClustering(distance_threshold=0, n_clusters=None, affinity='precomputed', linkage='average')
	cluster=agg.fit(mtx)

	linkage_matrix=LinkageMatrix(agg)

	plt.figure()
	dendrogram = hierarchy.dendrogram(linkage_matrix)
	plt.savefig("dendrogram.jaccard.nocut.png")
	plt.close()

	range_n_clusters = range(2,n)

	for i,n_clusters in enumerate(range_n_clusters):

		clusterer = AgglomerativeClustering(n_clusters=n_clusters, affinity='precomputed', linkage='average')
		cluster_labels = clusterer.fit_predict(mtx)
		silhouette_avg = silhouette_score(mtx, cluster_labels, metric='precomputed')
		result.append((n_clusters,silhouette_avg))

	k=sorted(result, key=lambda x: x[1], reverse=True)[0][0]

	#best

	plt.figure()
	dendrogram = hierarchy.dendrogram(linkage_matrix, truncate_mode = 'lastp', p = k)
	plt.savefig("dendrogram.jaccard.cut"+ str(k) + ".png")
	plt.close()

	agg=AgglomerativeClustering(n_clusters=k, affinity='precomputed', linkage='average')
	cluster=agg.fit(mtx)
	groups=set(cluster.labels_)

	result=[]

	for g in groups:
		
		group=list(np.take(mtx.columns.tolist(),np.where(cluster.labels_ == g))[0])
		result.append(group)

	d=dict()

	for g in result:

		for l in g:

			d[l]=[x for x in g if x != l]


	with open("dendrogram.jaccard.cut"+ str(k) + ".json", "w") as outfile:
		
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


def main():

	mtx=read_tsv("haplodiff.tsv")
	agglomerative(mtx,n=14)
