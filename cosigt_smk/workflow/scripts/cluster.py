#!/usr/bin/python3

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
	read tsv with distances between haplotypes from odgi similarity
	'''

	return 1-pd.read_table(tsv).set_index(['group.a', 'group.b'])['jaccard.similarity'].sort_index().unstack()

def agglomerative(mtx,prefix):

	'''
	perform agglomerative clustering and plot dendrogram
	'''

	agg=AgglomerativeClustering(distance_threshold=0, n_clusters=None, metric='precomputed', linkage='average')
	cluster=agg.fit(mtx)
	linkage_matrix=LinkageMatrix(agg)

	sil_score=[] #all silhouette scores

	range_n_clusters = range(2,mtx.shape[0])

	for i,n_clusters in enumerate(range_n_clusters):

		clusterer = AgglomerativeClustering(n_clusters=n_clusters, metric='precomputed', linkage='average')
		cluster_labels = clusterer.fit_predict(mtx)
		silhouette_avg = silhouette_score(mtx, cluster_labels, metric='precomputed')
		sil_score.append((n_clusters,silhouette_avg))

	#get best
	k=sorted(sil_score, key=lambda x: x[1], reverse=True)[0][0]

	#results
	hapdict=dict()
	agg=AgglomerativeClustering(n_clusters=k, metric='precomputed', linkage='average')
	cluster=agg.fit(mtx)
	groups=set(cluster.labels_)
	haplotype_labels=list()

	with open(prefix + '.clusters.tsv', 'w') as outfile:

		outfile.write('group\thaplotypes\n')	

		for i,g in enumerate(groups):
			
			group=list(np.take(mtx.columns.tolist(),np.where(cluster.labels_ == g))[0])

			for h in group:

				hapdict[h] = 'HG'+str(i)

			outfile.write('HG'+str(i) + '\t' + ','.join(group) + '\n')	
			haplotype_labels.append('HG'+str(i) + ' (' + str(len(group)) + ')')

	with open(prefix + '.clusters.json', 'w') as outfile:
		
		json.dump(hapdict, outfile, indent = 4)

	dendrogram = hierarchy.dendrogram(linkage_matrix, truncate_mode = 'lastp', p = k)
	hgl = {dendrogram['leaves'][ii]: haplotype_labels[ii] for ii in range(len(dendrogram['leaves']))}

	#plot best
	plt.figure(figsize=(10,6))
	hierarchy.dendrogram(
            linkage_matrix,
            truncate_mode='lastp',  # show only the last p merged clusters
            p=k,  # show only the last p merged clusters
            leaf_label_func=lambda x: hgl[x],
            leaf_rotation=45,
            leaf_font_size=5,
            show_contracted=True 
            )
	plt.xlabel('haplotype group (# haplotypes)')
	plt.ylabel('jaccard dissimilarity')
	plt.xticks(rotation = 45)
	plt.savefig(prefix + '.clusters.pdf')
	plt.close()



def LinkageMatrix(model):

	'''
	create scipy linkage matrix from sklearn model
	'''

	counts = np.zeros(model.children_.shape[0])
	n_samples = len(model.labels_)

	for i, merge in enumerate(model.children_):

		current_count = 0

		for child_idx in merge:

			if child_idx < n_samples:

				current_count += 1  #leaf node

			else:

				current_count += counts[child_idx - n_samples]

		counts[i] = current_count

	linkage_matrix = np.column_stack([model.children_, model.distances_,counts]).astype(float)

	return linkage_matrix


if __name__ == '__main__':

    mtx=read_tsv(sys.argv[1])
    agglomerative(mtx,sys.argv[2])