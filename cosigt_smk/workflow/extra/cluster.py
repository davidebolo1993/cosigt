#!/usr/bin/python3

import numpy as np
import json 
import sys
import os
import editdistance
from Bio import Phylo
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.cluster import AgglomerativeClustering,DBSCAN
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram,cut_tree
import matplotlib.pyplot as plt



def read_newick(nwk):

	'''
	Read tree in newick format
	'''

	tree = Phylo.read(nwk, "newick")

	return tree


def read_tsv(tsv):

	'''
	Read tsv with structure-to-name conversion
	'''

	d={}

	with open(tsv, 'r') as f:

		next(f) #skip header
		for line in f:

			key, value=line.split()
			
			if value in d:

				d[value].append(key)

			else:

				d[value] = [key]

	return d

def label_clustering_results(tsv_dict,cluster_results):

	'''
	Assign sequences to haplotypes
	'''

	res_d={}

	for i,groups in enumerate(cluster_results):

		labels=[]

		for g in groups:

			labels.extend(tsv_dict[g])

		res_d[i]=labels

	return res_d



def outcluster(cluster, outfile):

	'''
	Output file
	'''

	with open(outfile,'w') as f:

		for k in cluster:

			f.write(",".join(cluster[k]) + "\n")



def proximities(cluster):

	'''
	Proximities - from the same cluster
	'''

	d={}

	for k in cluster:

		ids=cluster[k]

		for id in ids:

			d[id] = [x for x in ids if x != id]

	return d


def TreeToLinkage(tree):
	
	'''
	Convert tree to Linkage matrix
	'''

	tree_height = max(tree.distance(c) for c in tree.find_clades(terminal=True))

	# add comment with id and span of terminal nodes:
	id_map = {}

	for i, c in enumerate(tree.find_clades(terminal=True)):

		c.comment = (i, 1) 
		id_map[i] = c.name

	# ancestor list orderred by distance from the root

	anc_lst = []

	for c in tree.find_clades(terminal=False):

		d = tree.distance(c)
		anc_lst.append((c, list(c), d))

	anc_lst.sort(key=lambda x:x[2], reverse=True)

	# running number of node

	nodes = len(list(tree.find_clades(terminal=True)))
	lnk_lst = []
	for anc,children, anc_d in anc_lst:

		n_children = len(children)
		assert n_children>=2
		child1 = children[0]

		for child2 in children[1:]:

			id1, n_leaves1 = child1.comment
			id2, n_leaves2 = child2.comment
			total_leaves = n_leaves1 + n_leaves2
			anc.comment = (nodes, total_leaves)
			distance = tree_height - anc_d
			nodes += 1
			row = [id1, id2, distance, total_leaves]
			lnk_lst.append(row)
			child1 = anc

	return np.array(lnk_lst), id_map


def GroupBundles(ids_old, idx_new):

	'''
	Re-group bundles
	'''

	m=dict()

	for i,idx in enumerate(idx_new):

		if idx in m:

			m[idx].append(ids_old[i])

		else:

			m[idx] = [ids_old[i]]

	return m


def GroupIds(idx_new):

	'''
	Re-group ids
	'''

	m=dict()

	for i,idx in enumerate(idx_new):

		if idx in m:

			m[idx].append(str(i))

		else:

			m[idx] = [str(i)]

	return m



def Recluster(groups):

	'''
	Re-cluster ids
	'''

	result=[]

	for key in groups:

		group=[]

		for val in groups[key]:

			group.append(val)

		result.append(group)

	return result



def main():

	'''
	Run
	'''

	#deal with input files
	nwk=os.path.abspath(sys.argv[1])
	tsv=os.path.abspath(sys.argv[2])
	outpath=os.path.abspath(sys.argv[3])

	#read newick
	tree=read_newick(nwk)
	dendro_labels=read_tsv(tsv)

	#no cluster
	outcluster(dendro_labels,outpath + "/dendrogram.nocut.tsv")
	prox_dendr=proximities(dendro_labels)

	with open(outpath + "/dendrogram.nocut.json", "w") as outfile:
		
		json.dump(prox_dendr, outfile,indent = 4)

	#re-build the linkage matrix to use with scipy and cut the dendrogram
	linkage,ids=TreeToLinkage(tree)
	plt.figure()
	dendrogram = hierarchy.dendrogram(linkage)
	plt.savefig(outpath + "/dendrogram.nocut.png")
	plt.close()

	#cut at predefined treshold - manual cut
	for d in reversed(range(2,15,1)):

		cut_linkage=cut_tree(linkage, n_clusters=d).reshape(-1,)
		cut_groups=GroupBundles(ids,cut_linkage)
		cut_ids=GroupIds(cut_linkage)

		plt.figure()
		dendrogram = hierarchy.dendrogram(linkage, truncate_mode = 'lastp', p = d)
		plt.savefig(outpath + "/dendrogram.cut"+ str(d) + ".png")
		plt.close()

		with open(outpath + "/dendrogram.cut"+ str(d) + ".ids.tsv", "w") as outfile:

			for ix in cut_ids:
				
				outfile.write(",".join(cut_ids[ix]) + "\n")

		cut_cluster=Recluster(cut_groups)
		cut_cluster_labels=label_clustering_results(dendro_labels,cut_cluster)
		outcluster(cut_cluster_labels, outpath + "/dendrogram.cut"+ str(d) + ".tsv")
		prox_clust=proximities(cut_cluster_labels)

		with open(outpath + "/dendrogram.cut"+ str(d) + ".json", "w") as outfile:
		
			json.dump(prox_clust, outfile, indent = 4)



if __name__ == '__main__':

	main()
