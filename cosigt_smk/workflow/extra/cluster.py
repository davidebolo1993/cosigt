#!/usr/bin/python3

import editdistance
from Bio import Phylo
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.cluster import DBSCAN
import numpy as np
import json 
import sys
import os

def editdist(x,y):

	'''
	Custom function for string clustering
	'''

	return int(editdistance.eval(data[int(x[0])], data[int(y[0])]))



def diffperc(x,y):

	'''
	Custom function for string clustering
	'''

	return 100*int(editdistance.eval(data[int(x[0])], data[int(y[0])]))/max(len(data[int(x[0])]),len(data[int(y[0])]))


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

		next(f)
		for line in f:

			key, value=line.split()
			
			if value in d:

				d[value].append(key)

			else:

				d[value] = [key]

	return d


def cluster_tree(tree,affinity):

	'''
	Cluster clades
	'''
	#output
	result=[]
	
	#input, make it global
	global data
	data=[]

	for clade in tree.get_terminals():
	
		key = clade.name
		data.append(key)

	X = np.arange(len(data)).reshape(-1, 1)
	metric= pairwise_distances(X, X, metric=diffperc)
	agg=DBSCAN(eps=100-affinity, min_samples=1,algorithm='brute', metric='precomputed')
	cluster_=agg.fit(metric)
	groups=set(cluster_.labels_)

	for g in groups:
		
		group=list(np.take(data,np.where(cluster_.labels_ == g))[0])
		result.append(group)

	return result


""" def label_dendrogram(tsv_dict,tree):

	'''
	Assign labels to original dendrogram
	'''

	d={}

	for clade in tree.get_terminals():
	
		key = clade.name
		d[key]=tsv_dict[key]

	return d
 """



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

	

def main():

	'''
	Run
	'''

	nwk=os.path.abspath(sys.argv[1])
	tsv=os.path.abspath(sys.argv[2])
	outpath=os.path.abspath(sys.argv[3])

	tree=read_newick(nwk)
	dendro_labels=read_tsv(tsv)
	outcluster(dendro_labels,outpath + "/dendrogram.clusters.tsv")
	prox_dendr=proximities(dendro_labels)

	with open(outpath + "/dendro.proxi.json", "w") as outfile:
		
		json.dump(prox_dendr, outfile,indent = 4)

	new_clust=cluster_tree(tree,90.0)
	cluster_labels=label_clustering_results(dendro_labels,new_clust)
	outcluster(cluster_labels, outpath + "/DBSCAN.clusters.tsv")
	prox_clust=proximities(cluster_labels)

	with open(outpath + "/DBSCAN.proxi.json", "w") as outfile:
		
		json.dump(prox_clust, outfile, indent = 4)



if __name__ == '__main__':

	main()
