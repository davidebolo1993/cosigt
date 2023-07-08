#!/usr/bin/python3

import json
import sys
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def infiles(infolder):

	'''
	find evaluation tables
	'''

	return glob.glob(infolder + '/*/evaluation.tsv')


def read_json(injson):

	'''
	read json with similarities between haplotypes
	'''

	try:

		with open(injson, 'r') as j:

			data=json.load(j)

		return data
	
	except:

		return None



def GetHaplotypes(sample,df):

	'''
	get possible haplotype names
	'''

	out=[]

	for k in set(df[3].tolist() + df[4].tolist()):

		if sample in k:

			out.append(k)

	return out



def eval_sample(infile,json):

	'''
	evaluate
	'''

	df=pd.read_table(infile, sep='\t', header=None)
	sample=list(df.head(1)[2])[0]
	haplos=GetHaplotypes(sample,df)
	vec=np.zeros(df.shape[0])

	h1 = [haplos[0]]

	if json is not None:

		if h1[0] in json:

			h1.extend(json[h1[0]])

	if len(haplos) > 1:

		h2=[haplos[1]]

		if json is not None:

			if h2[0] in json:

				h2.extend(json[h2[0]])

	else: #h2 can be anything

		h2=list(set(df[3].tolist() + df[4].tolist()))

	for i,row in df.iterrows():

		if (row[3] in h1 and row[4] in h2) or (row[3] in h2 and row[4] in h1):
		   
			vec[i:]=1
			break

		else:

			vec[i] = 0

	return (sample, vec)
	   

def main():

	'''
	calculate true positive rate
	'''

	eval_files=infiles(os.path.abspath(sys.argv[1]))
	jsonfile=read_json(os.path.abspath(sys.argv[2]))   
	vectors,samples=[],[]

	for f in eval_files:

		sample,eval_vec=eval_sample(f,jsonfile)
		vectors.append(eval_vec)
		samples.append(sample)

	result=np.zeros(len(vectors[0]))

	with open(os.path.abspath(sys.argv[3]).replace('.pdf', '.badsamples.tsv'), 'w') as outreport:

		for i in range(len(vectors[0])):

			nogood=[]
			sum=0

			for l,el in enumerate(vectors): #for each combination, store results plus names of no_good samples

				if el[i] != 1:

					nogood.append(samples[l])

				sum+=el[i]

			result[i] = sum/len(vectors)
			outreport.write(str(i+1) + '\t' + ','.join(nogood) + '\n')

	if jsonfile is not None:

		cut=os.path.basename(sys.argv[2]).split('.')[1]

	else:

		cut = 'no'

	plt.title('Evaluation - ' + cut + ' dendrogram')
	plt.ylim([0, 1])
	plt.plot(list(range(len(result))), result, marker='*')
	plt.xlabel('Combination index')
	plt.ylabel('True positive rate')
	plt.savefig(os.path.abspath(sys.argv[3]))



if __name__ == '__main__':

	main()
