#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
import itertools
from scipy import spatial
from operator import itemgetter


def main(bin_tsv,sam_tsv,outf):

    #load table
    binary_df=pd.read_table(bin_tsv).set_index("path.name", drop=True)
    #convert to numpy
    binary_ar=binary_df.to_numpy()
    #get haplotypes
    haplos=binary_df.index.tolist()
    #generate combinations indexes
    combos=list(itertools.combinations(enumerate(haplos),2))
    #initialize empty array
    #with no normalization, ok to use int
    matrix = np.zeros(shape=(len(combos),binary_df.shape[1]),dtype=int)
    # recode zeros in input vectors to this value
    zero_weight=-1e3
    #with normalization, we have to use floats
    #matrix = np.zeros(shape=(len(combos),binary_df.shape[1]),dtype=np.float32)

    #iterate over the combinations and get the sum of corresponding values
    for i,combo in enumerate(combos):
        matrix[i]=binary_ar[combo[0][0]] + binary_ar[combo[1][0]]
        matrix[i]=np.where(matrix[i]==0, zero_weight, matrix[i])
    #store as table
    combo_df = pd.DataFrame(matrix,columns=binary_df.columns.tolist())
    combo_df.index = [combo[0][1] + "/" +combo[1][1] for combo in combos]
    #combo_pairs = [[combo[0][1], combo[1][1]] for combo in combos]
    combo_df.to_csv(os.path.abspath(outf+ "/combo.tsv"), sep = "\t", index=True)
    #do the dot product
    cov_tsv=pd.read_table(sam_tsv).set_index("#sample", drop=True)
    
    #cosine similarity
    results=list()
    count=0
    matrix_cov=np.zeros(shape=(len(cov_tsv.index),len(combos)),dtype=np.float32)
    for idx,row in cov_tsv.iterrows():
        
        s_results=list()
        # zero weighting ~ -(max possible read coverage)
        row=np.where(row==0, zero_weight, row)

        for idx2,row2 in combo_df.iterrows():
            cos_sim=1-spatial.distance.cosine(row.tolist(), row2.tolist())
            s_results.append((idx,idx2,cos_sim))

        matrix_cov[count] = [x[2] for x in s_results]
        results.append(sorted(s_results,key=itemgetter(2),reverse=True)[0])
        #results.append(sorted(s_results,key=itemgetter(2),reverse=False)[0])
        count+=1
        
    #convert to df
    dot_df=pd.DataFrame(matrix_cov,columns=combo_df.index.tolist())
    dot_df.index=cov_tsv.index
    dot_df.to_csv(os.path.abspath(outf+"/genotype.tsv"), sep = "\t", index=True)
    
    #store best score and combo. If multiple highest scores, keep best
    #for dot product
    #perf_df=pd.concat([dot_df.max(axis=1), dot_df.idxmax(axis=1)], axis=1)
    #for cosine similarity
    perf_df=pd.DataFrame(results, columns=['#sample', 'best_genotype', 'best_score'])
    #for dot product
    #perf_df.to_csv(os.path.abspath(outf+ "/best_genotype.tsv"), sep = "\t", index=True)
    #for cosine similarity
    perf_df.to_csv(os.path.abspath(outf+ "/best_genotype.tsv"), sep = "\t", index=False)


if __name__ == '__main__':

    bin_tsv=os.path.abspath(sys.argv[1])
    sam_tsv=os.path.abspath(sys.argv[2])
    outfolder=os.path.abspath(sys.argv[3])
    os.makedirs(outfolder,exist_ok = True)
    main(bin_tsv,sam_tsv,outfolder)
