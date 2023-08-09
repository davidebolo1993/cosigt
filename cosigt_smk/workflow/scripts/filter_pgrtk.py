#!/usr/bin/python3 env

#standard libraries
import os
import sys
import pyfaidx

def main(input_fa,input_hit,output_fa):

    fasta=pyfaidx.Fasta(input_fa)

    with open(input_hit) as fin, open(output_fa, 'w') as fout:
        
        for row in fin:
            
            row = row.strip().split("\t")
            
            if row[0][0] == "#":
                
                continue
            
            (q_idx, q_name, query_bgn, query_end, q_len, aln_anchor_count, 
            src, ctg, ctg_bgn, ctg_end, orientation, out_seq_name) = row
            query_bgn = int(query_bgn)
            query_end = int(query_end)
            ctg_bgn = int(ctg_bgn)
            ctg_end = int(ctg_end)
            q_len = int(q_len)
            
            #maybe need to check why this is 12000 - but itworks for test regions
            if query_bgn > 12000 or q_len - query_end > 12000:
            
                continue
            
            if abs(ctg_end-ctg_bgn) > q_len * 3:
            
                continue
            
            key=[x for x in fasta.keys() if ctg+ '_' + str(ctg_bgn) in x]
            seq=fasta[key[0]][:].seq
            fout.write('>' + key[0].split('::')[1] + '\n' + seq + '\n')


if __name__ == '__main__':

    input_fa=os.path.abspath(sys.argv[1])
    input_hit=os.path.abspath(sys.argv[2])
    output_fa=os.path.abspath(sys.argv[3])
    main(input_fa,input_hit,output_fa)