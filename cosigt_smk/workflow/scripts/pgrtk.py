#!/usr/bin/python3

import pgrtk
import os
import sys

def main(agc_file,bed_file,padding,out_file):

    #load acg_file
    ref_db=pgrtk.AGCFile(agc_file)

    #load region
    with open(bed_file) as bed:

        for line in bed:

            fields=line.rstrip().split('\t')
            ref_file_name, roi_chr, roi_b, roi_e = fields[0],fields[1], int(fields[2]), int(fields[3])
            roi_len=roi_e - roi_b
            break 
            #not needed - but make sure we read just one line

    #sdb
    sdb = pgrtk.SeqIndexDB()
    sdb.load_from_agc_index(agc_file.replace('.agc', ''))

    #retrieve roi sequence
    roi_seq = ref_db.get_sub_seq(ref_file_name, roi_chr, roi_b-padding, roi_e+padding)

    #find hits in the pangenomic reference
    aln_range = pgrtk.query_sdb(sdb, roi_seq, merge_range_tol=100000) #force constant merging
    seq_list=[]

    for k in list(aln_range.keys()):

        ctg_name, source, _ = sdb.seq_info[k]
        rgns = aln_range[k].copy()

        rgns = pgrtk.merge_regions(rgns, tol=1000)

        for rgn in rgns:

            b, e,_, orientation, aln = rgn
            aln.sort()
    
            if aln[0][0][0] > padding or aln[-1][0][1] < padding + roi_len:
                continue
            
            if e-b < 0.75 * (roi_len + 2 * padding):             
                continue

            seq =  sdb.get_sub_seq(source, ctg_name, b, e)
    
            if orientation == 1:

                seq = pgrtk.rc_byte_seq(seq)

            seq_list.append(('{}_{}_{}_{}'.format(ctg_name, b, e, orientation), seq))

    f0 = open(out_file, 'w')
    
    for ctg, seq in seq_list:
            
        print('>{}'.format(ctg), file=f0)
        print(pgrtk.u8_to_string(seq), file=f0)
            
    f0.close()


if __name__ == '__main__':

    agc_file=os.path.abspath(sys.argv[1])
    bed_file=os.path.abspath(sys.argv[2])
    padding=int(sys.argv[3])
    out_file=os.path.abspath(sys.argv[4])
    main(agc_file,bed_file,padding,out_file)



