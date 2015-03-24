#!/usr/bin/env python

import os
import sys
import argparse
import subprocess
import struct
import shutil
import tempfile
import csv

from biolib import anabl_getContigsFromFASTA as readFASTA, translateCDS

"""
0: #ID
1: NN_Cmax_score
2: NN_Cmax_pos
3: NN_Cmax_pred
4: NN_Ymax_score
5: NN_Ymax_pos
6: NN_Ymax_pred
7: NN_Smax_score
8: NN_Smax_pos
9: NN_Smax_pred
10: NN_Smean_score
11: NN_Smean_pred
12: NN_D_score
13: NN_D_pred
14: HMM_type
15: HMM_Cmax_score
16: HMM_Cmax_pos
17: HMM_Cmax_pred
18: HMM_Sprob_score
19: HMM_Sprob_pred
"""
PRED_HEADER = ['#SeqID', 'NA_Seq', 'Pep_Seq', 'seqlen_NA', 'seqlen_Pep',
               'Pos_before_CleavageSite(NN)',
               'Pos_after_CleavageSite(NN)',
               'seqlen_Mature',
               'CleavageSite(NN)',
               'NA_Seq_w/o_SigP(NN)',
               'Signal peptide',
               'Mature peptide']



def doStuff(args):
    secreted = []
    for row in csv.reader(open(args.spoutFile), delimiter='\t', quotechar='"'):
        if row[0].startswith('#'):
            continue
        # for i, col in enumerate(row): out.write('%i: %s\n' % (i, col))
        if (float(row[7]) > args.Smax) or (float(row[4]) > args.Ymax) or (float(row[10]) > args.Smean):
            secreted.append([row[0], None, None, None, None, int(row[5]) - 1, int(row[5]), None, None, None, None, None])
    pepDict = dict(readFASTA(args.aaseqFile))
    for secPept in secreted:
        secPept[2] = pepDict.get(secPept[0], '')
        secPept[4] = len(secPept[2])
        secPept[7] = secPept[4] - secPept[5]
        secPept[8] = ''
        secPept[10] = secPept[2][:secPept[5]]
        secPept[11] = secPept[2][secPept[5]:]
    del pepDict

    with open(args.outputFile, 'wb') as out:
        for secPept in secreted:
            out.write('\n'.join(map(str, secPept)) + '\n')
        pass
    
    pass

def main(argv):
    
    descr = ''
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('--Cmax', default=0.32)
    parser.add_argument('--Ymax', default=0.33)
    parser.add_argument('--Smax', default=0.87)
    parser.add_argument('--Smean', default=0.48)
    parser.add_argument('--Dmax', default=0.43)


    parser.add_argument('spoutFile', type=str, help='.')
    parser.add_argument('naseqFile', type=str, help='.')
    parser.add_argument('aaseqFile', type=str, help='.')
    parser.add_argument('outputFile', type=str, help='.')
    
    try:
        args = parser.parse_args()
    except:
        sys.exit(1)

    if not os.path.exists(args.spoutFile):
        sys.stderr.write('Input file (%s) is missing.\n' % args.spoutFile)
        sys.exit(1)
    if not os.path.exists(args.naseqFile):
        sys.stderr.write('Input file (%s) is missing.\n' % args.naseqFile)
        sys.exit(1)
    if not os.path.exists(args.aaseqFile):
        sys.stderr.write('Input file (%s) is missing.\n' % args.aaseqFile)
        sys.exit(1)

    doStuff(args)
    # outargs = [args.Cmax, args.Ymax, args.Smax, args.Smean, args.Dmax]
    # open(args.outputFile, 'wb').write(','.join(map(str, outargs)))



    pass

if __name__ == '__main__': main(sys.argv[1:])
