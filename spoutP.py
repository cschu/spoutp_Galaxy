#!/usr/bin/env python

import os
import sys
import argparse
import subprocess
import struct
import shutil
import tempfile

from biolib import anabl_getContigsFromFASTA as readFASTA, translateCDS


def doStuff(fi, fo):
    with open(fo, 'wb') as out:
        for id_, seq in readFASTA(fi):
            # here take only peptides into account 
            # that start with MET
            # if seq[:3].upper().replace('U', 'T') == 'ATG':
            # that is actually done by translateCDS
            translated = translateCDS(seq)
            if translated:
                out.write('>%s\n%s\n' % (id_, translated))
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

    if not os.path.exists(args.inputFile):
        sys.stderr.write('Input file (%s) is missing.\n' % args.inputFile)
        sys.exit(1)

    # doStuff(args.inputFile, args.outputFile)
    open(args.outputFile, 'wb').write(','.join(map(float, [args.Cmax, args.Ymax, args.Smax, args.Smean, args.Dmax])))
    


    pass

if __name__ == '__main__': main(sys.argv[1:])
