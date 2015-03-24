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


def doStuff(args):
    with open(args.outputFile, 'wb') as out:
        for row in csv.reader(open(args.spoutFile), delimiter='\t', quotechar='"'):
            for i, col in enumerate(row): out.write('%i: %s\n' % (i, col))
    
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
