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
    uniquePeptides = set()
    with open(fo, 'wb') as out:
        for id_, seq in readFASTA(fi):
            # here take only peptides into account 
            # that start with MET
            # if seq[:3].upper().replace('U', 'T') == 'ATG':
            # that is actually done by translateCDS
            translated = translateCDS(seq)            
            if translated:
                translated = translated.split('*')[0]   
                if translated not in uniquePeptides:
                    uniquePeptides.add(translated)
                    out.write('>%s\n%s\n' % (id_, translated))
    pass

def main(argv):
    
    descr = ''
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('inputFile', type=str, help='The nucleic acid input file.')
    parser.add_argument('outputFile', type=str, help='The translated coding sequences.')
    
    try:
        args = parser.parse_args()
    except:
        sys.exit(1)

    if not os.path.exists(args.inputFile):
        sys.stderr.write('Input file (%s) is missing.\n' % args.inputFile)
        sys.exit(1)

    doStuff(args.inputFile, args.outputFile)


    pass

if __name__ == '__main__': main(sys.argv[1:])
