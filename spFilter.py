#!/usr/bin/env python

import os 
import sys
import csv
import argparse


def doStuff(args):	
	#print "zz"
	with open(args.filteredFile, 'wb') as out:
		reader = csv.reader(open(args.spoutFile), delimiter='\t', quotechar='"')
		#out.write('X\n')
		#out.write('\t'.join(reader.next()) + '\n')
		for row in reader:
			if row[0].startswith('#'):
				out.write('\t'.join(row) + '\n')
				continue
			# out.write(str(len(row)) + '\n')
			if len(row) == 12:				
				matureLength, signalLength = map(int, [row[7], row[5]])
			else:
				matureLength, signalLength = len(row[2]) - int(row[5]), int(row[5])
				
			if matureLength >= args.minMatureLength and signalLength >= args.minSignalLength:
				out.write('\t'.join(row) + '\n')
		pass
	pass



def main(argv):
	
	descr = ''
	parser = argparse.ArgumentParser(description=descr)
	parser.add_argument('--minSignalLength', type=int, default=5)
	parser.add_argument('--minMatureLength', type=int, default=10)

	parser.add_argument('spoutFile', type=str, help='.')
	parser.add_argument('filteredFile', type=str, help='.')

	try:
		args = parser.parse_args()
	except:
		# print "xx"
		sys.exit(1)

	if not os.path.exists(args.spoutFile):
		sys.stderr.write('Input file (%s) is missing.\n' % args.spoutFile)
		sys.exit(1)

	doStuff(args)


if __name__ == '__main__': main(sys.argv[1:])