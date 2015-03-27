#!/usr/bin/env python

import os 
import sys
import csv
import argparse


def doStuff(args):	
	def removeRedundantCells(row):
		del row[6]
		del row[1]
	uniquePeptides = set()
	with open(args.filteredFile, 'wb') as out:
		reader = csv.reader(open(args.spoutFile), delimiter='\t', quotechar='"')
		for row in reader:
			if row[0].startswith('#'):
				removeRedundantCells(row)
				out.write('\t'.join(row) + '\n')
				continue
			if row[2] in uniquePeptides:
				continue
			uniquePeptides.add(row[2])
			if len(row) == 12:				
				matureLength, signalLength = map(int, [row[7], row[5]])
			else:
				matureLength, signalLength = len(row[2]) - int(row[5]), int(row[5])

			if matureLength >= args.minMatureLength and signalLength >= args.minSignalLength:				
				removeRedundantCells(row)
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