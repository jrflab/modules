#!/usr/bin/env python
# Description: This script tries to align the read names of paired fastq files.  The fastq extraction programs do not check of this and often reads names are not matching properly.  Some aligners like gsnap complain about this refusing to run.  This script aims to generate a common prefix name for both paired-reads and then appends a /1 or /2 to the paired fastq respectively
# Authors: Fong Chun Chan <fongchunchan@gmail.com>
import argparse
import re
import sys

parser = argparse.ArgumentParser(description='Corrects the reads in the paired fastq files so that they match', epilog='By Fong Chun Chan')
parser.add_argument('fastqFile1', action='store', help='Fastq file 1')
parser.add_argument('fastqFile2', action='store', help='Fastq file 2')
parser.add_argument('fixedFastqFile1', action='store', help='Fixed fastq file 1')
parser.add_argument('fixedFastqFile2', action='store', help='Fixed fastq file 2')
args = parser.parse_args()

print 'Paired fastq files: ', args.fastqFile1, ' - ', args.fastqFile2
inFile1 = open(args.fastqFile1, 'rb')
inFile2 = open(args.fastqFile2, 'rb')

print 'Output to fixed paired fastq files: ', args.fixedFastqFile1, ' - ', args.fixedFastqFile2
outFile1 = open(args.fixedFastqFile1, 'wb')
outFile2 = open(args.fixedFastqFile2, 'wb')

i = 1
while True:
	line1 = inFile1.readline()
	line2 = inFile2.readline()

	if line1 == '' and line2 == '':
		print 'End of files reached.  Successfully fixed fastq read names'
		sys.exit(0)

	elif line1 == '' and line2 != '':
		print 'File 1 has reached the end before File 2'
		sys.exit(1)

	elif line1 != '' and line2 == '':
		print 'File 2 has reached the end before File 1'
		sys.exit(1)

	elif i == 1 or ( i % 4 ) == 1:
		try: 
			m = re.search('(.+/)\d', line1)
			readName = m.group(1)
		except AttributeError as e:
			readName = line1.rstrip() + '/'
		outFile1.write(readName + '1' + '\n')
		outFile2.write(readName + '2' + '\n')
	else : 
		outFile1.write(line1)
		outFile2.write(line2)
	i = i + 1
