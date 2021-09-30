#!/usr/bin/python

"""
The barcode file was not found!

Give a valid flowcell and lane specification.

Check the name given to the barcode file is in accord with the name of sequence file.

Usage:
    ./makeBarcodeSabre_V2.py <flowcell> <lane> <seqtype>
"""

import os, sys, fileinput

try:
	flow=sys.argv[1]
	lane=sys.argv[2]
	f = open('barcodes/barcodes_'+flow+'_'+lane,'r')
except:
	print(__doc__)
	sys.exit(1)

try:
	seqtype=sys.argv[3]
except:
	print('Specify a sequence type: SE or PE')
	sys.exit(1)

if seqtype=='SE':
	o = open('data/'+flow+'_'+lane+'_SE','w')
	line = f.readline()
	while line:
		list = line.split()
		o.write(list[0]+' '+flow+'_'+lane+'_'+list[1]+'.fq\n')
		line = f.readline()
	f.close()
	o.close()

if seqtype=='PE':
	o = open('data/'+flow+'_'+lane+'_PE','w')
	line = f.readline()
	while line:
		list = line.split()
		o.write(list[0]+' '+flow+'_'+lane+'_'+list[1]+'_R1.fq'+' '+flow+'_'+lane+'_'+list[1]+'_R2.fq\n')
		line = f.readline()
	f.close()
	o.close()
