#!/usr/bin/python
# args = fastq gzip files

import sys
from gzipfile import *

gz = GzipFile(sys.argv[1])
nl = gz.namelist()
f = [fn for fn in nl if "summary.txt" in fn]

if(len(f) != 1):
	print "Error: did not find summary.txt"
	exit(1)
	
su = gz.read(f[0]).splitlines()
gz.close()

if(len(sys.argv) > 1):
	header = "File"
	for i in su:
		header += "\t" + i.split("\t")[1]
	print header



for gzipfile in sys.argv[1:]:
	try:
		gz = GzipFile(zipfile)
	except:
		print "Problem with:",gzipfile
		raise
	nl = gz.namelist()
	f = [fn for fn in nl if "summary.txt" in fn]
	
	if(len(f) != 1):
		print "Error: did not find summary.txt"
		exit(1)
		
	su = z.read(f[0]).splitlines()
	z.close()
	out = su[0].split("\t")[2]
	for i in su:
		out += "\t" + i.split("\t")[0]
	print out
	
	
	il=@K00171:178:H7TL5BBXX:7:1101:25337:2422 1:N:0:CGACACAC
	print il.split(':', 1)
	@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1
ID=
SM=
PU=
PL=
LB=
split(3: )
