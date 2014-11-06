#!/usr/bin/python

import sys

f = open('Data3893.dat','w')
for line in sys.stdin:
	line = line.strip()
	i, j, val = line.split('\t',2)
	st = '%s\t%s\t%s\n' % (i,j,val)
	f.write(st)

f.close()

