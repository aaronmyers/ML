#!/usr/bin/python
import os
import sys
import string
import math

D = os.environ.get('poop')
D = float(D)

for line in sys.stdin:
	word, filename, freq, d = line.split('\t',3)
	d = float(d)
	freq = float(freq)
	print '%s\t%s\t%s' % (word, filename, freq*math.log(D/d)) 	
