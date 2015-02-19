#!/usr/bin/python

from operator import itemgetter
import sys
import string
import math
import os

current_filename = None
dummy = ''

for line in sys.stdin:
	line = line.strip()
	filename, word, tfidf = line.split('\t',2)
	
	if current_filename == filename:
		dummy += '\t' + word + ':' + tfidf		
	else:		
		if current_filename:
			print '%s' % (dummy)
		current_filename = filename	
		derp = filename.split('/')
		dummy = derp[-2]
		if dummy=='pos':
			dummy = '+1' + '\t' + word + ':' + tfidf
		else:
			dummy = '-1' + '\t' + word + ':' + tfidf

if current_filename:
	print '%s' % (dummy)
	
