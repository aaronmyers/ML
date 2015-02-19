#!/usr/bin/python

from operator import itemgetter
import sys
import string
import math
import os

current_word = None
fcount = 0
f = open('/home/aaron/Dict2.txt','r')

for line in sys.stdin:
	line = line.strip()
	word, filename, tfidf = line.split('\t',2)
	
	if current_word == word:
		print '%s\t%s\t%s' % (filename, fcount, tfidf)
	else:		
		#Add the iteration through the dictionary for the test data
		k=0
		while k < 1:
			p = f.tell()
			line2 = f.readline()
			if not line2:
				break
			line2 = line2.strip()
			fword, fcount = line2.split('\t',1)
			if word >= fword:
				if word == fword:
					k=1
					print '%s\t%s\t%s' % (filename, fcount, tfidf)
					current_word = word	
			else:
				k = 1
				f.seek(p)


f.close()


	
