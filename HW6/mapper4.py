#!/usr/bin/python

from operator import itemgetter
import sys
import string
import math
import os

current_word = None
countw = 0
for line in sys.stdin:
	line = line.strip()
	word, filename, tfidf = line.split('\t',2)
	
	if current_word == word:
		print '%s\t%s\t%s' % (filename, countw, tfidf)
	else:		
		#Add the iteration through the dictionary for the test data
		countw += 1
		print '%s\t%s\t%s' % (filename, countw, tfidf)
		current_word = word	




	
