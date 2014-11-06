#!/usr/bin/python

import sys
import string

current_word = None
countw=0
f = open('Dict2.txt','w')

for line in sys.stdin:
	line = line.strip()
	word, filename, id = line.split()
	if current_word != word and len(word)<20:
		countw +=1
		#print to file	word, countw
		st =  '%s\t%s\n' % (word.lower(), countw)
		f.write(st)
		current_word = word

f.close()
