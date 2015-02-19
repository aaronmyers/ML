#!/usr/bin/python

from operator import itemgetter
import sys

current_word = None
current_count = 0
word = None
current_filename = None
dummystr = ''
count = 1

for line in sys.stdin:
	line = line.strip()
	word, filename, freq = line.split('\t',2)

	try:
		count = float(count)
	except ValueError:
		continue

	if current_word == word:
		current_count += count
		dummystr += word +'\t' +filename +'\t'+ freq +'\n' 
	else:
		if current_word:
			dummysome = dummystr.splitlines()
			for line in dummysome:
				print '%s\t%s' % (line,current_count)
		current_count = count
		current_word = word
		dummystr = word +'\t' +filename +'\t'+ freq +'\n'  

if current_word == word:
	dummysome = dummystr.splitlines()
	for line in dummysome:
		print '%s\t%s' % (line,current_count)

