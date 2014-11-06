#!/usr/bin/python

from operator import itemgetter
import sys

current_word = None
current_count = 0
word = None
current_filename = None

for line in sys.stdin:
	line = line.strip()

	word, filename, count = line.split('\t', 2)

	try:
		count = float(count)
	except ValueError:
		continue

	if current_word == word and current_filename == filename:
		current_count += count
	else:
		if current_word:
			print '%s\t%s\t%f' % (current_word, current_filename, current_count)
		current_count = count
		current_word = word
		current_filename = filename

if current_word == word:
	print '%s\t%s\t%f' % (current_word, current_filename, current_count)

