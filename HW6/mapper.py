#!/usr/bin/python
#lkdf
import os
import sys
import string
filename = os.environ.get('mapreduce_map_input_file')
count = 0
dummystr = ''

for line in sys.stdin:
	line = line.strip()
	line = line.translate(None,'!@#$%^&*()\'.",/?{}\|;:=+-`~]><[_')
	words = line.split()
	for word in words:	
		dummystr += word + '\n'
		count+=1
		
# printing out the total word count
count = float(count)
dummysome = dummystr.splitlines()
for line in dummysome:
	print '%s\t%s\t%f' % (line, filename, 1/count)

