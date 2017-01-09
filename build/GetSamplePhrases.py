import sys
import re

if len(sys.argv) != 3:
	print 'Use: python ./GetSamplePhrases.py filename  output_filename'
	sys.exit('Error')

arg1 = sys.argv[1]
arg2 = sys.argv[2]
input_f = open(arg1,'r')
output_f = open(arg2, 'w')

length = 100
sample = 0

text = input_f.read()

n = len(text)
block  = n / 10000



while (sample < 10000):
    pos = sample * block
    s_phrase = text[pos:pos+length]
    s_phrase += '\0'
    output_f.write(s_phrase)
    sample += 1
	

	





