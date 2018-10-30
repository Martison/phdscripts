import os
from sys import argv

#open file

file = open('merged.txt', 'w')

print('Opening...')
handle = open('A2SAT-M_combined_R.fastq', 'r').readlines()
print('File opened')

for i in range (1,10000000000, 8):
    file.write(handle[i].rstrip()+handle[i+4])



print('Done')
