import os
import csv
from sys import argv

#open file

print('Opening...')
filename = input('Enter file name with fq extension to analyze >: ')
handle = open(os.path.expanduser(f"~/Desktop/{filename}"))
print('File opened')

w = csv.writer(open('output.csv', 'w'))


#change number for analyzing up to this # of reads

readstocount = int(input('Enter a max number of reads to count >: '))




def counter(num,fin):
    trimmed = list()
    trimming_count = 0
    counts = dict()
    handle.seek(0)
    for line in handle:
        line.rstrip()
        if line.startswith('SRN'):
            trimmed.append(line[num:fin])
            trimming_count += 1
            if trimming_count >= readstocount:
                break
            else:
                continue
    for barcode in trimmed:
            counts[barcode] = counts.get(barcode, 0) + 1
    print()
    print('***********')
    print('Position:', num,'-',fin)
    print(counts)
    w.writerow(['AAV2 pos', num + 474])
    for key, val in counts.items():
        if key != '*' and val > 50:
            w.writerow([key, val])
    print()

for x in range(1, 144, 1):
    counter(x, (x+1))



print('DONE')


#this loop generates an histogram of the frequency of each barcode
