import os
from sys import argv

#open file

print('Opening...')
filename = input('Enter file name with fq extension to analyze >: ')
handle = open(os.path.expanduser(f"~/Desktop/{filename}"))
print('File opened')




#change number for analyzing up to this # of reads

readstocount = int(input('Enter a max number of reads to count >: '))

#this loop trims the 72 first nucleotides of the sequence


def counter(pos1,pos1e):
    trimmed = list()
    trimming_count = 0
    counts = dict()
    handle.seek(0)
    for line in handle:
        line.rstrip()
        if line.startswith('CCTTT'):
            trimmed.append(line[pos1:pos1e])
            trimming_count += 1
            if trimming_count >= readstocount:
                break
            else:
                continue
    for barcode in trimmed:
            counts[barcode] = counts.get(barcode, 0) + 1

    for k,v in counts.items():
        if v >= 700:
            print(k)
            print(v)
    print()
    print('***********')
    print('Position:', pos1,'-',pos1e)
    print(counts)
    print()


counter(81, 85)



print('DONE')


#this loop generates an histogram of the frequency of each barcode
