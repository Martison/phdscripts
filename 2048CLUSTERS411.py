import os
from sys import argv

#open file

print('Opening...')
filename = input('Enter file name with fq extension to analyze >: ')
handle = open(os.path.expanduser(f"~/Desktop/{filename}"))
print('File opened')

codontable = {'AAA': 'F', 'GAA': 'F', 'TAA': 'L', 'CAA': 'L', 'AAG': 'L', 'GAG': 'L', 'TAG': 'L', 'CAG': 'L', 'AAT': 'I', 'GAT': 'I', 'TAT': 'I', 'CAT': 'V', 'AAC': 'V', 'GAC': 'V', 'TAC': 'V', 'AGA': 'S', 'GGA': 'S', 'TGA': 'S', 'CGA': 'S', 'AGG': 'P', 'GGG': 'P', 'TGG': 'P', 'CGG': 'P', 'AGT': 'T', 'GGT': 'T', 'TGT': 'T', 'CGT': 'T', 'AGC': 'A', 'GGC': 'A', 'TGC': 'A', 'CGC': 'A', 'ATA': 'Y', 'GTA': 'Y', 'TTA': '_', 'CTA': '_', 'ATG': 'H', 'GTG': 'H', 'TTG': 'Q', 'CTG': 'Q', 'ATT': 'N', 'GTT': 'N', 'TTT': 'K', 'CTT': 'K', 'ATC': 'D', 'GTC': 'D', 'TTC': 'E', 'CTC': 'E', 'ACA': 'C', 'GCA': 'C', 'TCA': '_', 'CCA': 'W', 'ACG': 'R', 'GCG': 'R', 'TCG': 'R', 'CCG': 'R', 'ACT': 'S', 'GCT': 'S', 'TCT': 'R', 'CCT': 'R', 'ACC': 'G', 'GCC': 'G', 'TCC': 'G', 'CCC': 'G'}



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
        if line.startswith('CAGA'):
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
    final_reads = list()
    for k,v in counts.items():
        if v >= 700:
            for dna, aa in codontable.items():
                if k == dna:
                    final_reads.append([aa, v])
                else:
                    continue
    for k,v in final_reads:
        if v >= 700:
            print(k)
            print(v)

print('POS 11')
counter(65, 68)






print('DONE')


#this loop generates an histogram of the frequency of each barcode
