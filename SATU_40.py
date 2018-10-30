import os
from sys import argv

#open file

print('Opening...')
filename = input('Enter file name with fq extension to analyze >: ')
handle = open(os.path.expanduser(f"~/Desktop/{filename}"))
print('File opened')

codontable = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M','ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T','AAC':'N', 'AAT':'N', 'AAA':'K','AAG':'K','AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R','CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',}



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
        if line.startswith('ACCT'):
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

counter(68, 71)







print('DONE')


#this loop generates an histogram of the frequency of each barcode
