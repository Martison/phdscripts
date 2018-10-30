import os
from sys import argv


#open file

print('Opening...')
filename = input('Enter file name with fq extension to analyze >: ')
handle = open(os.path.expanduser(f"~/Desktop/{filename}"))
print('File opened')

codontable = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M','ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T','AAC':'N', 'AAT':'N', 'AAA':'K','AAG':'K','AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R','CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',}

reverse_ct = {'AAA': 'F', 'GAA': 'F', 'TAA': 'L', 'CAA': 'L', 'AAG': 'L', 'GAG': 'L', 'TAG': 'L', 'CAG': 'L', 'AAT': 'I', 'GAT': 'I', 'TAT': 'I', 'CAT': 'M', 'AAC': 'V', 'GAC': 'V', 'TAC': 'V', 'AGA': 'S', 'GGA': 'S', 'TGA': 'S', 'CGA': 'S', 'AGG': 'P', 'GGG': 'P', 'TGG': 'P', 'CGG': 'P', 'AGT': 'T', 'GGT': 'T', 'TGT': 'T', 'CGT': 'T', 'AGC': 'A', 'GGC': 'A', 'TGC': 'A', 'CGC': 'A', 'ATA': 'Y', 'GTA': 'Y', 'TTA': '_', 'CTA': '_', 'ATG': 'H', 'GTG': 'H', 'TTG': 'Q', 'CTG': 'Q', 'ATT': 'N', 'GTT': 'N', 'TTT': 'K', 'CTT': 'K', 'ATC': 'D', 'GTC': 'D', 'TTC': 'E', 'CTC': 'E', 'ACA': 'C', 'GCA': 'C', 'TCA': '_', 'CCA': 'W', 'ACG': 'R', 'GCG': 'R', 'TCG': 'R', 'CCG': 'R', 'ACT': 'S', 'GCT': 'S', 'TCT': 'R', 'CCT': 'R', 'ACC': 'G', 'GCC': 'G', 'TCC': 'G', 'CCC': 'G'}


#change number for analyzing up to this # of reads

readstocount = int(input('Enter a max number of reads to count >: '))

#this loop trims the 72 first nucleotides of the sequence


def counter(pos1,pos1e, pos2, pos2e):
    trimmed = list()
    final_reads = list()
    trimming_count = 0
    counts = dict()
    handle.seek(0)
    for line in handle:
        line.rstrip()
        if line.startswith('TCTA'):
            trimmed.append(line[pos1:pos1e] + '_' + line[pos2:pos2e])
            trimming_count += 1
            if trimming_count >= readstocount:
                break
            else:
                continue
    for barcode in trimmed:
            counts[barcode] = counts.get(barcode, 0) + 1
    for k,v in counts.items():
        if v >= 1:
            print(k)
            print(v)
    for k,v in counts.items():
            for dna, aa in codontable.items():
                if k[0:3] == dna:
                    final_reads.append([aa + k[3:], v])
                else:
                    continue
    for k,v in final_reads:
            print(k)
            print(v)
    print('Position:', pos1,'-',pos1e)

    print(final_reads)

    print('**********************')
    print('**********************')
    print('**********************')
    print('**********************')

    #NEXT

    for k, v in final_reads:
        for dna, aa in reverse_ct.items():
            if k[2:] == dna:
                final_reads.append([k[0] + '_' + aa, v])
            else:
                continue

    print(final_reads)

    for k, v in final_reads:
        if v >= 1500:
            print(k,',',v)


counter(87, 90, 215, 218)









print('DONE')


#this loop generates an histogram of the frequency of each barcode
