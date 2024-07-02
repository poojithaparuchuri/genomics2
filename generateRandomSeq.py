import random
import matplotlib.pyplot as plt

l = 150 #read length
k = 10000
bases = ["A","T","C","G"]
def generateReads(l,k): #prints k reads of length l

    sequences = []
    for i in range(k):
        seq = ""
        for j in range(l):
            seq += bases[random.randint(0,3)]
        sequences.append(seq)
        print(">sequence " + str(i+1))
        print(seq)
    return sequences

def find_ORFs(sequence):
    start_codon = 'ATG'
    stop_codon = ['TAA', 'TAG','TGA']
    orfs = []
    i = 0
    while i < len(sequence) - 2:
        if sequence[i:i + 3] == start_codon:
            for j in range(i+3, len(sequence)-2, 3):
                if sequence[j:j +3] in stop_codon:
                    orfs.append(sequence[i:j+3])
                    break
            i += 3
        else:
            i += 1
    return orfs

def get_orf_lengths(sequences):
    orf_lengths = []
    for seq in sequences:
        orfs = find_ORFs(seq)
        orf_lengths.extend([len(orf) for orf in orfs])
    return orf_lengths


sequences = generateReads(l,k)
orf_lengths = get_orf_lengths(sequences)

plt.hist(orf_lengths, bins=range(min(orf_lengths), max(orf_lengths)+1,3), color= 'yellow' ,edgecolor='black')
plt.title('Distribution of ORF length')
plt.xlabel('ORF length(nucleotides)')
plt.ylabel('Frequency')
plt.show()



