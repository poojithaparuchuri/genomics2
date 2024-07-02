"""
Sree Poojitha Paruchuri - Problem Set-2
Command to execute the code- python Problemset2.py -i spaceSeq.fa -o Problemset2.fasta -k 13
"""
# import libraries
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


motif = [
        [.5, .5, .5, .5, 0, 0, 0, 0, 0, 2, -99, -99, .5],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, -99, 2, -99, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, -99, -99, -99, 0],
        [.5, .5, .5, .5, 0, 0, 0, 0, 0, .5, -99, 2, 0]
    ]

def scoreMotif(seq, motif):
    base_idx = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    score = 0
    for i, base in enumerate(seq):
        if base in base_idx:
            score += motif[base_idx[base]][i]
    return score

def scanSeq(sequence, k, motif):
    orf_sequences = [] #lists to store the ORF sequences
    start_positions = [] #list to store the start positions of the sequences
    orf_lengths = [] #list to store the lengths of the identified ORF sequences
    motif_score = [] #list to store the motif scores of the ORF
    orf_number = 0 #counter for the number of identified ORFs

    i = 0
    while i < len(sequence) - k:  #iterate over the sequence
        seq = sequence[i:i +k]    #extracting subsequence of length k
        score = scoreMotif(seq, motif)   #score the subsequence
        if score > 7.25:  # checking if score exceeds threshold
            for j in range(i, len(sequence), 3):
                if sequence[j:j + 3] in ["ATG", "GTG"]: # chech for start codons
                    orf, start = " ", j
                    for x in range(j, len(sequence), 3):
                        codon = sequence[x:x + 3]
                        orf += codon  # adding codon to ORF
                        if codon in ["TAA", "TAG", "TGA"] and len(orf) >= 60: # checking for stop codon and minimum length
                            orf_sequences.append(orf) # store ORF
                            start_positions.append(i) # store start positions
                            orf_lengths.append(j + 3 - i) # store ORF lengths
                            motif_score.append(score) # store score
                            orf_number += 1 # increment ORF counter
                            break
                    break
        i += 1

    return orf_sequences, start_positions, orf_lengths, motif_score, orf_number


def identify_ORF(input_file, output_file, k, motif):
    with open(output_file, 'w') as out: #open output file
        for sequences in SeqIO.parse(input_file, "fasta"): #parse input fasta
            sequence = str(sequences.seq).upper()  #convert sequence to uppercase
            orf_sequences, start_positions, orf_lengths, motif_score, orf_number = scanSeq(sequence, k, motif)
            # Write ORFs to output file
            for i, (seq, start, length, score) in enumerate(zip(orf_sequences, start_positions, orf_lengths, motif_score)):
                orf_record = SeqRecord(
                    Seq(seq),
                    id=f"{sequences.id}_ORF{i + 1}_start{start}_len{length}_score{score:.2f}",
                    description=""
                )
                SeqIO.write(orf_record, out, "fasta")

if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(description="Identify ORFs in DNA sequences")
    arg_parser.add_argument("-i", "--input", required=True, help="Input FASTA file with contigs")
    arg_parser.add_argument("-o", "--output", required=True, help="Output FASTA file with identified ORFs")
    arg_parser.add_argument("-k", "--kmer", type=int, default=13, help="K-mer length for motif scoring (default: 13)")
    args = arg_parser.parse_args()

    identify_ORF(args.input, args.output, args.kmer, motif)

