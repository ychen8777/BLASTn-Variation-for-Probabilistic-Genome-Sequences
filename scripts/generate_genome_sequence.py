
'''
Provided the following arguments, this script generates num genome sequences
and store in output_file:
  an FASTA file containing a genome sequence
  a file providing the probability of confidence for each nucleotide position in
  the FASTA file
  number of genome sequences to be generated
  length of each genome sequence
  path to output file containing the generated sequences
'''

import random
import argparse
from pathlib import Path







def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_file", help="path to FASTA file")
    parser.add_argument("prob_file", help="path to file for probability of confidence")
    parser.add_argument("-n", "--num_seq", help="number of sequences to be generated")
    parser.add_argument("-l", "--seq_length", help="length of each sequence to be generated")
    parser.add_argument("output", help="path to output file")

    args = parser.parse_args()

    fasta_file = Path(args.fasta_file)
    prob_file = Path(args.prob_file)

    if args.num_seq == None:
        num_seq = 100
    else:
        num_seq = int(args.num_seq)

    if args.seq_length == None:
        seq_length = 100
    else:
        seq_length = int(args.seq_length)

    output_file = Path(args.output)

    # print(fasta_file)
    # print(prob_file)
    # print("number", num_seq)
    # print("length", seq_length)
    # print(output_file)

if __name__ == '__main__':
    main()
