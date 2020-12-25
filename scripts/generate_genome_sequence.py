
'''
Provided the following arguments, this script generates num genome sequences
and store in output_file:
  an FASTA file containing a genome sequence
  a file providing the probability of confidence for each nucleotide position in
  the FASTA file
  number of genome sequences to be generated
  length of each genome sequence
  path to output file containing the generated sequences

If the number or length is not given, 100 sequences of length 100 will be genreted.

'''

import random
import argparse
from pathlib import Path

def create_prob(sequence_db, probability, nucleotides):
    ''' returns a list with probability of each nucleotide at each postion as
    a dictionary
    '''
    prob_db = []

    for i in range(len(sequence_db)):
        pred_nuc = sequence_db[i]
        confidence = float(probability[i])

        probs ={pred_nuc:confidence}
        other_confidence = (1-confidence) / (len(nucleotides) - 1)
        for nuc in nucleotides:
            if nuc not in probs:
                probs[nuc] = other_confidence

        prob_db.append(probs)

    return prob_db

def choose_nuc(nuc_prob):
    ''' returns a randomly chosen nucleotide based on the probability distribution
    '''
    nucs = []
    probs = []
    for nuc, prob in nuc_prob.items():
        nucs.append(nuc)
        probs.append(prob)

    #print(nuc_prob)
    #print(nucs)
    #print(probs)

    return random.choices(nucs, weights=probs, k=1)[0]



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

    with open(fasta_file) as f_file:
        sequence_db = f_file.readline()

    with open(prob_file) as p_file:
        probability = p_file.readline().split()

    nucleotides = set(list(sequence_db))
    prob_db = create_prob(sequence_db, probability, nucleotides)

    #print(len(sequence_db), len(prob_db))



if __name__ == '__main__':
    main()
