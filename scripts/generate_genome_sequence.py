
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

def generate_query_sequence(start, length, prob_db):
    ''' returns sequence generated based on probability in prob_db starting from
    start
    '''
    seq = ""
    pos = start

    for i in range(length):
        seq += choose_nuc(prob_db[pos+i])

    return seq
    #return (start, seq)

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

    # build sequence and probability db
    with open(fasta_file, "r") as f_file:
        sequence_db = f_file.readline()

    with open(prob_file, "r") as p_file:
        probability = p_file.readline().split()

    nucleotides = set(list(sequence_db))
    prob_db = create_prob(sequence_db, probability, nucleotides)

    #print(len(sequence_db), len(prob_db))

    # generate sequence and write to file
    db_length = len(sequence_db)
    with open(output_file, 'w') as o_file:
        for i in range(num_seq):
            start_pos = random.randint(0, db_length - seq_length - 1)
            seq = generate_query_sequence(start_pos, seq_length, prob_db)
            each_line = str(start_pos) + "," + seq + "\n"
            #print(each_line)
            o_file.write(each_line)

    print(f"Generated {num_seq} sequences of length {seq_length} in {output_file}")

if __name__ == '__main__':
    main()
