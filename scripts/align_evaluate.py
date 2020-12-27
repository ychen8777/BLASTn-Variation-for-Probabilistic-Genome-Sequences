'''
This script performs local alignment, examines if the start position of the alignment
matches the record from generation, and saves results to file

'''
import argparse
from pathlib import Path

def read_seq(input_file):
    ''' reads start position and seqence from input_file and returns a list of
    tuples (start_postion, sequence)
    '''

    test_data = []

    with open(input_file, 'r') as in_file:
        each_line = in_file.readline()
        while len(each_line) > 0:
            pos_seq = each_line.split(",")
            test_data.append((pos_seq[0], pos_seq[1][:-1]))

            each_line = in_file.readline()

    return test_data

def create_kmer(seq, k):
    '''return a list of kmer from chars in seq
    '''
    start = 0
    kmer_list = []

    while start+k-1 <= len(seq) - 1:
        kmer_list.append(seq[start:start+k])
        start += 1

    return kmer_list

def compute_prob_score(query_nuc, match, mismatch, pos, prob_db):
    ''' return score of single position given query_nuc and its postion
    '''

    match_prob = prob_db[pos][query_nuc]
    score = match_prob * match + (1 - match_prob) * mismatch

    return score

def compute_fixed_length_score(seq1, seq2, match, mismatch, start=None, prob_db=None):
    ''' return alignment score for seq1 and seq2
    non-probabilistic score is returned if start and prob_db are not provided.
    '''

    score = 0

    if start == None or prob_db == None:
        # treat as non-probalistic
        for i in range(len(seq1)):
            if seq1[i] == seq2[i]:
                score += match
            else:
                score += mismatch

        return score

    cur_pos = start
    for nuc in seq1:
        score += compute_prob_score(nuc, match, mismatch, cur_pos, prob_db)
        cur_pos += 1

    #print(cur_pos)
    return score


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("seq_file", help="path to file containing start position and sequence")
    parser.add_argument("output_file", help="file to save results")
    parser.add_argument("-c", "--num_choices", help="number of top choices to consider")

    args = parser.parse_args()
    seq_file = Path(args.seq_file)
    output_file = Path(args.output_file)
    num_choices = int(args.num_choices)

    #print(seq_file)
    #print(output_file)
    #print(num_choices)

    test_data = read_seq(seq_file)
    #print(test_data)

if __name__ == '__main__':
    main()
