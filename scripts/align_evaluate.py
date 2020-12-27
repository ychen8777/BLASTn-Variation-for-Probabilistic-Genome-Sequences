'''
This script performs local alignment, examines if the start position of the alignment
matches the record from generation, and saves results to file

'''
import re
import argparse
from pathlib import Path
from collections import Counter

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

def search_hit(query, sequence_db, match, mismatch, threshold, prob_db=None):
    ''' returns a Counter with start position of hit and corresponding
    alignment score
    '''
    len_query = len(query)
    len_seq_db = len(sequence_db)
    start = 0

    result = Counter()

    if prob_db == None:
        # do not consider probability
        while start+len_query <= len_seq_db:
            sub_sequence = sequence_db[start:start+len_query]
            #print(sub_sequence)

            alignment_score = compute_fixed_length_score(query, sub_sequence, match, mismatch)
            if alignment_score >= threshold:
                #result.append((start, alignment_score))
                result[start] = alignment_score

            #print(query, sub_sequence, alignment_score)
            start += 1

    else:
        # probabilistic sequence
        while start+len_query <= len_seq_db:

            alignment_score = compute_fixed_length_score(query, "", match, mismatch, start, prob_db)

            if alignment_score >= threshold:
                #result.append((start, alignment_score))
                result[start] = alignment_score

            start += 1

    return result

def extension(query_seq, kmer, sequence_db, start, score, match, mismatch, cutoff, prob_db=None):
    ''' returns the start and end positions in sequence_db, and score after
    extension step
    '''

    seq_length = len(query_seq)
    db_length = len(sequence_db)
    pos_in_query = [m.start() for m in re.finditer(f"(?={kmer})", query_seq)]
    #print(pos_in_query)

    #start_results = Counter()
    #end_results = Counter()

    results = Counter()

    for pos in pos_in_query:

        # extend left

        high_score = score
        cur_score = score
        high_start = start
        db_pos = start    # pos in database sequence
        query_pos = pos   # pos in query_seq

        while db_pos > 0 and query_pos > 0 and cur_score - high_score >= cutoff:
            db_pos -= 1
            query_pos -= 1

            query_nuc = query_seq[query_pos]

            if prob_db == None:
                # non-probabilistic
                if query_nuc == sequence_db[db_pos]:
                    incre_score = match
                else:
                    incre_score = mismatch

            else:
                # probabilistic
                incre_score = compute_prob_score(query_nuc, match, mismatch, db_pos, prob_db)

            cur_score += incre_score

            if cur_score > high_score:
                high_score = cur_score
                high_start = db_pos

            #print("left", cur_score, incre_score, high_score, high_start, query_pos)

        # extend right

        db_pos = start + len(kmer) - 1
        query_pos = pos + len(kmer) - 1
        high_end = db_pos

        while db_pos < db_length - 1 and query_pos < seq_length - 1 and cur_score - high_score >= cutoff:
            db_pos += 1
            query_pos += 1

            query_nuc = query_seq[query_pos]

            if prob_db == None:
                # non-probabilistic
                if query_nuc == sequence_db[db_pos]:
                    incre_score = match
                else:
                    incre_score = mismatch

            else:
                # probabilistic
                incre_score = compute_prob_score(query_nuc, match, mismatch, db_pos, prob_db)

            cur_score += incre_score

            if cur_score > high_score:
                high_score = cur_score
                high_end = db_pos

            #print("right", cur_score, incre_score, high_score, high_end, query_pos)

        results[(high_start, high_end)] = high_score
        #print("result", high_start, high_end-high_start+1, high_score)

    #return start_results.most_common()
    #return (high_start, length, high_score)
    return results

def find_local_alignment(query_seq, word_size, sequence_db, match, mismatch, threshold, extension_cutoff, prob_db=None):
    ''' return the start and end positions of local alignments in sequence_db
    '''

    # create list of kmers
    kmer_list = create_kmer(query_seq, word_size)

    # search for hits
    kmer_index = {}
    for word in kmer_list:
        kmer_index[word] = search_hit(word, sequence_db, match, mismatch, threshold, prob_db)
    #print(kmer_index)

    # extend hits
    local_alignments = Counter()
    for word in kmer_index:
        for db_start, score in kmer_index[word].items():
            extensions = extension(query_seq, word, sequence_db, db_start, score, match, mismatch, extension_cutoff, prob_db)
            for start_end, alignment_score in extensions.items():
                if start_end in local_alignments:
                    if alignment_score > local_alignments[start_end]:
                        local_alignments[start_end] = alignment_score
                else:
                    local_alignments[start_end] = alignment_score

    return local_alignments

def count_correctness(true_start_pos, seq_length, local_alignments, choices):
    ''' return increment for total_correct, total_start_correct, and
    total_end_correct numbers
    '''

    #choices = len(local_alignments)

    total_correct = [0] * choices
    total_start_correct = [0] * choices
    total_end_correct = [0] * choices

    true_end_pos = true_start_pos + seq_length - 1
    for i in range(len(local_alignments)):
        pre_start_pos = local_alignments[i][0][0]
        pre_end_pos = local_alignments[i][0][1]

        if pre_start_pos == true_start_pos and pre_end_pos == true_end_pos:
            for j in range(i, choices):
                total_correct[j] = 1

        if pre_start_pos == true_start_pos and pre_end_pos != true_end_pos:
            for j in range(i, choices):
                total_start_correct[j] = 1

        if pre_start_pos != true_start_pos and pre_end_pos == true_end_pos:
            for j in range(i, choices):
                total_end_correct[j] = 1
    for i in range(len(total_correct)):
        if total_correct[i] == 1:
            total_start_correct[i] = 0
            total_end_correct[i] = 0

    return total_correct, total_start_correct, total_end_correct

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
