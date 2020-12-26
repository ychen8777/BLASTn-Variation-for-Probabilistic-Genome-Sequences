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
