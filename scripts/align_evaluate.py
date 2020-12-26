'''
This script performs local alignment, examines if the start position of the alignment
matches the record from generation, and saves results to file

'''
import argparse
from pathlib import Path



def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("seq_file", help="path to file containing start position and sequence")
    parser.add_argument("output_file", help="file to save results")
    parser.add_argument("-c", "--num_choices", help="number of top choices to consider")

    args = parser.parse_args()
    seq_file = Path(args.seq_file)
    output_file = Path(args.output_file)
    num_choices = int(args.num_choices)

    print(seq_file)
    print(output_file)
    print(num_choices)


if __name__ == '__main__':
    main()
