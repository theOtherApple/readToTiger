import argparse
import csv


if __name__ == '__main__':
    # Set up argument parser
    argpar = argparse.ArgumentParser()
    # Positional (non-optional) arguments
    argpar.add_argument('-f', help='The fasta gap file')
    argpar.add_argument('-o', help='The output, matlab approved, file')

    arguments = argpar.parse_args()
    #TODO remove this testing line of code that overrides arguments
    arguments.f = '/Volumes/melements/tiger_pipeline-main/example/gap.txt'

    chr_a = []
    start_a = []
    end_a = []

    with open(arguments.f) as fasta_gap:
        tsv_file = csv.reader(fasta_gap, delimiter="\t")

        for line in tsv_file:
            chr_a.append(line[1])
            start_a.append(line[2])
            end_a.append(line[3])
            print('new line')
            print('0: ', line[0])
            print('1: ', line[1])
            print('2: ', line[2])
            print('3: ', line[3])
            print('4: ', line[4])

    #TODO write a matlab accepaed file
