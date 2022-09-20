import argparse
import csv
import time


def header_make(f, oct_v, loc):
    f.write("Created by Octave " + oct_v + ", " + str(time.time()) + " " + loc + '\n')
# Created by Octave 5.2.0, Tue Mar 29 17:40:12 2022 MDT <jgbaldwinbrown@jgbaldwinbrown-ThinkPad-T480>


def struct_make(f, name, struct_type, parameters):
    f.write("# name: " + name + '\n')
    f.write("# type: " + struct_type + '\n')
    if struct_type == "scalar struct":
        # parameters [0] = int(n_dimensions), [1] = str(structure length)
        n_dimensions = parameters[0]
        struc_len = parameters[1]
        f.write("# ndims: " + str(n_dimensions) + '\n')
        i = 0
        dim_written = "1"
        while i < n_dimensions - 1:
            dim_written += " 1"
            i = i + 1
        f.write(dim_written + '\n')
        f.write("# length: " + str(struc_len))
    elif struct_type == "cell":
        # parameters [0] = int(n_rows), [1] = int(n_cols)
        # Note that all elements within cell structs must be named "<cell-element>"
        rows = parameters[0]
        cols = parameters[1]
        f.write("# rows: " + str(rows) + '\n')
        f.write("# columns: " + str(cols))
    # TODO: make this an exception
    else:
        print("struct not recognized")


def scalar_make(f, name, scalar_type, parameters):
    f.write("# name: " + name + '\n')
    f.write("# type: " + scalar_type + '\n')
    if scalar_type == "sq_string":
        # parameters [0] = int(n_elements), [1] = str(content)
        n_elements = parameters[0]
        content = parameters[1]
        f.write("# elements: " + str(n_elements))
        for i in content:
            f.write('\n' +"# length: " + str(len(i)) + '\n' + i)
    elif scalar_type == "matrix":
        # parameters [0] = int(n_rows), [1] = int(n_cols), [3]
        # TODO: make other matrx sizes possible
        rows = parameters[0]
        cols = parameters[1]
        content = parameters[2]
        chrom = content[0]
        start = content[1]
        end = content[2]
        f.write("# rows: " + str(rows) + '\n')
        f.write("# columns: " + str(cols) + '\n')
        i = 0
        for j in chrom:
            f.write(str(chrom[i]) + '\t' + str(start[i]) + '\t' + str(end[i]) + '\n')
            i = i + 1

    # TODO: make this an exception
    else:
        print("scalar not recognized")

# # name: fasta
# # type: scalar struct
# # ndims: 2
#  1 1
# # length: 2
# # name: Header
# # type: sq_string
# # elements: 2
# # length: 2
# 3R
# # length: 2
# 3L
# # name: Sequence
# # type: sq_string
# # elements: 2
# # length: 5000
# acgggaccgagtatagtacc


if __name__ == '__main__':
    # Set up argument parser
    argpar = argparse.ArgumentParser()
    # Positional (non-optional) arguments
    argpar.add_argument('-f', help='The fasta gap file')
    argpar.add_argument('-o', help='The output, matlab approved, file')

    arguments = argpar.parse_args()
    # TODO remove this testing line of code that overrides arguments
    arguments_file = arguments.f
    #arguments_file = '/Volumes/melements/tiger_pipeline-main/example/gap.txt'
    arguments_file = '/Volumes/melements/dm6.fa_copy'

    chr_a = []
    start_a = []
    end_a = []
    chr_dict = {
        "chr2L": 20,
        "chr2R": 21,
        "chr3L": 30,
        "chr3R": 31,
        "chr4": 40,
        "chrX": 10,
        "chrY": 11,
        "chrX_": 12,
        "chrY_": 13,
        "chrUn_": 0
    }

    with open(arguments_file) as fasta_gap:
        tsv_file = csv.reader(fasta_gap, delimiter="\t")

        for line in tsv_file:
            chr_name = line[1]
            if len(chr_name) <= 5:
                if "chrX" in chr_name:
                    chr_name = "chrX"
                elif "chrY" in chr_name:
                    chr_name = "chrY"
                elif "chr4" in chr_name:
                    chr_name = "chr4"
                elif "chr3L" in chr_name:
                    chr_name = "chr3L"
                elif "chr3R" in chr_name:
                    chr_name = "chr3R"
                elif "chr2L" in chr_name:
                    chr_name = "chr2L"
                elif "chr2R" in chr_name:
                    chr_name = "chr2R"
            else:
                if "chrX" in chr_name:
                    chr_name = "chrX_"
                elif "chrY" in chr_name:
                    chr_name = "chrY_"
                else:
                    chr_name = "chrUn_"
            chr_a.append(chr_dict.get(chr_name))
            start_a.append(line[2])
            end_a.append(line[3])

    #matfile_path = '/Volumes/melements/tiger_pipeline-main/example/test.m'
    matfile_path = '/Volumes/melements/tiger_pipeline-main/example/full_test.m'

    if ((len(chr_a) != len(start_a)) or (len(end_a) != len(start_a))):
        print("error")
        # TODO make this an actual exception
    try:
        # matfile = open(arguments.o, "x")
        # TODO remove this testing line of code and replace it with one that actually gets the proper location on file
        matfile = open(matfile_path, "x")
    except FileExistsError:
        # matfile = open(arguments.o, "w")
        # TODO remove this testing line of code and replace it with one that actually gets the proper location on file
        matfile = open(matfile_path, "w")

    # TODO Remove all but the necessary lines
    # header_make(matfile, 'oct_v', 'loc')
    # matfile.write('\n')
    # struct_make(matfile, 'Gap', 'cell', [1, 3])
    # matfile.write('\n')
    # scalar_make(matfile, '<cell-element>', 'sq_string', [len(chr_a), chr_a])
    # matfile.write('\n')
    # scalar_make(matfile, '<cell-element>', 'matrix', [1, len(start_a), start_a])
    # matfile.write('\n')

    #scalar_make(matfile, 'Gap', 'matrix', [len(end_a), 3, [chr_a, start_a, end_a]])
    scalar_make(matfile, 'Sequence', 'matrix', [len(end_a), 3, [chr_a, start_a, end_a]])
