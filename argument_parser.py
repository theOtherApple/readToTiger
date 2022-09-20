import argparse


def read_args():
	# Set up argument parser
	argpar = argparse.ArgumentParser()
	# Positional (non-optional) arguments
	argpar.add_argument('-thread_n', type=int, nargs='+', required=True, help='the number of threads to use')
	argpar.add_argument('-trim_path', required=True, help='the pathway to Trimmoatic')
	argpar.add_argument('-in_fasta_f', required=True, help='The input of forward FASTA file')
	argpar.add_argument('-in_fasta_r', required=True, help='The input of reverse FASTA file')
	argpar.add_argument('-out_prefix', required=True, help='The pathway to The prefix that will be appended to the \
	output files. For ease of organization/alphabetization')
	# TODO Add default or required property for -out_loc
	argpar.add_argument('-out_loc', help='')
	argpar.add_argument('-adaptor_path', required=True, help='The FASTA file describing the adaptors for illumina\
	sequencing')
	argpar.add_argument('-ref_g', required=True, help='Reference genome for aligning')
	# Named (optional) arguments
	argpar.add_argument('-seed_mis', help='', default=2)
	argpar.add_argument('-pal_clip', help='', default=30)
	argpar.add_argument('-sim_clip', help='', default=10)
	argpar.add_argument('-min_adapt_len', help='', default=2)
	argpar.add_argument('-both_reads', help='', default=True)
	argpar.add_argument('-lead_q', help='', default=20)
	argpar.add_argument('-trail_q', help='', default=20)
	argpar.add_argument('-min_len', help='', default=30)

	# Read command line args
	args = argpar.parse_args()

	# Return command line args
	return args
