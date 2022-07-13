import argparse
import subprocess


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
	argpar.add_argument('-lead_q', help='', default=20)
	argpar.add_argument('-trail_q', help='', default=20)
	argpar.add_argument('-min_len', help='', default=30)

	# Read command line args
	args = argpar.parse_args()

	# Return command line args
	return args


def trimmomatics(pathToTrim, threadNum, in_fasta_f, in_fasta_r, out_prefix, out_loc, adaptor_path, seed_mis, pal_clip,
				 sim_clip, leadq, trailq, min_len):

	# Make the output files and put the names in a dictionary.
	out_pair1 = out_prefix + '_paired1.fq.gz'
	out_unpair1 = out_prefix + '_unpaired1.fq.gz'
	out_pair2 = out_prefix + '_paired2.fq.gz'
	out_unpair2 = out_prefix + '_unpaired2.fq.gz'

	#TODO Make this cleaner
	trimcommand = 'java -classpath '
	trimcommand += pathToTrim
	#trimcommand +=' org.usadellab.trimmomatic.TrimmomaticPE -threads '
	trimcommand += ' PE -threads '
	trimcommand += str(threadNum)
	trimcommand += ' -phred33 '
	trimcommand += in_fasta_f
	trimcommand += ' '
	trimcommand += in_fasta_r
	trimcommand += ' '
	trimcommand += out_pair1
	trimcommand += ' '
	trimcommand += out_unpair1
	trimcommand += ' '
	trimcommand += out_pair2
	trimcommand += ' '
	trimcommand += out_unpair2
	trimcommand += ' ILLUMINACLIP:'
	trimcommand += adaptor_path
	trimcommand += ':'
	trimcommand += str(seed_mis)
	trimcommand += ':'
	trimcommand += str(pal_clip)
	trimcommand += ':'
	trimcommand += str(sim_clip)
	trimcommand += ' LEADING:'
	trimcommand += str(leadq)
	trimcommand += ' TRAILING:'
	trimcommand += str(trailq)
	trimcommand += ' MINLEN:'
	trimcommand += str(min_len)

	print(trimcommand)

	# TODO Run java program with proper arguments
	# TODO figure out why org.usadellab.trimmomatic.TrimmomaticPE and PE is not found
	# Run trimmomatic
	# print(subprocess.run(trimcommand))
	print(subprocess.run(trimcommand))

	# Return the output files
	return out_pair1, out_unpair1 , out_pair2, out_unpair2


def run_bwa(thread_num, trimmed_for_inputs, trimmed_rev_inputs, ref_gen, out_prefix):
	#TODO Finish this
	subprocess.run('bwa', 'index', ref_gen, out_prefix)
	# TODO Figure out trimmed inputs
	subprocess.run('bwa', 'mem', '-a', '-C', '-H', '-M', '-P', '-t', thread_num, out_prefix, trimmed_for_inputs, trimmed_rev_inputs)
	# TODO Figure out parameters
	# NOTE Used all parameters given, I suspect only -a, -M, -P is needed


def run_tiger():
	#TODO Finish this
	subprocess.run('octave', 'TIGER_generate_processing_files.m')


if __name__ == '__main__':
	print("readtoTiger has started")
	arguments = read_args()
	# out_pair1, out_unpair1, out_pair2, out_unpair2 = \
	print(trimmomatics(pathToTrim=arguments.trim_path, threadNum=arguments.thread_n[0], in_fasta_f=arguments.in_fasta_f,
					   in_fasta_r=arguments.in_fasta_r,
					   out_prefix=arguments.out_prefix,
					   out_loc=arguments.out_loc,
					   adaptor_path=arguments.adaptor_path,
					   seed_mis=arguments.seed_mis,
					   pal_clip=arguments.pal_clip,
					   sim_clip=arguments.sim_clip,
					   leadq=arguments.lead_q,
					   trailq=arguments.trail_q,
					   min_len=arguments.min_len))
	# TODO Add BWA running
	# run_bwa(threadNum=arguments.threadNum, trimmed_for_inputs=out_pair1,
	# 		trimmed_rev_inputs=out_pair2, ref_gen=arguments.ref_g, out_prefix=arguments.out_prefix)
	# TODO Add TIGER running
