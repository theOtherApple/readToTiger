import argparse
import subprocess
import os


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


def create_output_directory(location, prefix):
	# Directory all other directories will be put into with check for duplicate names
	start_dir_path = os.path.join(location, prefix)
	i = 1
	if os.path.isdir(start_dir_path):
		while os.path.isdir(start_dir_path + '_' + str(i)):
			i = i + 1
		start_dir_path = start_dir_path + '_' + str(i)
	os.mkdir(start_dir_path)
	print("Directory made for files at", start_dir_path)
	# Directory for trimmomatic outputs
	trim_dir_path = os.path.join(start_dir_path, 'trim')
	os.mkdir(trim_dir_path)
	# Directory for bwa outputs
	bwa_dir_path = os.path.join(start_dir_path, 'bwa')
	os.mkdir(bwa_dir_path)
	# Directory for tiger outputs
	tiger_dir_path = os.path.join(start_dir_path, 'tiger')
	os.mkdir(tiger_dir_path)
	# Library of all the paths to directories created
	dir_paths = [start_dir_path, trim_dir_path, bwa_dir_path, tiger_dir_path]
	return dir_paths


def trimmomatics(path_to_trim, thread_num, in_fasta_f, in_fasta_r, out_prefix, out_loc, adaptor_path, seed_mis,
				pal_clip, sim_clip, min_adapt_len, both_reads, lead_q, trail_q, min_len):

	# Make the output files and put the names in a dictionary.
	out_pair1 = out_prefix + '_forward_paired.fq.gz'
	out_unpair1 = out_prefix + '_forward_unpaired.fq.gz'
	out_pair2 = out_prefix + '_reverse_paired.fq.gz'
	out_unpair2 = out_prefix + '_reverse_unpaired.fq.gz'

	# TODO keep/remove os.path stuff
	out_pair1 = os.path.join(out_loc, out_pair1)
	out_unpair1 = os.path.join(out_loc, out_unpair1)
	out_pair2 = os.path.join(out_loc, out_pair2)
	out_unpair2 = os.path.join(out_loc, out_unpair2)

	#TODO Make this cleaner
	#TODO remove/keep os.path stuff
	trim_command = 'java -jar '
	trim_command += path_to_trim
	#trim_command += ' org.usadellab.trimmomatic.TrimmomaticPE'
	trim_command += ' PE -threads '
	trim_command += str(thread_num)
	trim_command += ' -phred33 '
	trim_command += in_fasta_f
	trim_command += ' '
	trim_command += in_fasta_r
	trim_command += ' '
	#trim_command += os.path.join(out_loc, out_pair1)
	trim_command += out_pair1
	trim_command += ' '
	#trim_command += os.path.join(out_loc, out_unpair1)
	trim_command += out_unpair1
	trim_command += ' '
	#trim_command += os.path.join(out_loc, out_pair2)
	trim_command += out_pair2
	trim_command += ' '
	#trim_command += os.path.join(out_loc, out_unpair2)
	trim_command += out_unpair2
	trim_command += ' ILLUMINACLIP:'
	trim_command += adaptor_path
	trim_command += ':'
	trim_command += str(seed_mis)
	trim_command += ':'
	trim_command += str(pal_clip)
	trim_command += ':'
	trim_command += str(sim_clip)
	trim_command += ':'
	trim_command += str(min_adapt_len)
	trim_command += ':'
	trim_command += str(both_reads)
	trim_command += ' LEADING:'
	trim_command += str(lead_q)
	trim_command += ' TRAILING:'
	trim_command += str(trail_q)
	trim_command += ' MINLEN:'
	trim_command += str(min_len)

	print(trim_command)

	# Run trimmomatic
	subprocess.Popen(trim_command, shell=True).wait()
	print('output got')

	# Return the output files
	return out_pair1, out_unpair1 , out_pair2, out_unpair2


def run_bwa_index(ref_gen):
	gen_index_input = 'bwa index ' + ref_gen
	print(gen_index_input)
	subprocess.Popen(gen_index_input, shell=True).wait()


def run_bwa_mem(thread_num, trimmed_for_inputs, trimmed_rev_inputs, out_prefix, out_loc, ref_gen):
	bwa_command = 'bwa mem -a -C -H -M -P -t '
	bwa_command += str(thread_num)
	bwa_command += ' '
	# bwa_command += ' -k '
	# bwa_command += min_seed_len
	# bwa_command += '-w'
	# bwa_command += band_width
	# bwa_command += '-d'
	# bwa_command += z_dropoff
	# bwa_command += '-r'
	# bwa_command += seed_split_ratio
	# bwa_command += '-c'
	# bwa_command += maxOcc
	# bwa_command += '-A'
	# bwa_command += matchScore
	# bwa_command += '-B'
	# bwa_command += mmPenalty
	# bwa_command += '-O'
	# bwa_command += gapOpenPen
	# bwa_command += '-E'
	# bwa_command += gapExtPen
	# bwa_command += '-L'
	# bwa_command += clipPen
	# bwa_command += '-U'
	# bwa_command += unpairPen
	# bwa_command += '-R'
	# bwa_command += RGline
	# bwa_command += '-v'
	# bwa_command += verboseLevel
	# TODO add in db.prefix (genome)
	bwa_command += ref_gen + ' '
	bwa_command += trimmed_for_inputs + ' '
	bwa_command += trimmed_rev_inputs + ' > '
	bwa_command += out_loc + '/' + out_prefix + '_aln-se.sam'
	print(bwa_command)
	subprocess.Popen(bwa_command, shell=True).wait()


def run_tiger():
	#TODO Finish this
	tiger_command = 'generate_chromosome_mappability_mask.sh --ref '
	tiger_command += REF_GENOME_PATH
	tiger_command += '--k '
	tiger_command += KMER_LENGTH
	tiger_command += '--chr'
	tiger_command += CHR
	print(tiger_command)
	subprocess.Popen(tiger_command, shell=True).wait()


if __name__ == '__main__':
	print("readtoTiger has started")

	arguments = read_args()

	output_dir_loc = create_output_directory(location=arguments.out_loc, prefix=arguments.out_prefix)

	out_pair_f, out_unpair_f, out_pair_r, out_unpair_r = trimmomatics(path_to_trim=arguments.trim_path,
																		thread_num=arguments.thread_n[0],
																		in_fasta_f=arguments.in_fasta_f,
																		in_fasta_r=arguments.in_fasta_r,
																		out_prefix=arguments.out_prefix,
																		out_loc=output_dir_loc[1],
																		adaptor_path=arguments.adaptor_path,
																		seed_mis=arguments.seed_mis,
																		pal_clip=arguments.pal_clip,
																		sim_clip=arguments.sim_clip,
																		min_adapt_len=arguments.min_adapt_len,
																		both_reads=arguments.both_reads,
																		lead_q=arguments.lead_q,
																		trail_q=arguments.trail_q,
																		min_len=arguments.min_len)

	run_bwa_index(ref_gen=arguments.ref_g)

	run_bwa_mem(thread_num=arguments.thread_n[0],
				trimmed_for_inputs=out_pair_f,
				trimmed_rev_inputs=out_pair_r,
				out_prefix=arguments.out_prefix,
				out_loc=output_dir_loc[2],
				ref_gen=arguments.ref_g)

	# TODO Add TIGER running
