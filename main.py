import argparse
import subprocess
# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

##Add Tiger

def read_args():
	###############################
	# Set up argument parser
	argpar = argparse.ArgumentParser()
	# Positional (non-optional) arguments
	argpar.add_argument('-thread_n', type=int, nargs='+', required = True, help='the number of threads to use')
	argpar.add_argument('-trim_path', required = True, help='the pathway to Trimmoatic')
	argpar.add_argument('-in_fasta_f', required = True, help = 'The input of forward FASTA file')
	argpar.add_argument('-in_fasta_r', required = True, help = 'The input of reverse FASTA file')
	argpar.add_argument('-out_prefix', required = True, help='The pathway to The prefix that will be appended to the \
	output files. For ease of organization/alphabetization')
	argpar.add_argument('-out_loc', help = '') #Add default or required property
	argpar.add_argument('-adaptor_path', required = True, help = 'The FASTA file describing the adaptors for illumina\
	sequencing')
	argpar.add_argument('-ref_g', required= True, help='Reference genome for aligning')
	# Named (optional) arguments
	argpar.add_argument('-seed_mis', help = '', default=2)
	argpar.add_argument('-pal_clip', help = '', default=30)
	argpar.add_argument('-sim_clip', help = '', default=10)
	argpar.add_argument('-lead_q', help = '', default= 20)
	argpar.add_argument('-trail_q', help = '', default= 20)
	argpar.add_argument('-min_len', help = '', default= 30)

	###############################
	# Read command line args
	args = argpar.parse_args()

	###############################
	# Return command line args
	return args

def trimmomatics(pathToTrim=None, threadNum=None, in_fasta_f=None, in_fasta_r=None, out_prefix=None, out_loc=None,
				 adaptor_path=None, seed_mis=None, palClip=None, simClip=None,
				 leadQ=None, trailq=None, minlen=None):

	#Make the output files and put the names in a dictionary.
	out_pair1 = out_prefix + '_paired1.fq.gz'
	out_unpair1 = out_prefix + '_unpaired1.fq.gz'
	out_pair2 = out_prefix + '_paired2.fq.gz'
	out_unpair2 = out_prefix + '_unpaired2.fq.gz'
	out_files = [out_pair1, out_unpair1 , out_pair2, out_unpair2]

	#Run trimmomatic
	subprocess.run(['java', '-classpath', pathToTrim, 'org.usadellab.trimmomatic.TrimmomaticPE', '-threads', threadNum,
					'-phred33', in_fasta_f, in_fasta_r, out_pair1, out_unpair1, out_pair2, out_unpair2,
					'ILLUMINACLIP:', adaptor_path,':', seed_mis, ':', palClip,':', simClip,
					'LEADING:', leadQ,
					'TRAILING:', trailq,
					'MINLEN:', minlen])

	#Return the output files
	return out_pair1, out_unpair1 , out_pair2, out_unpair2

def run_bwa(thread_num, trimmed_for_inputs, trimmed_rev_inputs, ref_gen, out_prefix):
	subprocess.run('bwa', 'index', ref_gen, out_prefix)
	##Figure out trimmed_inputs
	subprocess.run('bwa', 'mem', '-a', '-C', '-H','-M',  '-P', '-t', thread_num, out_prefix, trimmed_for_inputs, trimmed_rev_inputs)
	##Figure out parameters
	#Used all parameters given, I suspect only -a, -M, -P is needed

def run_tiger():
	subprocess.run('octave', 'TIGER_generate_processing_files.m')

#def convert_to_


if __name__ == '__main__':
	print("readtoTiger has started")
	arguments = read_args()
	out_pair1, out_unpair1, out_pair2, out_unpair2 = trimmomatics(pathToTrim=arguments.trim_path,
																  threadNum=arguments.thread_n,
																  in_fasta_f=arguments.in_fasta_f,
																  in_fasta_r=arguments.in_fasta_r,
																  out_prefix=arguments.out_prefix,
																  out_loc=arguments.out_loc,
																  adaptor_path=arguments.adaptor_path,
																  seed_mis=arguments.seed_mis,
																  palClip=arguments.pal_clip,
																  simClip=arguments.sim_clip,
																  leadQ=arguments.lead_q, trailq=arguments.trail_q,
																  minlen=arguments.min_len)
	##Figure out which inputs are needed
	run_bwa(threadNum=arguments.threadNum, trimmed_for_inputs=out_pair1,
			trimmed_rev_inputs=out_pair2, ref_gen=arguments.ref_g, out_prefix=arguments.out_prefix)

