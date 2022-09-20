import argument_parser
import enviro_initialize
import call_trimmomatic
import call_bwa

if __name__ == '__main__':
	print("readtoTiger has started")

	arguments = argument_parser.read_args()

	output_dir_loc = enviro_initialize.create_output_directory(location=arguments.out_loc, prefix=arguments.out_prefix)

	out_pair_f, out_unpair_f, out_pair_r, out_unpair_r = call_trimmomatic.trimmomatics(path_to_trim=arguments.trim_path,
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

	call_bwa.run_bwa_index(ref_gen=arguments.ref_g)

	call_bwa.run_bwa_mem(thread_num=arguments.thread_n[0],
				trimmed_for_inputs=out_pair_f,
				trimmed_rev_inputs=out_pair_r,
				out_prefix=arguments.out_prefix,
				out_loc=output_dir_loc[2],
				ref_gen=arguments.ref_g)

	# TODO Add TIGER running
