import subprocess


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
    return out_pair1, out_unpair1, out_pair2, out_unpair2
