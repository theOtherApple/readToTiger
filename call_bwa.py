import subprocess


def run_bwa_index(ref_gen):
    gen_index_input = 'bwa index ' + ref_gen
    print(gen_index_input)
    subprocess.Popen(gen_index_input, shell=True).wait()


def run_bwa_mem(thread_num, trimmed_for_inputs, trimmed_rev_inputs, out_prefix, out_loc, ref_gen):
    print("bwa has started")
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
