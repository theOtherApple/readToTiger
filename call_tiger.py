import subprocess


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
