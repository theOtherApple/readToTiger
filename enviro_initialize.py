import os


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
