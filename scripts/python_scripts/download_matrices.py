import sys
import subprocess
from subprocess import Popen, PIPE, call, check_output
import os
from os import listdir
from os.path import isfile, join
import fnmatch
import numpy as np


def execute_command(command):
    print(">>>>> Executing command {}".format(command))
    process = subprocess.Popen(command, stdout=subprocess.PIPE,
            stderr = subprocess.STDOUT, shell=True,
            universal_newlines=True)
    return_code = 0 #process.wait()  FIXME: find a way to check output
    output = process.stdout.read()
    return output, return_code

def print_name(mat_list):
    for mat in mat_list:
        mat_split = mat.split('/')
        print(mat_split[1].strip(), end='')
        print(',', end='')
    print("")

def download_suitesparse_mats(mat_list, out_dir):
    err_files = []
    print('Downloading matrices ...')
    if not os.path.exists(out_dir):
        print('Creating a directory in {}'.format(out_dir))
        os.makedirs(out_dir)
    for mat in mat_list:
        cur_mat = 'https://sparse.tamu.edu/MM/' + mat + '.tar.gz'
        mat_splitted = mat.split('/')
        cur_out = out_dir+mat_splitted[len(mat_splitted)-1] + '.tar.gz'
        print('Downloading {} into {} ... '.format(cur_mat, cur_out))
        output, ret_code = execute_command("wget {} -O {}".format(cur_mat, cur_out))
        print(output, 'return code: ', ret_code)
        if ret_code != 0:
            err_files = np.concatenate((err_files, cur_mat))
    return err_files


def extract_suitesparse_mats(directory):
    '''
    Extract all tar.gz file in 'directory' and copy the matrix.mtx file there.
    :param directory:
    :return: False if directory does not exist
    '''
    if not os.path.exists(directory):
        print('The path is wrong!')
        return False
    only_files = [f for f in listdir(directory) if isfile(join(directory, f))]
    only_tars = fnmatch.filter(only_files, '*.tar.gz')
    for f in only_tars:
        cur_tar = join(directory, f)
        execute_command('tar xvf {} -C {}'.format(cur_tar, directory))
        tar_splitted = f.split('.')
        cur_mat = join(directory, tar_splitted[0], tar_splitted[0]+'.mtx')
        execute_command('cp {} {}'.format(cur_mat, directory))
    return True

def clean_up(directory):
    only_dirs = [f for f in listdir(directory) if not isfile(join(directory, f))]
    only_files = [f for f in listdir(directory) if isfile(join(directory, f))]
    only_tars = fnmatch.filter(only_files, '*.tar.gz')
    for d in only_dirs:
        execute_command('rm -f -r {}'.format(join(directory,d)))
    for f in only_tars:
        execute_command('rm -f {}'.format(join(directory, f)))


def main(argv):
    if len(argv) < 2:
        print('Missing matrix type argument')
        print('download_matrices.py -matrix_type   -output_directory')
        return
    mat_type = argv[0]
    output_dir = argv[1]
    mat_list = []
    if mat_type == '1':  # easy indefinit matrices
        mat_list = ['Oberwolfach/t2dal',
                    'GHS_indef/dixmaanl',
                    'Oberwolfach/rail_79841',
                    'GHS_indef/dawson5',
                    'Boeing/bcsstk39',
                    'GHS_indef/helm2d03',
                    'GHS_indef/copter2',
                    'Boeing/crystk03',
                    'Oberwolfach/filter3D',
                    'Boeing/pct20stif',
                    'Koutsovasilis/F2',
                    'Cunningham/qa8fk',
                    'Oberwolfach/gas_sensor',
                    'McRae/ecology1',
                    'Oberwolfach/t3dh',
                    'Lin/Lin',
                    'PARSEC/H2O',
                    'GHS_indef/sparsine',
                    'PARSEC/Ge99H100',
                    'PARSEC/Ga10As10H30',
                    'PARSEC/Ga19As19H42']
    elif mat_type == '2':  # difficult indefinite matrices
        mat_list = ['TSOPF/TSOPF_FS_b39_c7',
                    'TSOPF/TSOPF_FS_b162_c1',
                    'TSOPF/TSOPF_FS_b39_c19',
                    'QY/case39',
                    'TSOPF/TSOPF_FS_b39_c30',
                    'GHS_indef/cont-201',
                    'GHS_indef/stokes128',
                    'TSOPF/TSOPF_FS_b162_c3',
                    'TSOPF/TSOPF_FS_b162_c4',
                    'GHS_indef/ncvxqp1',
                    'GHS_indef/darcy003',
                    'GHS_indef/cont-300',
                    'GHS_indef/bratu3d',
                    'GHS_indef/cvxqp3',
                    'TSOPF/TSOPF_FS_b300',
                    'TSOPF/TSOPF_FS_b300_c1',
                    'GHS_indef/d_pretok',
                    'GHS_indef/turon_m',
                    'TSOPF/TSOPF_FS_b300_c2',
                    'TSOPF/TSOPF_FS_b300_c3',
                    'GHS_indef/ncvxqp5',
                    'GHS_indef/ncvxqp3',
                    'GHS_indef/ncvxqp7',
                    'Schenk_IBMNA/c-big' ]
    elif mat_type == '3':  # kkt matrices
        mat_list = ['GHS_indef/ncvxqp1',
                    'GHS_indef/ncvxqp5',
                    'GHS_indef/ncvxqp3',
                    'GHS_indef/ncvxqp7',
                    'GHS_indef/ncvxqp9',
                    'GHS_indef/boyd1',
                    'GHS_indef/boyd2',
                    'GHS_indef/cont-201',
                    'GHS_indef/cont-300',
                    'VDOL/dynamicSoaringProblem_1',
                    'VDOL/dynamicSoaringProblem_2',
                    'VDOL/dynamicSoaringProblem_3',
                    'VDOL/dynamicSoaringProblem_4',
                    'VDOL/dynamicSoaringProblem_5',
                    'VDOL/dynamicSoaringProblem_6',
                    'VDOL/dynamicSoaringProblem_7',
                    'VDOL/dynamicSoaringProblem_8',
                    'VDOL/freeFlyingRobot_1',
                    'VDOL/freeFlyingRobot_2',
                    'VDOL/freeFlyingRobot_3',
                    'VDOL/freeFlyingRobot_4',
                    'VDOL/freeFlyingRobot_5',
                    'VDOL/freeFlyingRobot_6',
                    'VDOL/freeFlyingRobot_7',
                    'VDOL/freeFlyingRobot_8',
                    'VDOL/freeFlyingRobot_9',
                    'VDOL/freeFlyingRobot_10',
                    'VDOL/freeFlyingRobot_11',
                    'VDOL/freeFlyingRobot_12',
                    'VDOL/freeFlyingRobot_13',
                    'VDOL/freeFlyingRobot_14',
                    'VDOL/freeFlyingRobot_15',
                    'VDOL/freeFlyingRobot_16',
                    'VDOL/goddardRocketProblem_1',
                    'VDOL/goddardRocketProblem_2',
                    'VDOL/hangGlider_1',
                    'VDOL/hangGlider_2',
                    'VDOL/hangGlider_3',
                    'VDOL/hangGlider_4',
                    'VDOL/hangGlider_5',
                    'Zaoui/kkt_power',
                    'VDOL/lowThrust_1',
                    'VDOL/lowThrust_2',
                    'VDOL/lowThrust_3',
                    'VDOL/lowThrust_4',
                    'VDOL/lowThrust_5',
                    'VDOL/lowThrust_6',
                    'VDOL/lowThrust_7',
                    'VDOL/lowThrust_8',
                    'VDOL/lowThrust_9',
                    'VDOL/lowThrust_10',
                    'VDOL/lowThrust_11',
                    'VDOL/lowThrust_12',
                    'VDOL/lowThrust_13',
                    'VDOL/orbitRaising_1',
                    'VDOL/orbitRaising_2',
                    'VDOL/orbitRaising_3',
                    'VDOL/orbitRaising_4',
                    'VDOL/reorientation_1',
                    'VDOL/reorientation_2',
                    'VDOL/reorientation_3',
                    'VDOL/reorientation_4',
                    'VDOL/reorientation_5',
                    'VDOL/reorientation_6',
                    'VDOL/reorientation_7',
                    'VDOL/reorientation_8',
                    'VDOL/spaceShuttleEntry_1',
                    'VDOL/spaceShuttleEntry_2',
                    'VDOL/spaceShuttleEntry_3',
                    'VDOL/spaceShuttleEntry_4',
                    'VDOL/spaceStation_1',
                    'VDOL/spaceStation_2',
                    'VDOL/spaceStation_3',
                    'VDOL/spaceStation_4',
                    'VDOL/spaceStation_5',
                    'VDOL/spaceStation_6',
                    'VDOL/spaceStation_7',
                    'VDOL/spaceStation_8',
                    'VDOL/spaceStation_9',
                    'VDOL/spaceStation_10',
                    'VDOL/spaceStation_11',
                    'VDOL/spaceStation_12',
                    'VDOL/spaceStation_13',
                    'VDOL/spaceStation_14',
                    'VDOL/tumorAntiAngiogenesis_1',
                    'VDOL/tumorAntiAngiogenesis_2',
                    'VDOL/tumorAntiAngiogenesis_3',
                    'VDOL/tumorAntiAngiogenesis_4',
                    'VDOL/tumorAntiAngiogenesis_5',
                    'VDOL/tumorAntiAngiogenesis_6',
                    'VDOL/tumorAntiAngiogenesis_7',
                    'VDOL/tumorAntiAngiogenesis_8'
                    ]

    print_name(mat_list)
    # err_files = download_suitesparse_mats(mat_list, output_dir)
    # if len(err_files) > 0:
    #     print('These files are not downloaded:', err_files)
    # extract_suitesparse_mats(output_dir)
    # clean_up(output_dir)


if __name__ == "__main__":
    main(sys.argv[1:])
