#!/usr/bin/python

"""
Script designed to pipeline a collection of NGS reads through various utilities in order to create a
reference assembly in FASTA format.

Author: Quinton Sirianni
Version: 0.1
"""

import argparse
import os
import subprocess
import sys


def check_env():
    """Check the execution environment to ensure all necessary utilities are installed"""

    required_utils = ['karect', 'bowtie2', 'samtools', 'bcftools']

    print 'Ensuring that all necessary utilities are installed...'

    with open(os.devnull, 'w') as null_out:
        for util in required_utils:
            # POSIX environment
            if os.name == 'posix':
                if subprocess.call(['which', util], stdout=null_out) == 0:
                    print util, 'found'
                else:
                    raise OSError('{} was not found. Ensure it is installed and in your PATH'.format(util))

            # Unsupported environments
            else:
                raise OSError('This script is designed to execute in a POSIX environment only')

    print 'All necessary utilities have been found. You are ready to assemble'


def get_args():
    """Return the arguments passed by the user to the script at runtime"""

    # Setup a parser object for user args
    parser = argparse.ArgumentParser(description='IonTorrent Pipeline Assembler', epilog='Use to understand why '
                                                                                         'yeast makes beer taste so '
                                                                                         'good!')

    parser.add_argument('-e', '--env', action='store_true', help='Checks the execution environment for the '
                                                                 'necessary utilites before executing the '
                                                                 'pipeline')

    parser.add_argument('-i', '--input', type=int, help='Specify an input file of NGS reads in BAM format. If '
                                                        'flag not present, stdin is used instead')

    parser.add_argument('-o', '--output', type=int, help='Specify an output file for the resulting genome in '
                                                         'FASTA format. If flag not present, stdout is used '
                                                         'instead')

    parser.add_argument('-t', '--threads', type=int, help='The maximum number of threads to use. If flag not '
                                                          'present, this number is determined dynamically '
                                                          'based off the number of available virtual cores')

    parser.add_argument('REF_GENOME', nargs='*', help='The reference genome used to align the read in FASTA format')

    # Retrieve the arguments and return them as a dictionary
    return vars(parser.parse_args())


def main():
    # Get user args
    args = get_args()

    # The user wishes to test the environment the script is executing in
    if args['env']:
        try:
            check_env()
        except OSError as e:
            print e.message
            sys.exit(1)

    # Convert the input BAM file into FASTQ if the user provided a REF_GENOME
    if args['REF_GENOME'] is not None:
        pass


if __name__ == '__main__':
    main()










