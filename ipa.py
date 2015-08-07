#!/usr/bin/python

"""
Script designed to pipeline a collection of NGS reads through various utilities in order to create a
reference assembly in FASTA format.

Author: Quinton Sirianni
"""

import argparse
import os
import subprocess
import sys


def check_env():
    """Check the execution environment to ensure all necessary utilities are installed"""

    required_utils = ['bamtofastq', 'fiona', 'bowtie2', 'samtools', 'bcftools', 'vcfutils.pl', 'fastq_to_fasta']

    print 'Ensuring that all necessary utilities are installed...'

    with open(os.devnull, 'w') as null:
        for util in required_utils:
            # POSIX environment
            if os.name == 'posix':
                if subprocess.call(['which', util], stdout=null) == 0:
                    print util, 'found'
                else:
                    raise OSError('{} was not found. Ensure it is installed and in your PATH'.format(util))

            # Windows environment
            elif os.name == 'nt':
                if subprocess.call(['where', '/Q','{}.exe'.format(util)], stdout=null) == 0:
                    print util, 'found'
                else:
                    raise OSError('{} was not found. Ensure it is installed and in your PATH'.format(util))

            # Unsupported environments
            else:
                raise OSError('This script cannot run in the current environment. Try using a POSIX or Windows NT '
                              'environment instead')

    print 'All necessary utilities have been found'


def get_args():
    """Return the arguments passed by the user to the script at runtime"""

    # Setup a parser object for user args
    parser = argparse.ArgumentParser(description='IonTorrent Pipeline Assembler', epilog='Use to understand why '
                                                                                         'S.cerevisae makes beer '
                                                                                         'taste so good!')

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
    """Process NGS reads into a reference assembly"""

    # Get user args
    args = get_args()

    #The user wishes to test the environment the script is executing in
    if args['env']:
        try:
            check_env()
        except OSError as e:
            print e.message
            sys.exit(1)


# Execute this script
if __name__ == '__main__':
    main()










