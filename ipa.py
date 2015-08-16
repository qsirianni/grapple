#!/usr/bin/python

"""
Script designed to pipeline a collection of NGS reads through various utilities in order to create a
reference assembly in FASTA format.

Author: Quinton Sirianni
"""

from __future__ import print_function

import argparse
import multiprocessing
import os
import subprocess
import sys
import tempfile


def check_env():
    """Check the execution environment to ensure all necessary utilities are installed"""

    required_utils = ['karect', 'bowtie2', 'samtools', 'bcftools']

    print('Ensuring that all necessary utilities are installed...')

    with open(os.devnull, 'w') as null_handle:
        for util in required_utils:
            # POSIX environment
            if os.name == 'posix':
                if subprocess.call(['which', util], stdout=null_handle) == 0:
                    print(util, 'found')
                else:
                    raise OSError('{} was not found. Ensure it is installed and in your PATH'.format(util))

            # Unsupported environments
            else:
                raise OSError('This script is designed to execute in a POSIX environment only')

    print('All necessary utilities have been found. You are ready to assemble')


def get_args():
    """Return the arguments passed by the user to the script at runtime"""

    # Setup a parser object for user args
    parser = argparse.ArgumentParser(description='IonTorrent Pipeline Assembler')

    parser.add_argument('-e', '--env', action='store_true', help='Checks the execution environment for the '
                                                                 'necessary utilites before executing the '
                                                                 'pipeline')

    parser.add_argument('-i', '--input', help='Specify an input file of NGS reads in BAM format. If flag not present, '
                                              'stdin is used instead')

    parser.add_argument('-o', '--output', help='Specify an output file for the resulting genome in FASTA format. If '
                                               'flag not present, stdout is used instead')

    parser.add_argument('REF_GENOME', nargs='?', help='The reference genome used to align the read in '
                                                      'FASTA format')

    # Retrieve the arguments and return them as a dictionary
    return vars(parser.parse_args())


def bam_to_fq(bam_file):
    """Convert the input file from BAM to FASTQ using samtools"""

    print('Converting input from BAM format into FASTQ...')

    # Create a temporary output file to place the FASTQ output into
    ofile = '{}/bam_to_fq_out_IPA.fq'.format(tempfile.gettempdir())
    with open(ofile, 'w') as ofile_handle:
        # Check if a BAM file location was provided for input
        ifile = None
        if bam_file:
            # Use BAM file for the input
            ifile = bam_file
        else:
            # Dump stdin into a temporary input file
            ifile = '{}/stdin_dump_IPA.bam'.format(tempfile.gettempdir())
            with open(ifile, 'wb') as ifile_handle:
                for line in sys.stdin:
                    ifile_handle.write(line)

        # Call the subprocess with the input and output files
        subprocess.call(['samtools', 'bam2fq', ifile], stdout=ofile_handle)

    return ofile




def main():
    """Executes the pipeline according to the user's arguments"""

    # Get user args
    args = get_args()

    # The user wishes to test the environment the script is executing in
    if args['env']:
        try:
            check_env()
        except OSError as e:
            print(e.message, file=sys.stderr)
            sys.exit(1)

    # Start the pipeline if the user provided a reference genome
    if args['REF_GENOME']:
        try:
            # Convert the input file containing the reads from BAM to FASTQ format
            raw_reads = bam_to_fq(args['input'])

            # Correct the reads
            #corr_reads = read_correction(raw_reads)

            # Align the reads
            #aligned_reads = read_alignment(corr_reads, args['REF_GENOME'][0])

            # Sort and index the aligned reads
            #sorted_reads = sort_and_index(aligned_reads)

            # Call the variants
            #variants = call_variants(sorted_reads, args['REF_GENOME'][0])

            # Create a consensus
            #consensus = create_consensus(variants, args['REF_GENOME'][0])

            # Print the consensus to the output file or stdout
            #print_consensus(consensus)

        # Temporary "catch all" except clause. Will be fine-tuned as development continues
        except Exception as e:
            print(e.message, file=sys.stderr)
            sys.exit(1)
    else:
        print('A reference genome was not provided so the pipeline cannot run')


if __name__ == '__main__':
    main()










