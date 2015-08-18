#!/usr/bin/env python

"""
Script designed to pipeline a collection of NGS reads through various utilities in order to create a
reference assembly in FASTA format.

Author: Quinton Sirianni
"""

from __future__ import print_function

import argparse
import os
import subprocess
import sys
import tempfile

import psutil


def check_env():
    """Check the execution environment to ensure all necessary utilities are installed."""

    required_utils = ['karect', 'bowtie2', 'bowtie2-build', 'samtools', 'bcftools']

    print('Ensuring that all necessary utilities are installed...')

    with open(os.devnull, 'w') as null_handle:
        for util in required_utils:
            # POSIX environment
            if os.name == 'posix':
                if subprocess.call(['which', util], stdout=null_handle) == 0:
                    print(util, 'found')
                else:
                    raise OSError(util + 'was not found. Ensure it is installed and in your PATH')
            # Unsupported environments
            else:
                raise OSError('This script is designed to execute in a POSIX environment only')

    print('All necessary utilities have been found. You are ready to assemble')


def get_args():
    """Return the arguments passed by the user to the script at runtime."""

    # Setup a parser object for user args
    parser = argparse.ArgumentParser(description='IonTorrent Pipeline Assembler')

    parser.add_argument('-d', '--disable_ec', action='store_true', help='Disable error correction')

    parser.add_argument('-e', '--env', action='store_true', help='Check the execution environment for the '
                                                                 'necessary utilites. If this option is selected, all '
                                                                 'other options are ignored and the main pipeline will '
                                                                 'not be executed')

    parser.add_argument('-i', '--input', help='Specify an input file of NGS reads in BAM format. If flag not present, '
                                              'stdin is used instead')

    parser.add_argument('-o', '--output', help='Specify an output file for the resulting genome in FASTA format. If '
                                               'flag not present, stdout is used instead')

    parser.add_argument('-r', '--ref', help='The reference genome used to align the read in FASTA format.')

    # Retrieve the arguments and return them as a dictionary
    return vars(parser.parse_args())


def bam_to_fq(read_file):
    """
    Convert the input file from BAM to FASTQ using SAMtools.

    bam_file = file containing the NGS reads in BAM format
    """

    print('Converting input from BAM format into FASTQ...')

    # Create a temporary output file to place the FASTQ output in
    ofile = os.path.join(tempfile.gettempdir(), 'bam_to_fq_out_IPA.fq')

    with open(ofile, 'w') as ofile_handle:
        subprocess.call(['samtools', 'bam2fq', read_file], stdout=ofile_handle)

    print('Read format converted!')

    return ofile


def read_correction(read_file, thread_number, memory_limit, cell_type='haploid', match_type='edit'):
    """
    Correct the raw reads using Karect.

    read_file = file containing the NGS reads in FASTQ format
    thread_number = the number of threads the subprocess can use
    memory_limit = the maximum amount of memory the subprocess can use
    cell_type = the type of cell the reads are from
    match_type = the correction mode to be used
    """

    print('Correcting the reads...')

    with open(os.devnull, 'w') as null_handle:
        subprocess.call(['karect', '-correct', '-inputfile=' + read_file, '-celltype=' + cell_type,
                         '-matchtype='+ match_type, '-threads=' + thread_number, '-memory=' + memory_limit,
                         'resultdir=' + tempfile.gettempdir(), '-tempdir=' + tempfile.gettempdir()], stdout=null_handle)

    print('Reads corrected!')

    # Return the location of the output file
    ifile_suffix = os.path.split(read_file)[1]
    ofile = os.path.join(tempfile.gettempdir(), 'karect_' + ifile_suffix)
    return ofile


def read_alignment(read_file, ref_genome_file, thread_number):
    """
    Align the reads to thre reference genome using Bowtie2.

    read_file - file containing the NGS reads to align in FASTQ format
    ref_genome_file - file containing the reference genome in FASTA format
    """

    print('Creating an index of the reference genome...')

    index_prefix = os.path.join(tempfile.gettempdir(), 'bt2_index_IPA')
    ofile = os.path.join(tempfile.gettempdir(), 'aligned_reads_IPA.sam')

    with open(os.devnull, 'w') as null_handle:
        # Create an index file from the reference genome
        subprocess.call(['bowtie2-build', ref_genome_file, index_prefix], stdout=null_handle)

    print('Reference genome indexed!')

    print('Aligning the reads against the reference...')

    with open(ofile, 'w') as ofile_handle:
        # Align the reads
        subprocess.call(['bowtie2', '-p', str(thread_number), '-x', index_prefix, '-U', read_file], stdout=ofile_handle)

    print('Reads aligned!')

    return ofile


def main():
    """Executes the pipeline according to the user's arguments."""

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
    elif args['ref']:
        # Get system information
        available_threads = psutil.cpu_count()
        available_memory = psutil.virtual_memory().total / 1000000000 / 2

        # Determine if the user has provided an input file or wishes to use stdin
        if args['input']:
            ifile = args['input']
        else:
            ifile = os.path.join(tempfile.gettempdir(), 'stdin_dump_IPA.bam')
            with open(ifile, 'wb') as ifile_handle:
                for line in sys.stdin:
                    ifile_handle.write(line)
        try:
            # Convert the input file containing the reads from BAM to FASTQ format
            raw_reads = bam_to_fq(ifile)

            # Correct the reads if error correction has not been disabled
            if not args['disable_ec']:
                corrected_reads = read_correction(raw_reads, available_threads, available_memory)
            else:
                corrected_reads = raw_reads

            # Align the reads
            aligned_reads = read_alignment(corrected_reads, args['ref'], available_threads)

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
        print('A reference genome was not provided so the pipeline cannot execute')


if __name__ == '__main__':
    main()










