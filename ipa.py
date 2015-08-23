#!/usr/bin/env python

"""
Script designed to pipeline a collection of NGS reads through various utilities in order to create a
reference assembly in FASTA format.

Author: Quinton Sirianni
"""

from __future__ import print_function

import os
import os.path
import subprocess
import sys
import tempfile
from subprocess import CalledProcessError

import psutil


def check_env():
    """Check the execution environment to ensure all necessary utilities are installed."""

    required_utils = ['karect', 'bowtie2', 'bowtie2-build', 'samtools', 'bcftools']

    print('Ensuring that all necessary utilities are installed...', file=sys.stderr)

    with open(os.devnull, 'w') as null_handle:
        for util in required_utils:
            # POSIX environment
            if os.name == 'posix':
                return_code = subprocess.call(['which', util], stdout=null_handle)

                # Ensure the subprocess executed correctly
                if return_code == 0:
                    print(util, 'found', file=sys.stderr)
                else:
                    raise CalledProcessError(util + 'was not found. Ensure it is installed and in your PATH')

            # Unsupported environments
            else:
                raise CalledProcessError('This script is designed to execute in a POSIX environment only')

    print('All necessary utilities have been found. You are ready to assemble', file=sys.stderr)


def bam_to_fq(read_file):
    """
    Convert the input file from BAM to FASTQ using SAMtools.

    bam_file = file containing the NGS reads in BAM format
    """

    # Create a temporary output file to place the FASTQ output in
    ofile = os.path.join(tempfile.gettempdir(), 'bam_to_fq_out_IPA.fq')

    print('Converting the input from BAM format to FASTQ format... ', end='', file=sys.stderr)

    with open(ofile, 'w') as ofile_handle, open(os.devnull, 'w') as null_handle:
        return_code = subprocess.call(['samtools', 'bam2fq', read_file], stdout=ofile_handle, stderr=null_handle)

        # Ensure the subprocess executed correctly
        if return_code != 0:
            raise CalledProcessError('The input could not be converted to FASTQ format')

    print('done!', file=sys.stderr)

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

    print('Correcting the reads... ', end='', file=sys.stderr)

    with open(os.devnull, 'w') as null_handle:
        return_code = subprocess.call(['karect', '-correct', '-inputfile=' + read_file, '-celltype=' + cell_type,
                                       '-matchtype='+ match_type, '-threads=' + thread_number, '-memory=' + memory_limit,
                                       'resultdir=' + tempfile.gettempdir(), '-tempdir=' + tempfile.gettempdir()],
                                        stdout=null_handle, stderr=null_handle)

        # Ensure the subprocess executed correctly
        if return_code != 0:
            raise CalledProcessError('The reads could not be corrected')

    print('done!', file=sys.stderr)

    # Return the location of the output file
    ifile_suffix = os.path.split(read_file)[1]
    ofile = os.path.join(tempfile.gettempdir(), 'karect_' + ifile_suffix)
    return ofile


def read_alignment(read_file, ref_genome_file, thread_number):
    """
    Align the reads to thre reference genome using Bowtie2.

    read_file - file containing the NGS reads to align in FASTQ format
    ref_genome_file - file containing the reference genome in FASTA format
    thread_number - number of threads to use during alignment
    """

    index_prefix = os.path.join(tempfile.gettempdir(), 'bt2_index_IPA')
    ofile = os.path.join(tempfile.gettempdir(), 'aligned_reads_IPA.sam')

    print('Aligning the reads... ', end='', file=sys.stderr)

    with open(os.devnull, 'w') as null_handle:
        # Create an index file from the reference genome
        return_code = subprocess.call(['bowtie2-build', ref_genome_file, index_prefix], stdout=null_handle,
                                        stderr=null_handle)

        # Ensure the subprocess executed successfully
        if return_code != 0:
            raise CalledProcessError('An index could not be constructed from the reference genome provided')

        with open(ofile, 'w') as ofile_handle:
            # Align the reads
            return_code = subprocess.call(['bowtie2', '-p', str(thread_number), '-x', index_prefix, '-U', read_file],
                                            stdout=ofile_handle, stderr=null_handle)

            # Ensure the subprocess executed successfully
            if return_code != 0:
                raise CalledProcessError('The reads could not be aligned to the reference genome')

    print('done!', file=sys.stderr)

    return ofile


def main(args):
    """Executes the pipeline according to the user's arguments."""

    # The user wishes to test the environment the script is executing in
    if args['env']:
        try:
            check_env()
        except CalledProcessError as e:
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

        except CalledProcessError as e:
            print(e.message, file=sys.stderr)
            sys.exit(1)

    else:
        print('A reference genome was not provided so the pipeline cannot execute', file=sys.stderr)
        sys.exit(1)


# Acquire arguments from the user and pass them to main function if this script is executed
if __name__ == '__main__':
    import argparse

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

    # Retrieve the arguments and pass them to the main function
    main(vars(parser.parse_args()))