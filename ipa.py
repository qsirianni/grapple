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
                    print('"{}"'.format(util), 'found', file=sys.stderr)
                else:
                    raise OSError('"{}"'.format(util) + 'was not found. Ensure it is installed and in your PATH')

            # Unsupported environments
            else:
                raise OSError('This script is designed to execute in a POSIX environment only')

    print('All necessary utilities have been found. You are ready to assemble', file=sys.stderr)


def bam_to_fq(read_file):
    """
    Convert the input file from BAM to FASTQ using samtools.

    bam_file - file containing the NGS reads in BAM format

    Returns the FASTQ file
    """

    # Create a temporary output file to place the FASTQ output in
    ofile = os.path.join(tempfile.gettempdir(), 'bam_to_fq_out.fq')

    print('Converting the input from BAM format to FASTQ format... ', end='', file=sys.stderr)

    with open(ofile, 'w') as ofile_handle, open(os.devnull, 'w') as null_handle:
        return_code = subprocess.call(['samtools', 'bam2fq', read_file], stdout=ofile_handle, stderr=null_handle)

        # Ensure the subprocess executed correctly
        if return_code != 0:
            raise OSError('The input could not be converted to FASTQ format')

    print('done!', file=sys.stderr)

    return ofile


def read_correction(read_file, thread_number, memory_limit, cell_type='haploid', match_type='edit'):
    """
    Correct the raw reads using Karect.

    read_file - file containing the NGS reads in FASTQ format
    thread_number - number of threads the subprocess can use
    memory_limit - maximum amount of memory the subprocess can use
    cell_type - type of cell the reads are from
    match_type - correction mode to be used

    Returns the FASTQ corrected file
    """

    print(read_file, file=sys.stderr)

    print('Correcting the reads... ', end='', file=sys.stderr)

    with open(os.devnull, 'w') as null_handle:
        return_code = subprocess.call(['karect', '-correct', '-inputfile=' + read_file, '-celltype=' + cell_type,
                                       '-matchtype=' + match_type, '-threads=' + str(thread_number),
                                       '-memory=' + str(memory_limit), '-resultdir=' + tempfile.gettempdir(),
                                       '-tempdir=' + tempfile.gettempdir()], stdout=null_handle, stderr=null_handle)

        # Ensure the subprocess executed correctly
        if return_code != 0:
            raise OSError('The reads could not be corrected')

    print('done!', file=sys.stderr)

    # Return the location of the output file
    ifile_suffix = os.path.split(read_file)[1]
    ofile = os.path.join(tempfile.gettempdir(), 'karect_' + ifile_suffix)

    print(ofile, file=sys.stderr)

    return ofile


def read_alignment(read_file, ref_genome_file, thread_number):
    """
    Align the reads to the reference genome using Bowtie2.

    read_file - file containing the NGS reads to align in FASTQ format
    ref_genome_file - file containing the reference genome in FASTA format
    thread_number - number of threads to use during alignment

    Returns the aligned FASTQ read file
    """

    index_prefix = os.path.join(tempfile.gettempdir(), 'bt2_index')
    ofile = os.path.join(tempfile.gettempdir(), 'aligned_reads.sam')

    print('Aligning the reads... ', end='', file=sys.stderr)

    with open(os.devnull, 'w') as null_handle:
        # Create an index file from the reference genome
        return_code = subprocess.call(['bowtie2-build', ref_genome_file, index_prefix], stdout=null_handle,
                                      stderr=null_handle)

        # Ensure the subprocess executed successfully
        if return_code != 0:
            raise OSError('An index could not be constructed from the reference genome provided')

        with open(ofile, 'w') as ofile_handle:
            # Align the reads
            return_code = subprocess.call(['bowtie2', '-p', str(thread_number), '-x', index_prefix, '-U', read_file],
                                          stdout=ofile_handle, stderr=null_handle)

            # Ensure the subprocess executed successfully
            if return_code != 0:
                raise OSError('The reads could not be aligned to the reference genome')

    print('done!', file=sys.stderr)

    return ofile


def sam_to_bam(read_file):
    """
    Convert the input file from SAM format to BAM format

    read_file - reads in SAM format

    Returns the converted BAM read file
    """

    ofile = os.path.join(tempfile.gettempdir(), 'aligned_reads.bam')

    print('Converting the aligned reads from SAM format to BAM format... ', end='', file=sys.stderr)

    with open(os.devnull, 'w') as null_handle:
        # Convert the read file format
        return_code = subprocess.call(['samtools', 'view', '-bo', ofile, read_file], stdout=null_handle,
                                      stderr=null_handle)

        # Ensure the subprocess executed successfully
        if return_code != 0:
            raise OSError('The reads could not be converted from SAM format to BAM format')

    print('done!', file=sys.stderr)

    return ofile


def sort_and_index(read_file, thread_number):
    """
    Sort and index the aligned reads

    read_file - aligned reads in SAM format
    thread_number - number of threads to use for sorting

    Returns the sorted and indexed SAM read file
    """

    temp_prefix = os.path.join(tempfile.gettempdir(), 'samtools_sorting')
    ofile = os.path.join(tempfile.gettempdir(), 'sorted_reads.bam')

    print('Sorting and indexing the reads... ', end='', file=sys.stderr)

    with open(os.devnull, 'w') as null_handle:
        # Sort the read file
        return_code = subprocess.call(['samtools', 'sort', '-o', ofile, '-@', str(thread_number), '-T', temp_prefix,
                                      read_file], stdout=null_handle, stderr=null_handle)

        # Ensure the subprocess executed successfully
        if return_code != 0:
            raise OSError('The read file could not be sorted')

        # Index the sorted file
        return_code = subprocess.call(['samtools', 'index', ofile], stdout=null_handle, stderr=null_handle)

        # Ensure the subprocess executed successfully
        if return_code != 0:
            raise OSError('The sorted reads could not be indexed')

    print('done!', file=sys.stderr)

    return ofile


def call_variants(read_file, ref_genome_file):
    """
    Call the variants in the read file using the reference genome

    read_file - sorted and indexed reads in BAM format
    ref_genome_file - reference genome in FASTA format

    Returns the call variants file in VCF format
    """

    pileup = os.path.join(tempfile.gettempdir(), 'pileup.vcf')
    variants = os.path.join(tempfile.gettempdir(), 'variants.vcf')
    ofile = os.path.join(tempfile.gettempdir(), 'consensus.fa')

    print('Calling the variants... ', end='', file=sys.stderr)

    with open(os.devnull, 'w') as null_handle:
        # Run mpileup
        return_code = subprocess.call(['samtools', 'mpileup', '-uf', ref_genome_file, '-o', pileup, read_file],
                                      stdout=null_handle, stderr=null_handle)

        # Ensure the process executed successfully
        if return_code != 0:
            raise OSError('The reads could not be processed before being called')

        # Call the variants
        return_code = subprocess.call(['bcftools', 'call', '-mv', '-Oz', '-o', variants, pileup], stdout=null_handle,
                                      stderr=null_handle)

        # Ensure the process executed successfully
        if return_code != 0:
            raise OSError('The variants could not be called')

        # Index the variants
        return_code = subprocess.call(['bcftools', 'index', variants], stdout=null_handle, stderr=null_handle)

        # Ensure the process executed successfully
        if return_code != 0:
            raise OSError('The variant index could not be constructed')

        # Generate a consensus
        return_code = subprocess.call(['bcftools', 'consensus', '-f', ref_genome_file, '-o', ofile, variants],
                                      stdout=null_handle, stderr=null_handle)

        # Ensure the process executed successfully
        if return_code != 0:
            raise OSError('A consensus could not be generated from the variants')

    print('done!', file=sys.stderr)

    return ofile


def main(args):
    """Executes the pipeline according to the user's arguments."""

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
            ifile = os.path.join(tempfile.gettempdir(), 'stdin_dump.bam')
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

            # Convert the aligned reads to BAM format from SAM format
            converted_aligned_reads = sam_to_bam(aligned_reads)

            # Sort and index the aligned reads
            sorted_reads = sort_and_index(converted_aligned_reads, available_threads)

            # Call the variants and generate a consensus
            consensus = call_variants(sorted_reads, args['ref'])

            # Determine if the user has provided an output file or wishes to use stdout
            with open(consensus) as consensus_handle:
                if args['output']:
                    with open(args['output'], 'w') as ofile_handle:
                        for line in consensus:
                            ofile_handle.write(line)

                else:
                    for line in consensus_handle:
                        sys.stdout.write(line)

            print('The reference genome has been successfully assembled!', file=sys.stderr)

        except OSError as e:
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
