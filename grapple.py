#!/usr/bin/env python

"""
Script designed to pipeline a collection of NGS reads through various utilities in order to create a
reference assembly in FASTA format.

Author: Quinton Sirianni
"""

from __future__ import print_function

import multiprocessing
import os.path
import random
import re
import subprocess
import sys
import tempfile
from subprocess import CalledProcessError

import psutil


def error(message):
    """
    Print an error message in stderr

    message - the error message to print
    """

    # Ensure the stderr is being fed into a terminal
    if sys.stderr.isatty():
        # Print the error message as red text using escape codes
        print('\x1B[31m', message, '\x1B[0m', sep='', file=sys.stderr)

    else:
        # Print the message with no formatting
        print(message, file=sys.stderr)

    # Exit the script with an error code
    sys.exit(1)


def status(message):
    """
    Print a status message in stderr

    message - the status message to print
    """

    # Ensure the stderr is being fed into a terminal
    if sys.stderr.isatty():
        # Print the status message as green text using escape codes
        print('\x1B[32m', message, '\x1B[0m', sep='', file=sys.stderr)

    else:
        # Print the message with no formatting
        print(message, file=sys.stderr)


def bam_to_fq(read_file, prefix_id='', verbose=False):
    """
    Convert the input file from BAM to FASTQ using samtools.

    read_file - file containing the NGS reads in BAM format
    prefix_id - prefix of all temp files
    verbose - verbosity of subprocess

    Returns the FASTQ file
    """

    # Ensure that the file passed is in the proper format
    if os.path.splitext(read_file)[1] != '.bam':
        raise ValueError('The read file is not in BAM format')

    # Create a temporary output file to place the FASTQ output in
    ofile = os.path.join(tempfile.gettempdir(), prefix_id + 'bam_to_fq_out.fq')

    status('Converting the input from BAM format to FASTQ format')

    with open(ofile, 'w') as ofile_handle, open(os.devnull, 'w') as null_handle:
        err_handle = sys.stderr if verbose else null_handle

        # Convert the input from BAM to FASTQ
        subprocess.check_call(['samtools', 'bam2fq', read_file], stdout=ofile_handle, stderr=err_handle)

    return ofile


def read_correction(read_file, cell_type='haploid', match_type='edit', verbose=False):
    """
    Correct the raw reads using Karect.

    read_file - file containing the NGS reads in FASTQ format
    cell_type - type of cell the reads are from
    match_type - correction mode to be used
    verbose - verbosity of subprocess

    Returns the FASTQ corrected file
    """

    # Ensure the file is in FASTQ format
    if not re.match(r'\.((fastq)|(fq))', os.path.splitext(read_file)[1]):
        raise ValueError('The read file is not in FASTQ format')

    # Ensure the file exists (karect doesn't return an error code if it doesn't)
    if not os.path.isfile(read_file):
        raise ValueError("The read file doesn't exist")

    # Ensure that karect's parameters are valid (karect doesn't return an error code if they are not)
    if not re.match(r'(haploid)|(diploid)', cell_type):
        raise ValueError('The cell type is not a valid value')

    if not re.match(r'(edit)|(hamming)|(insdel)', match_type):
        raise ValueError('The match type is not a valid value')

    # Get system parameters
    thread_number = multiprocessing.cpu_count()
    memory_limit = psutil.virtual_memory().available / 1000000000

    status('Correcting the reads')

    with open(os.devnull, 'w') as null_handle:
        err_handle = sys.stderr if verbose else null_handle

        # Correct the reads
        # Note: Karect uses stdout rather than stderr for user information so stdout is redirected to err_handle
        subprocess.check_call(['karect', '-correct', '-inputfile=' + read_file, '-celltype=' + cell_type,
                               '-matchtype=' + match_type, '-threads=' + str(thread_number),
                               '-memory=' + str(memory_limit), '-resultdir=' + tempfile.gettempdir(),
                               '-tempdir=' + tempfile.gettempdir()], stdout=err_handle, stderr=err_handle)

    # Return the location of the output file
    ifile_suffix = os.path.split(read_file)[1]
    ofile = os.path.join(tempfile.gettempdir(), 'karect_' + ifile_suffix)

    return ofile


def read_alignment(read_file, ref_genome_file, prefix_id='', verbose=False):
    """
    Align the reads to the reference genome using Bowtie2.

    read_file - file containing the NGS reads to align in FASTQ format
    ref_genome_file - file containing the reference genome in FASTA format
    prefix_id - prefix of all temp files
    verbose - verbosity of subprocess

    Returns the aligned FASTQ read file
    """

    # Ensure the passed files are in the appropriate formats
    if not re.match(r'\.((fq)|(fastq))', os.path.splitext(read_file)[1]):
        raise ValueError('The read file is not in FASTQ format')

    if not re.match(r'\.((fa)|(fna)|(fasta))', os.path.splitext(ref_genome_file)[1]):
        raise ValueError('The reference genome file is not in FASTA format')

    # Get system parameters
    thread_number = multiprocessing.cpu_count()

    index_prefix = os.path.join(tempfile.gettempdir(), prefix_id + 'bt2_index')
    ofile = os.path.join(tempfile.gettempdir(), prefix_id + 'aligned_reads.sam')

    status('Aligning the reads')

    with open(os.devnull, 'w') as null_handle:
        err_handle = sys.stderr if verbose else null_handle

        # Create an index file from the reference genome
        subprocess.check_call(['bowtie2-build', ref_genome_file, index_prefix], stdout=null_handle, stderr=err_handle)

        with open(ofile, 'w') as ofile_handle:
            # Align the reads
            subprocess.check_call(['bowtie2', '-p', str(thread_number), '-x', index_prefix, '-U', read_file],
                                  stdout=ofile_handle, stderr=err_handle)

    return ofile


def sam_to_bam(read_file, prefix_id='', verbose=False):
    """
    Convert the input file from SAM format to BAM format

    read_file - reads in SAM format
    prefix_id - prefix of all temp files
    verbose - verbosity of subprocess

    Returns the converted BAM read file
    """

    # Ensure the read file is in SAM format
    if os.path.splitext(read_file)[1] != '.sam':
        raise ValueError('The read file is not in SAM format')

    ofile = os.path.join(tempfile.gettempdir(), prefix_id + 'aligned_reads.bam')

    status('Converting the aligned reads from SAM format to BAM format')

    with open(os.devnull, 'w') as null_handle:
        err_handle = sys.stderr if verbose else null_handle

        # Convert the read file format
        subprocess.check_call(['samtools', 'view', '-bo', ofile, read_file], stdout=null_handle, stderr=err_handle)

    return ofile


def sort_and_index(read_file, prefix_id='', verbose=False):
    """
    Sort and index the aligned reads

    read_file - aligned reads in SAM format
    prefix_id - prefix of all temp files
    verbose - verbosity of subprocess

    Returns the sorted and indexed SAM read file
    """

    # Ensure read file is in BAM format
    if os.path.splitext(read_file)[1] != '.bam':
        raise ValueError('The read file is not in BAM format')

    # Get system parameters
    thread_number = multiprocessing.cpu_count()

    temp_prefix = os.path.join(tempfile.gettempdir(), prefix_id + 'samtools_sorting')
    ofile = os.path.join(tempfile.gettempdir(), prefix_id + 'sorted_reads.bam')

    status('Sorting and indexing the reads')

    with open(os.devnull, 'w') as null_handle:
        err_handle = sys.stderr if verbose else null_handle

        # Sort the read file
        subprocess.check_call(['samtools', 'sort', '-o', ofile, '-@', str(thread_number), '-T', temp_prefix, read_file],
                              stdout=null_handle, stderr=err_handle)

        # Index the sorted reads
        subprocess.check_call(['samtools', 'index', ofile], stdout=null_handle, stderr=err_handle)

    return ofile


def call_variants(read_file, ref_genome_file, prefix_id='', verbose=False):
    """
    Call the variants in the read file using the reference genome

    read_file - sorted and indexed reads in BAM format
    ref_genome_file - reference genome in FASTA format
    prefix_id - prefix of all temp files
    verbose - verbosity of subprocess

    Returns the call variants file in VCF format
    """

    # Ensure the files are in the appropriate format
    if os.path.splitext(read_file)[1] != '.bam':
        raise ValueError('The read file is not in BAM format')

    if not re.match(r'\.((fa)|(fna)|(fasta))', os.path.splitext(ref_genome_file)[1]):
        raise ValueError('The reference genome file is not in FASTA format')

    pileup = os.path.join(tempfile.gettempdir(), prefix_id + 'pileup.vcf')
    variants = os.path.join(tempfile.gettempdir(), prefix_id + 'variants.vcf')
    ofile = os.path.join(tempfile.gettempdir(), prefix_id + 'consensus.fa')

    status('Calling the variants')

    with open(os.devnull, 'w') as null_handle:
        err_handle = sys.stderr if verbose else null_handle

        # Run mpileup
        subprocess.check_call(['samtools', 'mpileup', '-uf', ref_genome_file, '-o', pileup, read_file],
                              stdout=null_handle, stderr=err_handle)

        # Call the variants
        subprocess.check_call(['bcftools', 'call', '-mv', '-Oz', '-o', variants, pileup], stdout=err_handle,
                              stderr=null_handle)

        # Index the variants
        subprocess.check_call(['bcftools', 'index', variants], stdout=null_handle, stderr=err_handle)

        # Generate a consensus
        subprocess.check_call(['bcftools', 'consensus', '-f', ref_genome_file, '-o', ofile, variants],
                              stdout=null_handle, stderr=err_handle)

    return ofile


def format_consensus(consensus_file, prefix_id=''):
    """
    Edit the consenses file to ensure it is formatted correctly

    consensus_file - Consensus file to format
    prefix_id - prefix of all temp files
    """

    # Ensure the files are in the appropriate format
    if not re.match(r'\.((fa)|(fna)|(fasta))', os.path.splitext(consensus_file)[1]):
        raise ValueError('The consensus file is not in FASTA format')

    formatted_file = os.path.join(tempfile.gettempdir(), prefix_id + 'formatted_consensus.fa')

    status('Formatting the consensus')

    # Format the consensus reads
    with open(consensus_file) as consensus_handle, open(formatted_file, mode='w') as formatted_handle:
        for line in consensus_handle:
            if re.match(r'^>', line):
                formatted_handle.write(line)

            else:
                formatted_handle.write(line.upper())

    return formatted_file


def main(args):
    """Executes the pipeline according to the user's arguments."""

    try:
        # Start the pipeline if the user provided a reference genome
        if args['ref']:
            # Ensure the reference file exists
            if not os.path.isfile(args['ref']):
                raise IOError()

            else:
                # Generate a random identifier to label the temp files with
                prefix_id = str(random.getrandbits(32)) + '_'

                # Determine if the user has provided an input file or wishes to use stdin
                if args['input']:
                    # Ensure the file exists
                    if os.path.isfile(args['input']):
                        ifile = args['input']
                    else:
                        raise IOError()

                else:
                    ifile = os.path.join(tempfile.gettempdir(), prefix_id + 'stdin_dump.bam')
                    with open(ifile, 'wb') as ifile_handle:
                        for line in sys.stdin:
                            ifile_handle.write(line)

                # Convert the input file containing the reads from BAM to FASTQ format
                raw_reads = bam_to_fq(ifile, prefix_id, args['verbose'])

                # Correct the reads if error correction has not been disabled
                if not args['disable_ec']:
                    # Transform the command line arguments into values Karect can use
                    ploidy = 'haploid' if args['ploidy'] == 'n' else 'diploid'

                    if args['mode'] == 'equal':
                        mode = 'edit'

                    elif args['mode'] == 'indel':
                        mode = 'insdel'

                    else:
                        mode = 'hamming'

                    # Run Karect
                    corrected_reads = read_correction(raw_reads, ploidy, mode, args['verbose'])

                else:
                    corrected_reads = raw_reads

                # Align the reads
                aligned_reads = read_alignment(corrected_reads, args['ref'], prefix_id, args['verbose'])

                # Convert the aligned reads to BAM format from SAM format
                converted_aligned_reads = sam_to_bam(aligned_reads, prefix_id, args['verbose'])

                # Sort and index the aligned reads
                sorted_reads = sort_and_index(converted_aligned_reads, prefix_id, args['verbose'])

                # Call the variants and generate a consensus
                consensus = call_variants(sorted_reads, args['ref'], prefix_id, args['verbose'])

                # Clean up the consensus formatting
                cleaned_consensus = format_consensus(consensus, prefix_id)

                # Determine if the user has provided an output file or wishes to use stdout
                with open(cleaned_consensus) as consensus_handle:
                    if args['output']:
                        with open(args['output'], 'w') as ofile_handle:
                            for line in consensus_handle:
                                ofile_handle.write(line)

                    else:
                        for line in consensus_handle:
                            sys.stdout.write(line)

                status('The reference genome has been successfully assembled!')

        else:
            # Raise an exception since no reference was provided
            raise ValueError('A reference genome was not provided so the pipeline cannot execute')

    except KeyboardInterrupt:
        # Exit the script cleanly if interrupted by user
        error('')

    except ValueError as e:
        # Print the error message before exiting the script
        error(e)

    except EnvironmentError:
        # Inform the user something is wrong with the execution environment
        error('An error has occurred. Please ensure the input and reference files exist and all of the '
              'required utilities are installed in your PATH')

    except CalledProcessError as e:
        # Print an error message depending on which process failed
        if e.cmd[0] == 'samtools':
            if e.cmd[1] == 'bam2fq':
                error('The reads could not be converted from BAM to FASTQ format')

            elif e.cmd[1] == 'view':
                error('The reads could not be converted from SAM format to BAM format')

            elif e.cmd[1] == 'sort':
                error('The read file could not be sorted')

            elif e.cmd[1] == 'index':
                error('The sorted reads could not be indexed')

            elif e.cmd[1] == 'mpileup':
                error('The reads could not be processed by mpileup before being called')

        elif e.cmd[0] == 'karect':
            error('The reads could not be corrected')

        elif e.cmd[0] == 'bowtie2-build':
            error('An index could not be constructed from the reference genome provided')

        elif e.cmd[0] == 'bowtie2':
            error('The reads could not be aligned to the reference genome')

        elif e.cmd[0] == 'bcftools':
            if e.cmd[1] == 'call':
                error('The variants could not be called')

            elif e.cmd[1] == 'index':
                error('The variant index could not be constructed')

            elif e.cmd[1] == 'consensus':
                error('A consensus could not be generated from the variants')


if __name__ == '__main__':
    import argparse

    # Setup a parser object for user args
    parser = argparse.ArgumentParser(prog='grapple', description='Genome Reference Assembly Pipeline', add_help=False)

    parser.add_argument('-d', '--disable_ec', action='store_true', help='Disable error correction')

    parser.add_argument('-h', '--help', action='help', help='Display this help screen')

    parser.add_argument('-i', '--input', help='Specify an input file of NGS reads in BAM format. If this flag is '
                                              'not present, stdin is used instead')

    parser.add_argument('-o', '--output', help='Specify an output file for the resulting genome in FASTA format. If '
                                               'this flag is not present, stdout is used instead')

    parser.add_argument('-r', '--ref', help='The reference genome used to align the read in FASTA format')

    parser.add_argument('-V', '--version', action='version', version='{} 0.2.1'.format(parser.prog.title()),
                        help='Show the current version of the software.')

    parser.add_argument('-v', '--verbose', action='store_true', help='Output more information about each subprocess '
                                                                     'being executed')

    parser.add_argument('--ploidy', choices=['n', '2n'], default='n', help='Specify the ploidy of the cells from which '
                                                                           'the reads came from. If error correction '
                                                                           'is disabled, this option is ignored. '
                                                                           'Default value = n')

    parser.add_argument('--mode', choices=['equal', 'indel', 'subs'], default='equal',
                        help='Specify the type of errors that should be favoured during read correction '
                             'The equal option weighs all types of errors equally. If error correction is disabled, '
                             'this option is ignored. Default value = equal')

    # Retrieve the arguments and pass them to the main function
    main(vars(parser.parse_args()))
