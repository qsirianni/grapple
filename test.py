#!/usr/bin/env python

"""Contains unit tests for Grapple."""

import os.path
import unittest
from subprocess import CalledProcessError
from unittest import TestCase

import grapple


class TestError(TestCase):
    """Tests involving error()"""

    def test_exit(self):
        """Should raise an exception to exit the script when called with any message"""

        with self.assertRaises(SystemExit):
            grapple.error('')


class TestBamToFq(TestCase):
    """Tests involving bam_to_fq()"""

    def setUp(self):
        """Setup code for test cases"""

        # Available test file
        self._test_file = os.path.join('test_files', 'lambda_iontorrent.bam')

    def test_valid_file(self):
        """Should not raise an exception when a valid bam file is used"""

        try:
            grapple.bam_to_fq(self._test_file)

        except Exception as e:
            self.fail(e)

    def test_invalid_file(self):
        """Should raise an exception when the wrong type of file is used"""

        with self.assertRaises(ValueError):
            grapple.bam_to_fq(os.path.join('test_files', 'lambda_ref.fa'))

    def test_absent_file(self):
        """Should raise an exception when the file does not exist"""

        with self.assertRaises(CalledProcessError):
            grapple.bam_to_fq('this_file_does_not_exist.bam')

    def test_no_file(self):
        """Should raise an exception when no file is given"""

        with self.assertRaises(ValueError):
            grapple.bam_to_fq('')

    def test_none_file(self):
        """Should raise an exception when None is given"""

        with self.assertRaises(AttributeError):
            grapple.bam_to_fq(None)

    def test_bad_prefix(self):
        """Should raise an exception when a type that cannot be converted to a string is used as the prefix"""

        with self.assertRaises(TypeError):
            grapple.bam_to_fq(self._test_file, prefix_id=0)

    def test_none_prefix(self):
        """Should raise an exception when None is used as the prefix"""

        with self.assertRaises(TypeError):
            grapple.bam_to_fq(self._test_file, prefix_id=None)


class TestReadCorrection(TestCase):
    """Tests involving read_correction()"""

    def setUp(self):
        """Setup code before tests"""

        # Available test file
        self._test_file = os.path.join('test_files', 'lambda_reads.fq')

    def test_valid_file(self):
        """Should not raise an exception when a valid FASTQ file is used"""

        try:
            grapple.read_correction(self._test_file)

        except Exception as e:
            self.fail(e)

    def test_invalid_file(self):
        """Should raise an exception when a file with the wrong format is used"""

        with self.assertRaises(ValueError):
            grapple.read_correction(os.path.join('test_files', 'lambda_ref.fa'))

    def test_absent_file(self):
        """Should raise an exception when the file is not found"""

        with self.assertRaises(ValueError):
            grapple.read_correction('this_file_does_not_exit.fq')

    def test_none_file(self):
        """Should raise an exception when None is passed as for the file"""

        with self.assertRaises(AttributeError):
            grapple.read_correction(None)

    def test_invalid_cell_type(self):
        """Should raise an exception when an invalid cell type is given"""

        with self.assertRaises(ValueError):
            grapple.read_correction(self._test_file, cell_type='not_a_cell_type')

    def test_none_cell_type(self):
        """Should raise an exception when None is passed for the cell type"""

        with self.assertRaises(TypeError):
            grapple.read_correction(self._test_file, cell_type=None)

    def test_invalid_match_type(self):
        """Should raise an exception when an invalid match type is given"""

        with self.assertRaises(ValueError):
            grapple.read_correction(self._test_file, match_type='not_a_match_type')

    def test_none_match_type(self):
        """Should raise an exception when None is given as a match type"""

        with self.assertRaises(TypeError):
            grapple.read_correction(self._test_file, match_type=None)


class TestAlignment(TestCase):
    """Unit tests for read_alignment()"""

    def setUp(self):
        """Setup code for tests"""

        # Available test file
        self._test_file = os.path.join('test_files', 'lambda_reads.fq')

        # Available reference file
        self._ref_file = os.path.join('test_files', 'lambda_ref.fa')

    def test_valid_files(self):
        """Should not raise an exception when valid files are given as input"""

        try:
            grapple.read_alignment(self._test_file, self._ref_file)

        except Exception as e:
            self.fail(e)

    def test_invalid_read_file(self):
        """Should raise an exception when the read file is formatted wrong"""

        with self.assertRaises(ValueError):
            grapple.read_alignment(self._ref_file, self._ref_file)

    def test_absent_read_file(self):
        """Should raise an exception when the read file is absent"""

        with self.assertRaises(CalledProcessError):
            grapple.read_alignment('this_file_does_not_exist.fq', self._ref_file)

    def test_none_read_file(self):
        """Should raise an exception when None is passed as the read_file"""

        with self.assertRaises(AttributeError):
            grapple.read_alignment(None, self._ref_file)

    def test_invalid_ref_file(self):
        """Should raise an exception when the reference file is in the wrong format"""

        with self.assertRaises(ValueError):
            grapple.read_alignment(self._test_file, self._test_file)

    def test_absent_ref_file(self):
        """Should raise an exception when the reference file doesn't exist"""

        with self.assertRaises(CalledProcessError):
            grapple.read_alignment(self._test_file, 'this_file_does_not_exist.fa')

    def test_none_ref_file(self):
        """Should raise an exception when None is given as the reference file"""

        with self.assertRaises(AttributeError):
            grapple.read_alignment(self._test_file, None)

    def test_bad_prefix(self):
        """Should raise an exception when a type that cannot be converted to a string is given for the prefix"""

        with self.assertRaises(TypeError):
            grapple.read_alignment(self._test_file, self._ref_file, prefix_id=0)

    def test_none_prefix(self):
        """Should raise an exception when None is used as the prefix"""

        with self.assertRaises(TypeError):
            grapple.read_alignment(self._test_file, self._ref_file, prefix_id=None)


class TestSamToBam(TestCase):
    """Test cases for sam_to_bam() conversion"""

    def setUp(self):
        """Setup code for test cases"""

        # Available test files
        self._test_file = os.path.join('test_files', 'aligned_lambda.sam')

    def test_valid_file(self):
        """Should not raise an exception when a valid file is to be converted"""

        try:
            grapple.sam_to_bam(self._test_file)

        except Exception as e:
            self.fail(e)

    def test_invalid_file(self):
        """Should raise an exception when a file with the wrong format is converted"""

        with self.assertRaises(ValueError):
            grapple.sam_to_bam(os.path.join('test_files', 'lambda_ref.fa'))

    def test_absent_file(self):
        """Should raise an exception when the read file does not exist"""

        with self.assertRaises(CalledProcessError):
            grapple.sam_to_bam('this_file_does_not_exist.sam')

    def test_none_file(self):
        """Should raise an exception when None is passed as the read file"""

        with self.assertRaises(AttributeError):
            grapple.sam_to_bam(None)

    def test_bad_prefix(self):
        """Should raise an exception when a type that cannot be converted to a string is given for the prefix"""

        with self.assertRaises(TypeError):
            grapple.sam_to_bam(self._test_file, prefix_id=0)

    def test_none_prefix(self):
        """Should raise an exception when None is used as the prefix"""

        with self.assertRaises(TypeError):
            grapple.sam_to_bam(self._test_file, prefix_id=None)


class TestSortAndIndex(TestCase):
    """Test cases for sort_and_index()"""

    def setUp(self):
        """Setup code for test cases"""

        # Available test files
        self._test_file = os.path.join('test_files', 'aligned_lambda.bam')

    def test_valid_file(self):
        """Should not raise an exception when a valid file is provided"""

        try:
            grapple.sort_and_index(self._test_file)

        except Exception as e:
            self.fail(e)

    def test_invalid_file(self):
        """Should raise an exception when a file with the wrong format is used"""

        with self.assertRaises(ValueError):
            grapple.sort_and_index(os.path.join('test_files', 'lambda_ref.fa'))

    def test_absent_file(self):
        """Should raise an exception when a read file does not exist"""

        with self.assertRaises(CalledProcessError):
            grapple.sort_and_index('this_file_does_not_exist.bam')

    def test_none_file(self):
        """Should raise an exception when None is passed as the read file"""

        with self.assertRaises(AttributeError):
            grapple.sort_and_index(None)

    def test_bad_prefix(self):
        """Should raise an exception when a type that cannot be converted to a string is given for the prefix"""

        with self.assertRaises(TypeError):
            grapple.sort_and_index(self._test_file, prefix_id=0)

    def test_none_prefix(self):
        """Should raise an exception when None is used as the prefix"""

        with self.assertRaises(TypeError):
            grapple.sort_and_index(self._test_file, prefix_id=None)


class TestCallVariants(TestCase):
    """Test cases for call_variants()"""

    def setUp(self):
        """Setup code for test cases"""

        # Available test file
        self._test_file = os.path.join('test_files', 'sorted_lambda.bam')

        # Available reference file
        self._ref_file = os.path.join('test_files', 'lambda_ref.fa')

    def test_valid_files(self):
        """Should not raise an exception when supplying valid files"""

        try:
            grapple.call_variants(self._test_file, self._ref_file)

        except Exception as e:
            self.fail(e)

    def test_invalid_read_file(self):
        """Should raise an exception when supplying a read file with the wrong format"""

        with self.assertRaises(ValueError):
            grapple.call_variants(self._ref_file, self._ref_file)

    def test_absent_read_file(self):
        """Should raise an exception when supplying a read file that doesn't exist"""

        with self.assertRaises(CalledProcessError):
            grapple.call_variants('this_file_does_not_exist.bam', self._ref_file)

    def test_none_read_file(self):
        """Should raise an exception when None is used as the read file"""

        with self.assertRaises(AttributeError):
            grapple.call_variants(None, self._ref_file)

    def test_invalid_ref_file(self):
        """Should raise an exception when the reference file is not in the right format"""

        with self.assertRaises(ValueError):
            grapple.call_variants(self._test_file, os.path.join('test_files', 'lambda_reads.fq'))

    def test_absent_ref_file(self):
        """Should raise an exception when the reference file doesn't exist"""

        with self.assertRaises(CalledProcessError):
            grapple.call_variants(self._test_file, 'this_file_does_not_exist.fa')

    def test_none_ref_file(self):
        """Should raise an exception when None is passed as the reference file"""

        with self.assertRaises(AttributeError):
            grapple.call_variants(self._test_file, None)

    def test_bad_prefix(self):
        """Should raise an exception when a type that cannot be converted to a string is given for the prefix"""

        with self.assertRaises(TypeError):
            grapple.call_variants(self._test_file, self._ref_file, prefix_id=0)

    def test_none_prefix(self):
        """Should raise an exception when None is used as the prefix"""

        with self.assertRaises(TypeError):
            grapple.call_variants(self._test_file, self._ref_file, prefix_id=None)


class TestFormatConsensus(TestCase):
    """Test cases for format_consensus()"""

    def setUp(self):
        """Setup for test cases"""

        # Available test file
        self._test_file = os.path.join('test_files', 'lambda_consensus.fa')

    def test_valid_file(self):
        """Should not raise an exception when supplying valid files"""

        try:
            grapple.format_consensus(self._test_file)

        except Exception as e:
            self.fail(e)

    def test_invalid_file(self):
        """Should raise an exception when the file is not in FASTA format"""

        with self.assertRaises(ValueError):
            grapple.format_consensus(os.path.join('test_files', 'lambda_reads.fq'))

    def test_absent_file(self):
        """Should raise an exception when the file does not exist"""

        with self.assertRaises(IOError):
            grapple.format_consensus('this_file_does_not_exist.fa')

    def test_none_file(self):
        """Should raise an exception when None is passed as a file"""

        with self.assertRaises(AttributeError):
            grapple.format_consensus(None)

    def test_bad_prefix(self):
        """Should raise an exception when a type that cannot be converted to a string is given for the prefix"""

        with self.assertRaises(TypeError):
            grapple.format_consensus(self._test_file, prefix_id=0)

    def test_none_prefix(self):
        """Should raise an exception when None is used as the prefix"""

        with self.assertRaises(TypeError):
            grapple.format_consensus(self._test_file, prefix_id=None)


if __name__ == '__main__':
    unittest.main()
