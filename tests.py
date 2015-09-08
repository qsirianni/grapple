#!/usr/bin/python

"""
Contains unit tests for the IPA script

Author: Quinton Sirianni
"""

import os
import os.path
import unittest
from subprocess import CalledProcessError
from unittest import TestCase

import ipa


class TestCheckEnv(TestCase):
    """Tests involving check_env()"""

    def setUp(self):
        """Setup code for the tests"""

        # Operating system name needs to be recorded from os module
        self._os_name = os.name

    def test_utils_exist(self):
        """Should not raise an exception when checking for the echo utility"""

        try:
            os.name = 'posix'
            ipa.check_env(['echo'])

        except Exception as e:
            self.fail(e.message)

    def test_utils_absent(self):
        """Should raise an exception when checking for a non-existent utility"""

        os.name = 'posix'

        with self.assertRaises(CalledProcessError):
            ipa.check_env(['this_utility_does_not_exist'])

    def test_utils_none(self):
        """Should raise an exception when given None as the required utilities"""

        os.name ='posix'

        with self.assertRaises(TypeError):
            ipa.check_env(None)

    def test_utils_empty(self):
        "Should not raise an exception when given an empty list of required utilities"

        try:
            os.name = 'posix'
            ipa.check_env([])

        except CalledProcessError:
            self.fail('check_env() raised a CalledProcessError unexpectedly')

    def test_env_is_bad(self):
        """Should raise an exception when environment is not approved"""

        os.name = 'bad_env'
        with self.assertRaises(OSError):
            ipa.check_env(['echo'])

    def tearDown(self):
        """Tear down code for tests"""

        # Operating system type needs to be restored to os module
        os.name = self._os_name


class TestBamToFq(TestCase):
    """Tests involving bam_to_fq()"""

    def test_valid_file(self):
        """Should not raise an exception when a valid bam file is used"""

        try:
            ipa.bam_to_fq(os.path.join('test_files', 'lambda_iontorrent.bam'))

        except Exception as e:
            self.fail(e.message)

    def test_invalid_file(self):
        """Should raise an exception when the wrong type of file is used"""

        with self.assertRaises(ValueError):
            ipa.bam_to_fq(os.path.join('test_files', 'lambda_ref.fa'))

    def test_absent_file(self):
        """Should raise an exception when the file does not exist"""

        with self.assertRaises(CalledProcessError):
            ipa.bam_to_fq('this_file_does_not_exist.bam')

    def test_no_file(self):
        """Should raise an exception when no file is given"""

        with self.assertRaises(ValueError):
            ipa.bam_to_fq('')

    def test_none_file(self):
        """Should raise an exception when None is given"""

        with self.assertRaises(AttributeError):
            ipa.bam_to_fq(None)


class TestReadCorrection(TestCase):
    """Tests involving read_correction()"""

    def setUp(self):
        """Setup code before tests"""

        # Available test file
        self._test_file = os.path.join('test_files', 'lambda_reads.fq')

    def test_valid_file(self):
        """Should not raise an exception when a valid FASTQ file is used"""

        try:
            ipa.read_correction(self._test_file)

        except Exception as e:
            self.fail(e.message)

    def test_invalid_file(self):
        """Should raise an exception when a file with the wrong format is used"""

        with self.assertRaises(ValueError):
            ipa.read_correction(os.path.join('test_files', 'lambda_ref.fa'))

    def test_absent_file(self):
        """Should raise an exception when the file is not found"""

        with self.assertRaises(ValueError):
            ipa.read_correction('this_file_does_not_exit.fq')

    def test_none_file(self):
        """Should raise an exception when None is passed as for the file"""

        with self.assertRaises(AttributeError):
            ipa.read_correction(None)

    def test_invalid_cell_type(self):
        """Should raise an exception when an invalid cell type is given"""

        with self.assertRaises(ValueError):
            ipa.read_correction(self._test_file, cell_type='not_a_cell_type')

    def test_none_cell_type(self):
        """Should raise an exception when None is passed for the cell type"""

        with self.assertRaises(TypeError):
            ipa.read_correction(self._test_file, cell_type=None)

    def test_invalid_match_type(self):
        """Should raise an exception when an invalid match type is given"""

        with self.assertRaises(ValueError):
            ipa.read_correction(self._test_file, match_type='not_a_match_type')

    def test_none_match_type(self):
        """Should raise an exception when None is given as a match type"""

        with self.assertRaises(TypeError):
            ipa.read_correction(self._test_file, match_type=None)


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
            ipa.read_alignment(self._test_file, self._ref_file)

        except Exception as e:
            self.fail(e.message)

    def test_invalid_read_file(self):
        """Should raise an exception when the read file is formatted wrong"""

        with self.assertRaises(ValueError):
            ipa.read_alignment(self._ref_file, self._ref_file)

    def test_absent_read_file(self):
        """Should raise an exception when the read file is absent"""

        with self.assertRaises(CalledProcessError):
            ipa.read_alignment('this_file_does_not_exist.fq', self._ref_file)

    def test_none_read_file(self):
        """Should raise an exception when None is passed as the read_file"""

        with self.assertRaises(AttributeError):
            ipa.read_alignment(None, self._ref_file)

    def test_invalid_ref_file(self):
        """Should raise an exception when the reference file is in the wrong format"""

        with self.assertRaises(ValueError):
            ipa.read_alignment(self._test_file, self._test_file)

    def test_absent_ref_file(self):
        """Should raise an exception when the reference file doesn't exist"""

        with self.assertRaises(CalledProcessError):
            ipa.read_alignment(self._test_file, 'this_file_does_not_exist.fa')

    def test_none_ref_file(self):
        """Should raise an exception when None is given as the reference file"""

        with self.assertRaises(AttributeError):
            ipa.read_alignment(self._test_file, None)


class TestSamToBam(TestCase):
    """Test cases for sam_to_bam() conversion"""

    def test_valid_file(self):
        """Should not raise an exception when a valid file is to be converted"""

        try:
            ipa.sam_to_bam(os.path.join('test_files', 'aligned_lambda.sam'))

        except Exception as e:
            self.fail(e.message)

    def test_invalid_file(self):
        """Should raise an exception when a file with the wrong format is converted"""

        with self.assertRaises(ValueError):
            ipa.sam_to_bam(os.path.join('test_files', 'lambda_ref.fa'))

    def test_absent_file(self):
        """Should raise an exception when the read file does not exist"""

        with self.assertRaises(CalledProcessError):
            ipa.sam_to_bam('this_file_does_not_exist.sam')

    def test_none_file(self):
        """Should raise an exception when None is passed as the read file"""

        with self.assertRaises(AttributeError):
            ipa.sam_to_bam(None)


class TestSortAndIndex(TestCase):
    """Test cases for sort_and_index()"""

    def test_valid_file(self):
        """Should not raise an exception when a valid file is provided"""

        try:
            ipa.sort_and_index(os.path.join('test_files', 'aligned_lambda.bam'))

        except Exception as e:
            self.fail(e.message)

    def test_invalid_file(self):
        """Should raise an exception when a file with the wrong format is used"""

        with self.assertRaises(ValueError):
            ipa.sort_and_index(os.path.join('test_files', 'lambda_ref.fa'))

    def test_absent_file(self):
        """Should raise an exception when a read file does not exist"""

        with self.assertRaises(CalledProcessError):
            ipa.sort_and_index('this_file_does_not_exist.bam')

    def test_none_file(self):
        """Should raise an exception when None is passed as the read file"""

        with self.assertRaises(AttributeError):
            ipa.sort_and_index(None)


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
            ipa.call_variants(self._test_file, self._ref_file)

        except Exception as e:
            self.fail(e.message)

    def test_invalid_read_file(self):
        """Should raise an exception when supplying a read file with the wrong format"""

        with self.assertRaises(ValueError):
            ipa.call_variants(self._ref_file, self._ref_file)

    def test_absent_read_file(self):
        """Should raise an exception when supplying a read file that doesn't exist"""

        with self.assertRaises(CalledProcessError):
            ipa.call_variants('this_file_does_not_exist.bam', self._ref_file)

    def test_none_read_file(self):
        """Should raise an exception when None is used as the read file"""

        with self.assertRaises(AttributeError):
            ipa.call_variants(None, self._ref_file)

    def test_invalid_ref_file(self):
        """Should raise an exception when the reference file is not in the right format"""

        with self.assertRaises(ValueError):
            ipa.call_variants(self._test_file, os.path.join('test_files', 'lambda_reads.fq'))

    def test_absent_ref_file(self):
        """Should raise an exception when the reference file doesn't exist"""

        with self.assertRaises(CalledProcessError):
            ipa.call_variants(self._test_file, 'this_file_does_not_exist.fa')

    def test_none_ref_file(self):
        """Should raise an exception when None is passed as the reference file"""

        with self.assertRaises(AttributeError):
            ipa.call_variants(self._test_file, None)


if __name__ == '__main__':
    unittest.main()
