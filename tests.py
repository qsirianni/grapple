#!/usr/bin/python

"""
Contains unit tests for the IPA script

Author: Quinton Sirianni
"""

import filecmp
import os
import os.path
import unittest
from subprocess import CalledProcessError
from unittest import TestCase

import ipa


class TestCheckEnv(TestCase):
    """Tests involving the check_env function"""

    def test_check_utils(self):
        # Util exists
        try:
            ipa.check_env(['echo'])

        except CalledProcessError:
            self.fail('check_env() raised a CalledProcessError unexpectedly')

        # Util doesn't exist
        with self.assertRaises(CalledProcessError):
            ipa.check_env(['this_utility_should_not_exist'])

    def test_check_OS(self):
        # Supported environment
        if os.name == 'posix':
            try:
                ipa.check_env([])

            except OSError:
                self.fail('check_env() raised an OSError unexpectedly')

        # Unsupported environment
        else:
            with self.assertRaises(OSError):
                ipa.check_env([])


class TestBamToFq(TestCase):
    """Tests involving the bam_to_fq function"""

    def test_conversion(self):
        try:
            raw_reads = os.path.join('test_files', 'lambda_iontorrent.bam')
            ref_reads = os.path.join('test_files', 'lambda_reads.fq')
            converted_reads = ipa.bam_to_fq(raw_reads)

            # Ensure conversion is correct
            self.assertTrue(filecmp.cmp(ref_reads, converted_reads))

        except CalledProcessError:
            self.fail('bam_to_fq() raised a CalledProcessError unexpectedly')


if __name__ == '__main__':
    unittest.main()
