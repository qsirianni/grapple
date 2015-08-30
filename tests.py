#!/usr/bin/env python

"""
Contains unit tests for the IPA script

Author: Quinton Sirianni
"""

import os
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

        else:
            with self.assertRaises(OSError):
                ipa.check_env([])


class TestBamToFq(TestCase):
    """Tests involving the bam_to_fq function"""

    def test_conversion(self):
        pass


if __name__ == '__main__':
    unittest.main()