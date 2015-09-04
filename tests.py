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

    def tearDown(self):
        """Tear down code for tests"""

        # Operating system type needs to be restored to os module
        os.name = self._os_name

    def test_utils_exist(self):
        """Should not raise an error when checking for the echo utility"""

        try:
            os.name = 'posix'
            ipa.check_env(['echo'])

        except CalledProcessError:
            self.fail('check_env() raised a CalledProcessError unexpectedly')

    def test_utils_absent(self):
        """Should raise an error when checking for a non-existent utility"""

        os.name = 'posix'

        with self.assertRaises(CalledProcessError):
            ipa.check_env(['this_utility_should_not_exist'])

    def test_env_is_valid(self):
        """Should not raise an exception when environment is approved"""

        try:
            os.name = 'posix'
            ipa.check_env(['echo'])

        except OSError:
            self.fail('check_env() raised an OSError unexpectedly')

    def test_env_is_bad(self):
        """Should raise an exception when environment is not approved"""

        os.name = 'nt'
        with self.assertRaises(OSError):
            ipa.check_env(['echo'])


if __name__ == '__main__':
    unittest.main()
