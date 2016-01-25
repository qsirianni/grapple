#!/usr/bin/env python

"""Setup script for distribution."""

from distutils.core import setup


if __name__ == '__main__':
    setup(
        name='grapple',
        description='Genome Reference Assembly Pipeline',
        author='Quinton Sirianni',
        version='0.2.3',
        license='MIT',
        scripts=['grapple.py'])