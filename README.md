Grapple [![Build Status](https://travis-ci.org/qsirianni/grapple.svg)](https://travis-ci.org/qsirianni/grapple)
===============================================================================================================

Genome Reference Assembly Pipeline
----------------------------------

Grapple is a small Python script intended to allow users to create consensus genomes from unpaired NGS reads.
Grapple utilizes common bioinformatics tools such as samtools and bowtie2 in order to pipeline the inputted NGS reads
into a consensus FASTA file.

Supported Platforms
-------------------

Grapple was developed using Python 2.7 on OS X. However, it will likely work with little to no modification on other
platforms with other Python runtimes (including Python 3 runtimes). The main requirements for the script to execute
are its dependencies (enumerated below).

Using Grapple
-------------

To run Grapple, clone this repository and run the main script *grapple.py*.
You will need to install Grapple's dependencies from PyPI before the script can be used. Try the command:

    sudo pip install -r requirements.txt

Grapple also requires the following utilities to be installed and in your PATH:

* samtools
* karect
* bowtie2
* bcftools

Homebrew Formula
----------------

If you would like to automate the installation of Grapple, a Homebrew formula is available. Try the following commands:

    brew tap qsirianni/bioinformatics
    brew install grapple
