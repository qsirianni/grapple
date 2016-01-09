Grapple [![Build Status](https://travis-ci.org/qsirianni/grapple.svg)](https://travis-ci.org/qsirianni/grapple)
===============================================================================================================

Genome Reference Assembly Pipeline
----------------------------------

Grapple is a small Python script intended to allow users to create consensus genomes from unpaired NGS reads.
Grapple utilizes common bioinformatics tools such as samtools and bowtie2 in order to pipeline the inputted NGS reads
into a consensus FASTA file.

Using Grapple
-------------

To run Grapple, clone this repository and run the main script *grapple.py* with a Python runtime
(the script has been successfully executed on both CPython version 2.7 as well as CPython versions 3.2 and greater).
You will need to install Grapple's dependencies from PyPI before the script can be used. Try the command:

    sudo pip install -r requirements.txt

Grapple also requires the following utilities to be installed and in your PATH:

* samtools
* karect
* bowtie2
* bcftools

Homebrew Formula
----------------

If you would like to automate the installation of Grapple, a homebrew formula is available. Simply enter the following 
commands:

    brew tap qsirianni/bioinformatics
    brew install grapple
