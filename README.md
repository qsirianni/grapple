Grapple
=======

Genome Reference Assembly Pipeline
----------------------------------

Grapple is a small Python script intended to allow users to create consensus genomes from unpaired NGS reads on 
POSIX platforms. Grapple utilizes common bioinformatics tools such as samtools and bowtie2 in order to pipeline the 
inputted NGS reads into a consensus FASTA file.

Using Grapple
-------------

To run Grapple, clone this repository and run the main script *grapple.py* with a Python interpreter
(it was designed to work on version 2.7). You will need to install Grapple's dependencies from
PyPI before the script can be used. Try the command:

    sudo pip install -r requirements.txt

Grapple also requires the following utilities to be installed and in your PATH:

* samtools
* karect
* bcftools
* bowtie2

[![Build Status](https://travis-ci.org/qsirianni/ipa.svg?branch=master)](https://travis-ci.org/qsirianni/ipa)
