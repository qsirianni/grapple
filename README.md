IPA
===

Iontorrent Pipeline Assembler
-----------------------------

IPA is a small Python script intended to allow users to create consensus genomes from NGS reads on POSIX platforms. IPA utilizes common bioinformatics tools such as samtools and bowtie2 in order to pipeline the inputted NGS reads into a consensus FASTA file.

Currently, the script is focused on Iontorrent NGS reads. However, it has the potential to accept other forms of NGS reads
with a bit more development.

Using IPA
---------

To run IPA, clone this repository and run the main script *ipa.py* with a Python interpreter
(it was designed to work on version 2.7). You will need to install IPA's dependencies from
PyPI before the script can be used. Try the command:

    sudo pip install -r requirements.txt

IPA also requires the following utilities to be installed and in your PATH:

* samtools
* karect
* bcftools
* bowtie2

[![Build Status](https://travis-ci.org/qsirianni/ipa.svg?branch=master)](https://travis-ci.org/qsirianni/ipa)
