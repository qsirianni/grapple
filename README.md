IPA
===

Iontorrent Pipeline Assembler
-----------------------------

IPA is a small Python script intended to allow users to create consensus genomes from NGS reads. IPA utilizes common
bioinformatics tools such as samtools and bowtie2 in order to pipeline the inputted NGS reads into a consensus FASTA file.

Currently, the script is focused on Iontorrent NGS reads. However, it has the potential to accept other forms of NGS reads
with a bit more development.

Using IPA
---------

To run IPA, clone this respository and run the main script *ipa.py* with a Python interpreter 
(it should work on version 2.7 and greater). You will need to install the dependency *psutil* from 
PyPI before the script can be used. Try the command:

    sudo pip install psutil
  
in order to install the dependency.
