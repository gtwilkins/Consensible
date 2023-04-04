# Consensible

## Overview
Consensible is an application for the analysis of sedimentary ancient DNA.

## Limitations
Consensible is designed for use with short read sequence data, namely those produced by Illumina

## Requirements
* gcc

## Installation
The install directory can be specified with the following command (if this omitted, Consensible is installed to /usr/local/bin/):

	./configure --prefix=/path/to/directory/

Consensible can be installed with the following commands:

	make
	sudo make install

## Use
When running consensible, input, output and temporary files must be specified with the following arguments:
* -i	Input shotgun sequence file(s).
* -p	Prefix for indexed shotgun sequence files. See notes for details.
* -q	Query fasta file containing one or more query sequences.
* -o	Output filename prefix.

An example command should look as follows:
	consensible -i /myinputs/project101_data.fa -p /mytempdata/project101 -q /myinput/interesting_gene.fa -o /myoutput/interesting_gene_result.fa

Note that when running successive queries on a dataset, it does not need to be re-indexed. Hence, only the prefix for the pre-existing index files need be supplied, and the raw input data can be omitted as follows:
	consensible -p /mytempdata/project101 -q /myinput/another_interesting_gene.fa -o /myoutput/another_interesting_gene_result.fa
