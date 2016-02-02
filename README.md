# PAst

OVERVIEW

The <i>Pseudomonas aeruginosa</i> serotyper (PAst) is a command-line-tool for fast and reliable <i>in silico</i> serotyping of <i>P. aeruginosa</i> isolates, based on whole genome sequence assembly data. 

PAst is distributed as a Perl script as well as a web service tool hosted by the Center for Genomic Epidemiology at Center for Biological Sequence Analysis DTU: https://cge.cbs.dtu.dk/services/PAst-1.0/

USAGE

Program:      PAst (Pseudomonas aeruginosa serotyper)

Version:      1.0 (using NCBI blast+ 2.2.29)

Usage:        ./PAserotyper.pl <path/to/BLASTbin> <path/to/output/directory> <path/to/input/directory> <path/to/OSAdb>

BLAST bin:    Path to the bin containing blasts executables (ex. /usr/bin/ncbi-blast-2.2.29+/bin/)
Output dir:   Directory where blast reports, OSA multifasta files and result summary (serotyping.txt) are placed
Input dir:    Directory containing assembled input genomes in multifasta format (and nothing else)
OSA database: Path to OSA database file (downloadable from Github)

SETUP

Download the Perl script and OSA database file. PAst is dependent on the installation of Perl and blastn (blast+: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

DESCRIPTION OF OUTPUT FILES
