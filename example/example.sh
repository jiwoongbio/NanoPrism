#!/bin/bash

set -euo pipefail

threads=16

# Download Zymo sequencing reads
wget --no-verbose ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR315/004/ERR3152364/ERR3152364.fastq.gz
wget --no-verbose ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR315/006/ERR3152366/ERR3152366.fastq.gz

# Taxonomic classification by Kraken 2
export KRAKEN2=kraken2                 # Kraken 2 executable
export KRAKEN2_DB=kraken/k2_pluspf     # Kraken 2 database directory

# Only the Kraken 2 report is used by NanoPrism, so read-level output is discarded.
time $KRAKEN2 --db $KRAKEN2_DB --threads $threads --report ERR3152364.kraken2.report.txt --gzip-compressed ERR3152364.fastq.gz > /dev/null
time $KRAKEN2 --db $KRAKEN2_DB --threads $threads --report ERR3152366.kraken2.report.txt --gzip-compressed ERR3152366.fastq.gz > /dev/null

# Prepare coding sequences using Kraken 2 taxonomic profiles
time perl ../NanoPrism.pl -p $threads -k ERR3152364.kraken2.report.txt -k ERR3152366.kraken2.report.txt NanoPrism.fasta

# Profile gene abundances
time perl ../NanoPrism.pl -p $threads NanoPrism.fasta ERR3152364.fastq.gz > ERR3152364.NanoPrism.abundance.txt
time perl ../NanoPrism.pl -p $threads NanoPrism.fasta ERR3152366.fastq.gz > ERR3152366.NanoPrism.abundance.txt
