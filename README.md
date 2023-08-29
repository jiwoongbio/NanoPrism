## Requirements

1. Perl - https://www.perl.org
2. DIAMOND - https://github.com/bbuchfink/diamond
3. minimap2 - https://github.com/lh3/minimap2
4. Linux commands: sort, gzip, wget - https://www.gnu.org/software/wget/


## Install

If you already have Git (https://git-scm.com) installed, you can get the latest development version using Git.

```
git clone https://github.com/jiwoongbio/NanoPrism.git
```


## Usage

```
Usage:   perl NanoPrism_data.pl [options] [CDS.fasta [read.fastq ...]]

Options: -h       display this help message
         -r       redownload data
         -p INT   number of threads [1]
         -t STR   target NCBI taxonomy IDs or file
         -x STR   minimap2 preset [map-ont]
         -c FLOAT minimum coverage [0.9]
         -P FILE  minimap2 PAF file
         -S STR   stranded, "f" or "r"
         -M DIR   MetaPrism directory [MetaPrism]
```


## Usage example

1. Create a non-redundant coding sequence database file using the given NCBI taxonomy IDs as input.

```
threads=8
taxonomy_id=287
cds_fasta_file=Pseudomonas_aeruginosa.CDS.fasta

perl NanoPrism.pl -p $threads -t $taxonomy_id $cds_fasta_file
```

2. Align sequencing reads to the coding sequences and determine the abundance of genes.

```
threads=8
taxonomy_id=287
cds_fasta_file=Pseudomonas_aeruginosa.CDS.fasta
read_fastq_file=Pseudomonas_aeruginosa.nanopore_reads.fastq.gz

perl NanoPrism.pl -p $threads -t $taxonomy_id $cds_fasta_file $read_fastq_file
```
