## Requirements

1. Perl - https://www.perl.org
2. Perl module "Bio::DB::Taxonomy" - https://bioperl.org
3. DIAMOND - https://github.com/bbuchfink/diamond
4. minimap2 - https://github.com/lh3/minimap2
5. MetaPrism - https://github.com/jiwoongbio/MetaPrism
6. Linux commands: sort, gzip, wget - https://www.gnu.org/software/wget/
7. Kraken 2 - https://github.com/DerrickWood/kraken2 (You can use another taxonomic identification tool.)


## Install

If you already have Git (https://git-scm.com) installed, you can get the latest development version using Git.

```
git clone https://github.com/jiwoongbio/NanoPrism.git

cd NanoPrism
git clone https://github.com/jiwoongbio/MetaPrism.git
```


## Usage

```
Usage:   perl NanoPrism_data.pl [options] [CDS.fasta [read.fastq ...]]

Options: -h       display this help message
         -r       redownload data
         -p INT   number of CPU threads [1]
         -t STR   comma-separated target NCBI taxonomy IDs or files
         -k STR   comma-separated Kraken 2 report files
         -x STR   minimap2 preset [map-ont]
         -c FLOAT minimum coverage for redundancy [0.9]
         -i FLOAT minimum identity for read mapping [0.5]
         -b STR   base abundance orthologies [K02950,K02874,K02946,K02948,K02867,K02952,K02886,K02988,K02992,K02965]
         -P FILE  minimap2 PAF file
         -S STR   stranded, "f" or "r"
         -M DIR   MetaPrism directory [MetaPrism]
```


## Example: Analysis of single species

1. Create a non-redundant coding sequence database file using the given NCBI taxonomy IDs as input.

```
threads=8
taxonomy_id=287
cds_fasta_file=Pseudomonas_aeruginosa.CDS.fasta # output file

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


## Example: Analysis of metagenomes

1. Download Zymo sequencing reads

```
wget --no-verbose ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR315/004/ERR3152364/ERR3152364.fastq.gz
wget --no-verbose ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR315/006/ERR3152366/ERR3152366.fastq.gz
```

2. Taxonomic classification by Kraken 2

```
export KRAKEN2=kraken2                 # Kraken 2 executable
export KRAKEN2_DB=kraken/k2_pluspf     # Kraken 2 database directory

# Only the Kraken 2 report is used by NanoPrism, so read-level output is discarded.
time $KRAKEN2 --db $KRAKEN2_DB --threads $threads --report ERR3152364.kraken2.report.txt --gzip-compressed ERR3152364.fastq.gz > /dev/null
time $KRAKEN2 --db $KRAKEN2_DB --threads $threads --report ERR3152366.kraken2.report.txt --gzip-compressed ERR3152366.fastq.gz > /dev/null
```

3. Prepare coding sequences using Kraken 2 taxonomic profiles

```
time perl NanoPrism.pl -p $threads -k ERR3152364.kraken2.report.txt -k ERR3152366.kraken2.report.txt NanoPrism.fasta # NanoPrism.fasta is the output file.
```

4. Profile gene abundances

```
time perl NanoPrism.pl -p $threads NanoPrism.fasta ERR3152364.fastq.gz > ERR3152364.NanoPrism.abundance.txt
time perl NanoPrism.pl -p $threads NanoPrism.fasta ERR3152366.fastq.gz > ERR3152366.NanoPrism.abundance.txt
```


By default, the abundances are normalized by the following genes:

| Orthology | Definition |
| --- | --- |
| K02950 | RP-S12, MRPS12, rpsL; small subunit ribosomal protein S12 |
| K02874 | RP-L14, MRPL14, rplN; large subunit ribosomal protein L14 |
| K02946 | RP-S10, MRPS10, rpsJ; small subunit ribosomal protein S10 |
| K02948 | RP-S11, MRPS11, rpsK; small subunit ribosomal protein S11 |
| K02867 | RP-L11, MRPL11, rplK; large subunit ribosomal protein L11 |
| K02952 | RP-S13, rpsM; small subunit ribosomal protein S13 |
| K02886 | RP-L2, MRPL2, rplB; large subunit ribosomal protein L2 |
| K02988 | RP-S5, MRPS5, rpsE; small subunit ribosomal protein S5 |
| K02992 | RP-S7, MRPS7, rpsG; small subunit ribosomal protein S7 |
| K02965 | RP-S19, rpsS; small subunit ribosomal protein S19 |
