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


## Example: Analysis of single species

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


## Example: Analysis of metagenomes

1. Download zymo sequencing reads

```
wget --no-verbose ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR315/004/ERR3152364/ERR3152364.fastq.gz
wget --no-verbose ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR315/006/ERR3152366/ERR3152366.fastq.gz
```

2. Taxonomic classification by Centrifuge

```
CF_DIR=centrifuge # Centrifuge directory
CF_IDX=centrifuge/data/p_compressed # Centrifuge index filename prefix

export PATH=$CF_DIR${PATH:+:${PATH}}

sample=ERR3152364
centrifuge -p $threads -x $CF_IDX -U $sample.fastq.gz --report-file $sample.centrifuge.report.txt | gzip > $sample.centrifuge.txt.gz

sample=ERR3152366
centrifuge -p $threads -x $CF_IDX -U $sample.fastq.gz --report-file $sample.centrifuge.report.txt | gzip > $sample.centrifuge.txt.gz
```

3. Filter 1% abundant species

```
sample=ERR3152364
awk -F'\t' '(NR == 1 || $3 == "species")' $sample.centrifuge.report.txt > $sample.centrifuge.report.species.txt
awk -F'\t' -va=$(awk -F'\t' '(NR > 1) {a += $6} END {print a}' $sample.centrifuge.report.species.txt) '(NR == 1 || $6 / a >= 0.01)' $sample.centrifuge.report.species.txt > $sample.centrifuge.report.species.filtered.txt

sample=ERR3152366
awk -F'\t' '(NR == 1 || $3 == "species")' $sample.centrifuge.report.txt > $sample.centrifuge.report.species.txt
awk -F'\t' -va=$(awk -F'\t' '(NR > 1) {a += $6} END {print a}' $sample.centrifuge.report.species.txt) '(NR == 1 || $6 / a >= 0.01)' $sample.centrifuge.report.species.txt > $sample.centrifuge.report.species.filtered.txt

for sample in ERR3152364 ERR3152366; do awk -F'\t' '(NR > 1) {print $2}' $sample.centrifuge.report.species.filtered.txt; done | sort -nu > species.filtered.txt
```

4. Prepare coding sequences

```
time perl NanoPrism.pl -p $threads -t species.filtered.txt NanoPrism.fasta
```

5. Profile gene abundances

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
