threads=16


# Download zymo sequencing reads
wget --no-verbose ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR315/004/ERR3152364/ERR3152364.fastq.gz
wget --no-verbose ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR315/006/ERR3152366/ERR3152366.fastq.gz


# Taxonomic classification by Centrifuge

CF_DIR=centrifuge # Centrifuge directory
CF_IDX=centrifuge/data/p_compressed # Centrifuge index filename prefix

export PATH=$CF_DIR${PATH:+:${PATH}}

sample=ERR3152364
centrifuge -p $threads -x $CF_IDX -U $sample.fastq.gz --report-file $sample.centrifuge.report.txt | gzip > $sample.centrifuge.txt.gz

sample=ERR3152366
centrifuge -p $threads -x $CF_IDX -U $sample.fastq.gz --report-file $sample.centrifuge.report.txt | gzip > $sample.centrifuge.txt.gz


# Filter 1% abundant species

sample=ERR3152364
awk -F'\t' '(NR == 1 || $3 == "species")' $sample.centrifuge.report.txt > $sample.centrifuge.report.species.txt
awk -F'\t' -va=$(awk -F'\t' '(NR > 1) {a += $6} END {print a}' $sample.centrifuge.report.species.txt) '(NR == 1 || $6 / a >= 0.01)' $sample.centrifuge.report.species.txt > $sample.centrifuge.report.species.filtered.txt

sample=ERR3152366
awk -F'\t' '(NR == 1 || $3 == "species")' $sample.centrifuge.report.txt > $sample.centrifuge.report.species.txt
awk -F'\t' -va=$(awk -F'\t' '(NR > 1) {a += $6} END {print a}' $sample.centrifuge.report.species.txt) '(NR == 1 || $6 / a >= 0.01)' $sample.centrifuge.report.species.txt > $sample.centrifuge.report.species.filtered.txt

for sample in ERR3152364 ERR3152366; do awk -F'\t' '(NR > 1) {print $2}' $sample.centrifuge.report.species.filtered.txt; done | sort -nu > species.filtered.txt


# Prepare coding sequences
time perl ../NanoPrism.pl -p $threads -t species.filtered.txt NanoPrism.fasta


# Profile gene abundances
time perl ../NanoPrism.pl -p $threads NanoPrism.fasta ERR3152364.fastq.gz > ERR3152364.NanoPrism.abundance.txt
time perl ../NanoPrism.pl -p $threads NanoPrism.fasta ERR3152366.fastq.gz > ERR3152366.NanoPrism.abundance.txt
