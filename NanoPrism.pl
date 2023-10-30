#!/usr/bin/env perl
# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use IPC::Open2;
use List::Util qw(sum);
use Bio::DB::Taxonomy;
use Getopt::Long qw(:config no_ignore_case);

chomp(my $directory = `readlink -f $0`);
$directory =~ s/\/[^\/]*$//;
my $dataDirectory = "$directory/data";
system("mkdir -p $dataDirectory");

my @refseqCategoryList = ('reference genome', 'representative genome', 'na');
my @assemblyLevelList = ('Complete Genome', 'Chromosome', 'Scaffold', 'Contig');

my @assemblySummaryFilePathList = (
	'ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt',
	'ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt',
);

my $sortedAssemblySummaryFile = "$dataDirectory/assembly_summary.sorted.txt";
my $maximumNumber = 100;

my @pafMandatoryColumnList = ('query_name', 'query_length', 'query_start', 'query_end', 'strand', 'target_name', 'target_length', 'target_start', 'target_end', 'match', 'alignment_length', 'mapping_quality');

chomp(my $hostname = `hostname`);

my $baseAbundanceOrthologies = 'K02950,K02874,K02946,K02948,K02867,K02952,K02886,K02988,K02992,K02965';

my @targetTaxonomyIdList = ();
GetOptions(
	'h' => \(my $help = ''),
	'r' => \(my $redownload = ''),
	'p=i' => \(my $threads = 1),
	't=s' => \@targetTaxonomyIdList,
	'x=s' => \(my $preset = 'map-ont'),
	'c=f' => \(my $minimumCoverage = 0.9),
	'b=s' => \$baseAbundanceOrthologies,
	'P=s' => \(my $pafFile = ''),
	'S=s' => \(my $stranded = ''),
	'M=s' => \(my $MetaPrismDirectory = "$directory/MetaPrism"),
	'fastaLineLength=i' => \(my $fastaLineLength = 80),
);
if($help) {
	die <<EOF;

Usage:   perl NanoPrism_data.pl [options] [CDS.fasta [read.fastq ...]]

Options: -h       display this help message
         -r       redownload data
         -p INT   number of CPU threads [$threads]
         -t STR   comma-separated target NCBI taxonomy IDs or file
         -x STR   minimap2 preset [$preset]
         -c FLOAT minimum coverage [$minimumCoverage]
         -b STR   base abundance orthologies [$baseAbundanceOrthologies]
         -P FILE  minimap2 PAF file
         -S STR   stranded, "f" or "r"
         -M DIR   MetaPrism directory [$MetaPrismDirectory]

EOF
}
@targetTaxonomyIdList = map {split(/,/, $_)} @targetTaxonomyIdList;

my @optionList = ();
push(@optionList, "-t $threads");
push(@optionList, "-x $preset") if($preset ne '');

foreach my $file ("$MetaPrismDirectory/MetaPrism_data.gene.pl", "$MetaPrismDirectory/MetaPrism.pl") {
	die "'$file' is not readable.\n" unless(-r $file);
}
system("ln -sf $dataDirectory $MetaPrismDirectory/data") unless(-d "$MetaPrismDirectory/data");

lockFile($dataDirectory);

if(not -r "$dataDirectory/nodes.dmp" or not -r "$dataDirectory/names.dmp" or $redownload) {
	my $URL = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";
	my $file = "$dataDirectory/taxdump.tar.gz";
	system("wget --no-verbose -O $file $URL") if(not -r $file or $redownload);
	system("cd $dataDirectory; tar -zxf taxdump.tar.gz nodes.dmp");
	system("cd $dataDirectory; tar -zxf taxdump.tar.gz names.dmp");
	system("rm -f $dataDirectory/$_") foreach('nodes', 'parents', 'names2id', 'id2names');
}
my $db = Bio::DB::Taxonomy->new(-source => 'flatfile', -directory => $dataDirectory, -nodesfile => "$dataDirectory/nodes.dmp", -namesfile => "$dataDirectory/names.dmp");

my %refseqCategoryIndexHash = map {$refseqCategoryList[$_] => $_} 0 .. $#refseqCategoryList;
my %assemblyLevelIndexHash = map {$assemblyLevelList[$_] => $_} 0 .. $#assemblyLevelList;

if(not -r $sortedAssemblySummaryFile) {
	my @assemblySummaryFileList = ();
	foreach my $assemblySummaryFilePath (@assemblySummaryFilePathList) {
		(my $assemblySummaryFile = $assemblySummaryFilePath) =~ s/^.*\///;
		$assemblySummaryFile = "$dataDirectory/$assemblySummaryFile";
		system("wget --no-verbose -O $assemblySummaryFile $assemblySummaryFilePath") if(not -r $assemblySummaryFile or $redownload);
		push(@assemblySummaryFileList, $assemblySummaryFile);
	}
	my $pid = open2(my $reader, my $writer, "sort -t '\t' -k1,1n -k2,2n -k3,3n -k4,4n -k5,5nr | cut -f6-");
	my @headerLineList = ();
	my @outputColumnList = ();
	foreach my $assemblySummaryFileIndex (0 .. $#assemblySummaryFileList) {
		my $assemblySummaryFile = $assemblySummaryFileList[$assemblySummaryFileIndex];
		open(my $reader, ($assemblySummaryFile =~ /\.gz$/ ? "gzip -dc $assemblySummaryFile |" : $assemblySummaryFile)) or die "Can't read '$assemblySummaryFile': $!";
		my $line;
		chomp($line = <$reader>);
		push(@headerLineList, $line) if($assemblySummaryFileIndex == 0);
		chomp($line = <$reader>);
		push(@headerLineList, $line) if($assemblySummaryFileIndex == 0);
		$line =~ s/^# ?//;
		my @columnList = split(/\t/, $line, -1);
		@outputColumnList = @columnList if($assemblySummaryFileIndex == 0);
		while($line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{@columnList} = split(/\t/, $line, -1);
			next if($tokenHash{'ftp_path'} eq 'na');
			my $refseqCategoryIndex = $refseqCategoryIndexHash{$tokenHash{'refseq_category'}};
			my $assemblyLevelIndex = $assemblyLevelIndexHash{$tokenHash{'assembly_level'}};
			(my $assemblyAccessionNumber = $tokenHash{'assembly_accession'}) =~ s/^[^0-9]*//;
			my $assemblyAccessionNumberVersion = 0;
			if($assemblyAccessionNumber =~ s/\.([0-9]+)$//) {
				$assemblyAccessionNumberVersion = $1;
			}
			print $writer join("\t", $refseqCategoryIndex, $assemblyLevelIndex, $assemblySummaryFileIndex, $assemblyAccessionNumber, $assemblyAccessionNumberVersion, @tokenHash{@outputColumnList}), "\n";
		}
		close($reader);
	}
	close($writer);
	{
		open(my $writer, "> $sortedAssemblySummaryFile");
		foreach my $line (@headerLineList) {
			print $writer $line, "\n";
		}
		my %assemblyAccessionHash = ();
		while(my $line = <$reader>) {
			chomp($line);
			my %tokenHash = ();
			@tokenHash{@outputColumnList} = split(/\t/, $line, -1);
			unless($assemblyAccessionHash{$tokenHash{'assembly_accession'}}) {
				print $writer join("\t", @tokenHash{@outputColumnList}), "\n";
			}
			$assemblyAccessionHash{$tokenHash{'assembly_accession'}} = 1;
			$assemblyAccessionHash{$tokenHash{'gbrs_paired_asm'}} = 1 if($tokenHash{'paired_asm_comp'} eq 'identical');
		}
		close($writer);
	}
	close($reader);
	waitpid($pid, 0);
}

MetaPrism_data();

system("mkdir -p $dataDirectory/assembly.protein.MetaPrism");
system("mkdir -p $dataDirectory/assembly.CDS");

unlockFile($dataDirectory);

my ($fastaFile, @fastqFileList) = @ARGV;

my @cdsFastaFileList = ();
if(@targetTaxonomyIdList) {
	my %targetTaxonomyIdCountHash = ();
	foreach my $targetTaxonomyId (@targetTaxonomyIdList) {
		if(-r $targetTaxonomyId) {
			chomp(my @taxonomyIdList = `cat $targetTaxonomyId`);
			$targetTaxonomyIdCountHash{$_} = 0 foreach(@taxonomyIdList);
		} else {
			$targetTaxonomyIdCountHash{$targetTaxonomyId} = 0;
		}
	}
	my %taxonomyIdTargetTaxonomyIdHash = ();
	foreach my $targetTaxonomyId (keys %targetTaxonomyIdCountHash) {
		if(my $taxon = $db->get_taxon(-taxonid => $targetTaxonomyId)) {
			$taxonomyIdTargetTaxonomyIdHash{$_}->{$targetTaxonomyId} = 1 foreach($targetTaxonomyId, map {$_->id} $db->get_all_Descendents($taxon));
		}
	}
	open(my $reader, $sortedAssemblySummaryFile) or die "Can't read '$sortedAssemblySummaryFile': $!";
	my $line;
	chomp($line = <$reader>);
	chomp($line = <$reader>);
	$line =~ s/^# ?//;
	my @columnList = split(/\t/, $line, -1);
	while($line = <$reader>) {
		chomp($line);
		my %tokenHash = ();
		@tokenHash{@columnList} = split(/\t/, $line, -1);

		my $taxonomyId = $tokenHash{'taxid'};
		next unless(defined($taxonomyIdTargetTaxonomyIdHash{$taxonomyId}));

		my @targetTaxonomyIdList = keys %{$taxonomyIdTargetTaxonomyIdHash{$taxonomyId}};
		next unless(grep {$targetTaxonomyIdCountHash{$_} < $maximumNumber} @targetTaxonomyIdList);

		my $taxonomyIdLineage = getTaxonomyIdLineage($taxonomyId);

		my $accession = $tokenHash{'assembly_accession'};
		my $ftpPath = $tokenHash{'ftp_path'};
		(my $prefix = $ftpPath) =~ s/^.*\///;

		my $proteinMetaPrismFile = "$dataDirectory/assembly.protein.MetaPrism/$accession.txt";
		if(not -r $proteinMetaPrismFile or $redownload) {
			lockFile($proteinMetaPrismFile);
			system("(wget --no-verbose -O - https://cdc.biohpc.swmed.edu/NanoPrism/data/assembly.protein.MetaPrism/$accession.txt.gz | gzip -d > $proteinMetaPrismFile) 2> /dev/null");
			system("wget --no-verbose -O - $ftpPath/${prefix}_protein.faa.gz | gzip -d | perl $MetaPrismDirectory/MetaPrism.pl -p $threads -P - > $proteinMetaPrismFile") if(-z $proteinMetaPrismFile);
			unlockFile($proteinMetaPrismFile);
		}

		my $cdsFastaFile = "$dataDirectory/assembly.CDS/$accession.fasta";
		if(not -r $cdsFastaFile or $redownload) {
			lockFile($cdsFastaFile);
			my %proteinOrthologyHash = ();
			{
				open(my $reader, $proteinMetaPrismFile) or die "Can't read '$proteinMetaPrismFile': $!";
				my $line;
				chomp($line = <$reader>);
				my @columnList = split(/\t/, $line, -1);
				while($line = <$reader>) {
					chomp($line);
					my %tokenHash = ();
					@tokenHash{@columnList} = split(/\t/, $line, -1);
					my $proteinId = $tokenHash{'input'};
					my $orthology = $tokenHash{'gene'};
					$proteinOrthologyHash{$proteinId} = $orthology;
				}
				close($reader);
			}
			my %chromosomeSequenceHash = ();
			{
				my $chromosome = '';
				open(my $reader, "wget --no-verbose -O - $ftpPath/$prefix\_genomic.fna.gz | gzip -d |");
				while(my $line = <$reader>) {
					chomp($line);
					next if($line =~ /^>(\S*)/ && ($chromosome = $1));
					$chromosomeSequenceHash{$chromosome} .= $line;
				}
				close($reader);
			}
			my %nameSequenceHash = ();
			{
				open(my $reader, "wget --no-verbose -O - $ftpPath/$prefix\_genomic.gff.gz | gzip -d |");
				while(my $line = <$reader>) {
					chomp($line);
					next if($line =~ /^#/);
					my %tokenHash = ();
					@tokenHash{'chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'} = split(/\t/, $line, -1);
					my %attributeHash = ();
					$attributeHash{$1} = $2 while($tokenHash{'attribute'} =~ m/([^;= ]+)=([^;]+)(;|$)/g);
					$attributeHash{$1} = $2 while($tokenHash{'attribute'} =~ m/([^;" ]+) +"([^;"]+)"(;|$)/g);
					if($tokenHash{'feature'} eq 'CDS' && defined(my $proteinId = $attributeHash{'protein_id'})) {
						my $locusTag = $attributeHash{'locus_tag'};
						$locusTag = '' unless(defined($locusTag));
						my $name = join('|', $taxonomyIdLineage, $accession, $locusTag, $proteinId);
						if(defined(my $orthology = $proteinOrthologyHash{$proteinId})) {
							$name = "$name=>$orthology";
						}
						my ($chromosome, $start, $end, $strand) = @tokenHash{'chromosome', 'start', 'end', 'strand'};
						my $chromosomeSequence = $chromosomeSequenceHash{$chromosome};
						my $chromosomeLength = length($chromosomeSequence);
						my $sequence = '';
						if($start > $chromosomeLength) {
							$start = $start - $chromosomeLength;
							$end = $end - $chromosomeLength;
						}
						if($end > $chromosomeLength) {
							$sequence .= uc(substr($chromosomeSequence, $start - 1));
							$sequence .= uc(substr($chromosomeSequence, 0, $end - $chromosomeLength));
						} else {
							$sequence .= uc(substr($chromosomeSequence, $start - 1, $end - $start + 1));
						}
						$nameSequenceHash{$name} = '' unless(defined($nameSequenceHash{$name}));
						$nameSequenceHash{$name} = $nameSequenceHash{$name} . $sequence if($strand eq '+');
						$nameSequenceHash{$name} = getReverseComplementarySequence($sequence) . $nameSequenceHash{$name} if($strand eq '-');
					}
				}
				close($reader);
			}
			open(my $writer, "> $cdsFastaFile");
			foreach my $name (sort keys %nameSequenceHash) {
				my $sequence = $nameSequenceHash{$name};
				print $writer ">$name\n";
				for(my $index = 0; $index < length($sequence); $index += $fastaLineLength) {
					print $writer substr($sequence, $index, $fastaLineLength), "\n";
				}
			}
			close($writer);
			unlockFile($cdsFastaFile);
		}
		$targetTaxonomyIdCountHash{$_} += 1 foreach(@targetTaxonomyIdList);
		push(@cdsFastaFileList, $cdsFastaFile);
	}
	close($reader);
}

if(defined($fastaFile)) {
	foreach my $cdsFastaFile (@cdsFastaFileList) {
		my %excludeNameHash = ();
		my %excludeReferenceNameHash = ();
		my %nameUpdatedHash = ();
		if(-r $fastaFile) {
			open(my $reader, "minimap2 @optionList $fastaFile $cdsFastaFile |") or die "minimap2 error: $!";
			while(my $line = <$reader>) {
				chomp($line);
				next if($line =~ /^@/);
				my %tokenHash = ();
				(@tokenHash{@pafMandatoryColumnList}, my @tagTypeValueList) = split(/\t/, $line, -1);
				$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
				next unless($tokenHash{'strand'} eq '+');
				if($tokenHash{'match'} / $tokenHash{'query_length'} >= $minimumCoverage) {
					if($tokenHash{'query_name'} =~ /=>(.*)$/) {
						my $orthology = $1;
						if($tokenHash{'target_name'} =~ /=>(.*)$/) {
							if($1 eq $orthology) {
								$excludeNameHash{$tokenHash{'query_name'}} = 1;
							} else {
								print STDERR "different orthology $tokenHash{'target_name'} $tokenHash{'query_name'}\n";
							}
							my ($taxonomyIdLineage1, $accession1, $locusTag1, $proteinId1) = split(/\|/, $tokenHash{'target_name'});
							my ($taxonomyIdLineage2, $accession2, $locusTag2, $proteinId2) = split(/\|/, $tokenHash{'query_name'});
							$nameUpdatedHash{$tokenHash{'target_name'}} = join('|', getCommonTaxonomyIdLineage($taxonomyIdLineage1, $taxonomyIdLineage2), $accession1, $locusTag1, $proteinId1);
						} else {
							print STDERR "exclude $tokenHash{'target_name'}\n";
							$excludeReferenceNameHash{$tokenHash{'target_name'}} = 1;
						}
					} else {
						$excludeNameHash{$tokenHash{'query_name'}} = 1;
					}
				}
			}
			close($reader);
		}
		if(%excludeReferenceNameHash) {
			open(my $writer, "> $fastaFile.updated") or die "Can't write '$fastaFile.updated': $!";
			open(my $reader, "$fastaFile") or die "Can't read '$fastaFile': $!";
			my $add = '';
			while(my $line = <$reader>) {
				chomp($line);
				if($line =~ /^>(\S*)/) {
					my $name = $1;
					if($excludeReferenceNameHash{$name}) {
						$add = '';
					} else {
						$add = 1;
						if(defined($nameUpdatedHash{$name})) {
							print $writer ">$nameUpdatedHash{$name}\n";
						} else {
							print $writer ">$name\n";
						}
					}
				} elsif($add) {
					print $writer "$line\n";
				}
			}
			close($reader);
			close($writer);
			system("mv $fastaFile.updated $fastaFile");
		}
		{
			open(my $writer, ">> $fastaFile") or die "Can't write '$fastaFile': $!";
			open(my $reader, "$cdsFastaFile") or die "Can't read '$cdsFastaFile': $!";
			my $add = '';
			my $addCount = 0;
			my $count = 0;
			while(my $line = <$reader>) {
				chomp($line);
				if($line =~ /^>(\S*)/) {
					my $name = $1;
					if($excludeNameHash{$name}) {
						$add = '';
					} else {
						$add = 1;
						$addCount += 1;
						print $writer ">$name\n";
					}
					$count += 1;
				} elsif($add) {
					print $writer "$line\n";
				}
			}
			print STDERR "add $cdsFastaFile $addCount / $count\n";
			close($reader);
			close($writer);
		}
	}
}

if(@fastqFileList) {
	my %orthologyCoverageHash = ();
	foreach my $fastqFile (@fastqFileList) {
		die "'$fastqFile' is not readable.\n" unless(-r $fastqFile);
	}
	my $writer;
	if($pafFile ne '') {
		open($writer, "| gzip > $pafFile") or die "Can't write '$pafFile': $!";
	}
	foreach my $index (0 .. $#fastqFileList) {
		my $fastqFile = $fastqFileList[$index];
		open(my $reader, "minimap2 @optionList $fastaFile $fastqFile |") or die "minimap2 error: $!";
		while(my $line = <$reader>) {
			chomp($line);
			if($line =~ /^@/) {
				if($index == 0) {
					print $writer $line, "\n" if(defined($writer));
				}
				next;
			}
			print $writer $line, "\n" if(defined($writer));
			my %tokenHash = ();
			(@tokenHash{@pafMandatoryColumnList}, my @tagTypeValueList) = split(/\t/, $line, -1);
			$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
			if($stranded eq 'f' || $stranded eq 'forward') {
				next unless($tokenHash{'strand'} eq '+');
			}
			if($stranded eq 'r' || $stranded eq 'reverse') {
				next unless($tokenHash{'strand'} eq '-');
			}
			if($tokenHash{'target_name'} =~ /=>(.*)$/) {
				my $orthology = $1;
				$orthologyCoverageHash{$orthology} = $tokenHash{'match'} / $tokenHash{'target_length'};
			}
		}
		close($reader);
	}
	close($writer) if(defined($writer));
	if($baseAbundanceOrthologies ne '') {
		my @baseAbundanceOrthologyList = split(/,/, $baseAbundanceOrthologies);
		my $baseAbundance = median(@orthologyCoverageHash{@baseAbundanceOrthologyList});
		foreach my $orthology (sort keys %orthologyCoverageHash) {
			$orthologyCoverageHash{$orthology} = $orthologyCoverageHash{$orthology} / $baseAbundance;
		}
	}
	foreach my $orthology (sort keys %orthologyCoverageHash) {
		print join("\t", $orthology, $orthologyCoverageHash{$orthology}), "\n";
	}
}

sub MetaPrism_data {
	if(-r "$MetaPrismDirectory/data/database") {
		chomp(my $database = `cat $MetaPrismDirectory/data/database`);
		return if(-r "$MetaPrismDirectory/data/$database.dmnd");
	}
	system("perl $MetaPrismDirectory/MetaPrism_data.gene.pl");
}

sub lockFile {
	my ($file) = @_;
	while(-r "$file.lock") {
		my ($hostname, $pid) = `cat $file.lock`;
		chomp($hostname);
		chomp($pid);
		chomp(my $stat = `ssh $hostname ps -p $pid --no-headers -o stat`);
		last if($stat eq '' || $stat =~ /Z/);
		sleep(1);
	}
	open(my $writer, "> $file.lock");
	print $writer $hostname, "\n";
	print $writer $$, "\n";
	close($writer);
}

sub unlockFile {
	my ($file) = @_;
	system("rm -f $file.lock");
}

sub getReverseComplementarySequence {
	my ($sequence) = @_;
	($sequence = reverse($sequence)) =~ tr/ACGT/TGCA/;
	return $sequence;
}

sub getCommonTaxonomyIdLineage {
	my ($taxonomyIdLineage1, $taxonomyIdLineage2) = @_;
	return $taxonomyIdLineage1 if($taxonomyIdLineage1 eq $taxonomyIdLineage2);
	my @taxonomyIdList1 = split(/,/, $taxonomyIdLineage1);
	my @taxonomyIdList2 = split(/,/, $taxonomyIdLineage2);
	my $index = 0;
	$index += 1 while($index < scalar(@taxonomyIdList1) && $index < scalar(@taxonomyIdList2) && $taxonomyIdList1[$index] eq $taxonomyIdList2[$index]);
	return join(',', @taxonomyIdList1[0 .. $index - 1]);
}

sub getTaxonomyIdLineage {
	my ($taxonomyId) = @_;
	my @taxonomyIdList = ();
	my $taxon = $db->get_taxon(-taxonid => $taxonomyId);
	while(defined($taxon)) {
		@taxonomyIdList = ($taxon->id, @taxonomyIdList);
		$taxon = $taxon->ancestor;
	}
	return join(',', @taxonomyIdList);
}

sub medianIndexList {
	my ($length) = @_;
	if($length % 2 == 0) {
		return ($length / 2 - 1, $length / 2);
	} else {
		return (($length - 1) / 2);
	}
}

sub median {
	my @tokenList = @_;
	@tokenList = sort {$a <=> $b} @tokenList;
	return mean(@tokenList[medianIndexList(scalar(@tokenList))]);
}

sub mean {
	my @tokenList = @_;
	return sum(@tokenList) / scalar(@tokenList);
}
