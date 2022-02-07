# Iso-seq-tools
Iso-seq data analysis tools, including predict polycistronic transcript, find full length chemic reads,adjest fusion genes
#2022.2.7
#version1
Make sure the following software is installed before running the program:
python versionn 3
conda install -c bioconda gmap
conda install -c bioconda cufflinks
conda install -c bioconda transdecoder
conda install -c bioconda gffcompare
conda install -c bioconda blat

################################################################################
example commandline 1
```
python predict-polycistronic-transcript.py -percent 50.0 -gtf ./test.gtf -ref ./ref.gff3 -genome ./genome.fa --out ./test/

-percent
	the polycistronic transcript coverage CDS of reference gff3 file. This option is used to filter out false positive results with low overlap.
-gtf
	the gtf file generated by ToFu
-ref
	the gff3 file for reference annotation
-genome
	the genome fasta file
-o
	outputdir please note the outpudir must endwith '/'

final output file name 'Final-predict-result.txt'
```
################################################################################
example commandline 2
```
python find-FLC-transcript.py -percent 50 -fasta ./test.fasta -genome ./genome.fa -o ./test/

-percent
	the full length chimeric sequences coverage itself. This option is used to filter out false positive results with low overlap.
-genome
	the genome fasta file
-fasta
	the fasta file for filter full length chimeric sequences. you can use 'gffread(cufflinks)' obtained from gtf file.

-o
	outputdir please note the outpudir must endwith '/'
```

final output file name 'Final-File-Filter.fasta'
################################################################################
example commandline 3
```
python find-adjacent-transcript.py -gtf ./gtf -ref ./ref.gff3 -o ./test/

-gtf
	the gtf file generated by ToFu
-ref
	the gff3/gtf file for reference annotation
-o
	outputdir please note the outpudir must endwith '/'
```
final output file name 'Final-adjacent-transcript-list.txt'
	
